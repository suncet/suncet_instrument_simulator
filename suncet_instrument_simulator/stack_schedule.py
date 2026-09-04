"""
SHDR observation timing: map onboard integration stacks to MHD radiance map indices.
"""
from dataclasses import dataclass
import math

import astropy.units as u
from astropy.time import Time, TimeDelta
import numpy as np
import sunpy.map


@dataclass(frozen=True)
class MapContribution:
    model_index: int
    weight: float


@dataclass(frozen=True)
class Observation:
    output_index: int
    start_seconds: float


@dataclass
class ExposureStackSchedule:
    start_seconds: float
    short_members: list
    long_members: list

    @property
    def unique_model_indices(self):
        indices = set()
        for member in self.short_members + self.long_members:
            for contribution in member:
                indices.add(contribution.model_index)
        return sorted(indices)


def _to_seconds(quantity):
  return float(quantity.to_value(u.s)) if hasattr(quantity, 'to_value') else float(quantity)


def integration_window_seconds(exposure_time, num_integrations):
  return _to_seconds(exposure_time) * num_integrations


def observation_window_seconds(config):
  short_window = integration_window_seconds(
    config.exposure_time_short, config.num_short_exposures_to_stack)
  long_window = integration_window_seconds(
    config.exposure_time_long, config.num_long_exposures_to_stack)
  return max(short_window, long_window)


def output_stride_model_steps(config):
  model_dt = _to_seconds(config.model_timestep)
  window = observation_window_seconds(config)
  return window / model_dt


def build_observation_sequence(config):
  """Build instrument-cadence output starts within the configured model bounds."""
  # The third timesteps_to_process value selects maps during radiance-map
  # generation; it is not an instrument output stride.
  first_model_index, last_model_index, _radiance_map_generation_step = (
    config.timesteps_to_process)
  model_dt = _to_seconds(config.model_timestep)
  cadence = observation_window_seconds(config)
  if cadence <= 0:
    raise ValueError('Observation cadence must be positive.')

  first_start = first_model_index * model_dt
  last_start = last_model_index * model_dt
  num_intervals = int(math.floor((last_start - first_start) / cadence + 1e-12))
  return [
    Observation(output_index=index, start_seconds=first_start + index * cadence)
    for index in range(num_intervals + 1)
  ]


def model_indices_to_load_at_time(start_seconds, config):
  model_dt = _to_seconds(config.model_timestep)
  window = observation_window_seconds(config)
  end_time = start_seconds + window
  first_index = int(math.floor(start_seconds / model_dt))
  last_index = int(math.floor((end_time - 1e-12) / model_dt))
  return list(range(first_index, last_index + 1))


def model_indices_to_load(t0, config):
  """Backward-compatible wrapper for a start expressed as a model index."""
  return model_indices_to_load_at_time(t0 * _to_seconds(config.model_timestep), config)


def _overlap_contributions_at_time(
    observation_start_seconds, integration_index, exposure_time, model_timestep):
  exposure_s = _to_seconds(exposure_time)
  model_dt = _to_seconds(model_timestep)
  t_start = observation_start_seconds + integration_index * exposure_s
  t_end = t_start + exposure_s

  contributions = []
  n_start = int(math.floor(t_start / model_dt))
  n_end = int(math.floor((t_end - 1e-12) / model_dt))

  for n in range(n_start, n_end + 1):
    bin_start = n * model_dt
    bin_end = (n + 1) * model_dt
    overlap_start = max(t_start, bin_start)
    overlap_end = min(t_end, bin_end)
    if overlap_end > overlap_start:
      weight = (overlap_end - overlap_start) / exposure_s
      contributions.append(MapContribution(model_index=n, weight=weight))

  if not contributions:
    contributions.append(MapContribution(model_index=n_start, weight=1.0))

  return contributions


def _overlap_contributions(t0, integration_index, exposure_time, model_timestep):
  """Backward-compatible wrapper for a start expressed as a model index."""
  model_dt = _to_seconds(model_timestep)
  return _overlap_contributions_at_time(
    t0 * model_dt, integration_index, exposure_time, model_timestep)


def _pad_contributions(contributions, available_indices):
  if not available_indices:
    raise ValueError('No radiance maps available to pad stack schedule.')

  last_in_window = max(available_indices)
  padded = []
  for contribution in contributions:
    if contribution.model_index in available_indices:
      padded.append(contribution)
    else:
      padded.append(MapContribution(model_index=last_in_window, weight=contribution.weight))
  return padded


def build_exposure_stack_members(start_seconds, exposure_time, num_integrations, model_timestep,
                                 available_indices):
  members = []
  for integration_index in range(num_integrations):
    contributions = _overlap_contributions_at_time(
      start_seconds, integration_index, exposure_time, model_timestep)
    members.append(_pad_contributions(contributions, available_indices))
  return members


def build_stack_schedule_at_time(start_seconds, config, available_indices=None):
  if available_indices is None:
    available_indices = model_indices_to_load_at_time(start_seconds, config)
  available_set = set(available_indices)

  short_members = build_exposure_stack_members(
    start_seconds,
    config.exposure_time_short,
    config.num_short_exposures_to_stack,
    config.model_timestep,
    available_set,
  )
  long_members = build_exposure_stack_members(
    start_seconds,
    config.exposure_time_long,
    config.num_long_exposures_to_stack,
    config.model_timestep,
    available_set,
  )

  return ExposureStackSchedule(
    start_seconds=start_seconds,
    short_members=short_members,
    long_members=long_members,
  )


def build_stack_schedule(t0, config, available_indices=None):
  """Backward-compatible wrapper for a start expressed as a model index."""
  start_seconds = t0 * _to_seconds(config.model_timestep)
  return build_stack_schedule_at_time(start_seconds, config, available_indices)


def combine_radiance_for_member(contributions, radiance_by_model_index, *,
                                start_seconds, model_timestep):
  normalized = _normalize_contributions(contributions)
  reference_index = normalized[0].model_index
  reference_map = _get_reference_wavelength_map(radiance_by_model_index, reference_index)

  combined_by_wavelength = {}
  for wavelength, reference in reference_map.items():
    combined_data = np.zeros_like(reference.data, dtype=np.float64)
    for contribution in normalized:
      source_map = radiance_by_model_index[contribution.model_index][wavelength]
      combined_data += contribution.weight * source_map.data
    # The model timestamp already includes reference_index * model_timestep.
    # Add only the integration's offset from that frame. This also timestamps
    # integrations between model frames and integrations padded with the final
    # frame correctly, without advancing the model epoch a second time.
    timestamp = reference.meta.get('DATE-OBS')
    if timestamp in (None, '', 'N/A'):
      raise ValueError('Radiance map does not contain a valid DATE-OBS.')
    time_scale = str(reference.meta.get('TIMESYS', 'UTC')).strip().lower()
    reference_time = Time(timestamp, format='isot', scale=time_scale, precision=6).utc
    offset_seconds = start_seconds - reference_index * _to_seconds(model_timestep)
    integration_start = reference_time + TimeDelta(offset_seconds, format='sec')
    metadata = reference.meta.copy()
    metadata['DATE-OBS'] = integration_start.isot
    metadata['TIMESYS'] = 'UTC'
    combined_by_wavelength[wavelength] = sunpy.map.Map(combined_data, metadata)

  return combined_by_wavelength


def build_radiance_by_stack_member(stack_members, radiance_by_model_index, *,
                                   start_seconds, exposure_time, model_timestep):
  """Combine each integration's scene and timestamp its actual start time."""
  radiance_by_member = {}
  for member_index, contributions in enumerate(stack_members):
    radiance_by_member[member_index] = combine_radiance_for_member(
      contributions, radiance_by_model_index,
      start_seconds=start_seconds + member_index * _to_seconds(exposure_time),
      model_timestep=model_timestep)
  return radiance_by_member


def _normalize_contributions(contributions):
  total_weight = sum(contribution.weight for contribution in contributions)
  if total_weight <= 0:
    raise ValueError('Integration contributions must have positive total weight.')
  return [
    MapContribution(contribution.model_index, contribution.weight / total_weight)
    for contribution in contributions
  ]


def _get_reference_wavelength_map(radiance_by_model_index, model_index):
  if model_index not in radiance_by_model_index:
    raise KeyError('Radiance map for model index {} is not loaded.'.format(model_index))
  return radiance_by_model_index[model_index]
