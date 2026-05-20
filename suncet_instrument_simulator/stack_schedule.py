"""
SHDR observation timing: map onboard integration stacks to MHD radiance map indices.
"""
from dataclasses import dataclass
import math

import astropy.units as u
import numpy as np
import sunpy.map


@dataclass(frozen=True)
class MapContribution:
    model_index: int
    weight: float


@dataclass
class ExposureStackSchedule:
    t0: int
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
  if not config.filter_out_particle_hits:
    return 1
  model_dt = _to_seconds(config.model_timestep)
  window = observation_window_seconds(config)
  return max(1, int(math.ceil(window / model_dt)))


def model_indices_to_load(t0, config):
  model_dt = _to_seconds(config.model_timestep)
  window = observation_window_seconds(config)
  t0_seconds = t0 * model_dt
  end_time = t0_seconds + window
  indices = []
  n = t0
  while n * model_dt < end_time:
    indices.append(n)
    n += 1
  return indices


def _overlap_contributions(t0, integration_index, exposure_time, model_timestep):
  exposure_s = _to_seconds(exposure_time)
  model_dt = _to_seconds(model_timestep)
  t0_seconds = t0 * model_dt
  t_start = t0_seconds + integration_index * exposure_s
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
    contributions.append(MapContribution(model_index=max(t0, n_start), weight=1.0))

  return contributions


def _pad_contributions(contributions, available_indices, t0):
  if not available_indices:
    raise ValueError('No radiance maps available to pad stack schedule.')

  last_in_window = max(i for i in available_indices if i >= t0)
  padded = []
  for contribution in contributions:
    if contribution.model_index in available_indices:
      padded.append(contribution)
    else:
      padded.append(MapContribution(model_index=last_in_window, weight=contribution.weight))
  return padded


def build_exposure_stack_members(t0, exposure_time, num_integrations, model_timestep,
                                 available_indices):
  members = []
  for integration_index in range(num_integrations):
    contributions = _overlap_contributions(
      t0, integration_index, exposure_time, model_timestep)
    members.append(_pad_contributions(contributions, available_indices, t0))
  return members


def build_stack_schedule(t0, config, available_indices=None):
  if available_indices is None:
    available_indices = model_indices_to_load(t0, config)
  available_set = set(available_indices)

  short_members = build_exposure_stack_members(
    t0,
    config.exposure_time_short,
    config.num_short_exposures_to_stack,
    config.model_timestep,
    available_set,
  )
  long_members = build_exposure_stack_members(
    t0,
    config.exposure_time_long,
    config.num_long_exposures_to_stack,
    config.model_timestep,
    available_set,
  )

  return ExposureStackSchedule(t0=t0, short_members=short_members, long_members=long_members)


def combine_radiance_for_member(contributions, radiance_by_model_index):
  normalized = _normalize_contributions(contributions)
  reference_index = normalized[0].model_index
  reference_map = _get_reference_wavelength_map(radiance_by_model_index, reference_index)

  combined_by_wavelength = {}
  for wavelength, reference in reference_map.items():
    combined_data = np.zeros_like(reference.data, dtype=np.float64)
    for contribution in normalized:
      source_map = radiance_by_model_index[contribution.model_index][wavelength]
      combined_data += contribution.weight * source_map.data
    combined_by_wavelength[wavelength] = sunpy.map.Map(combined_data, reference.meta)

  return combined_by_wavelength


def build_radiance_by_stack_member(stack_members, radiance_by_model_index):
  radiance_by_member = {}
  for member_index, contributions in enumerate(stack_members):
    radiance_by_member[member_index] = combine_radiance_for_member(
      contributions, radiance_by_model_index)
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
