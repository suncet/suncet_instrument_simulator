import math
from pathlib import Path

import astropy.units as u
import pytest

from suncet_instrument_simulator import config_parser, stack_schedule


CONFIG_DIR = Path(__file__).resolve().parents[1] / 'suncet_instrument_simulator' / 'config_files'


class _ConfigStub:
    def __init__(self):
        self.exposure_time_short = 0.0035 * u.s
        self.exposure_time_long = 15 * u.s
        self.model_timestep = 10 * u.s
        self.num_short_exposures_to_stack = 9
        self.num_long_exposures_to_stack = 4
        self.filter_out_particle_hits = True
        self.timesteps_to_process = [0, 360, 1]


def test_short_stack_uses_single_model_index():
    config = _ConfigStub()
    schedule = stack_schedule.build_stack_schedule(4, config, available_indices=[4, 5, 6, 7, 8, 9])

    assert len(schedule.short_members) == 9
    for member in schedule.short_members:
        assert len(member) == 1
        assert member[0].model_index == 4
        assert member[0].weight == pytest.approx(1.0)


def test_long_first_integration_spans_two_model_maps():
    config = _ConfigStub()
    schedule = stack_schedule.build_stack_schedule(0, config, available_indices=list(range(8)))

    first_member = schedule.long_members[0]
    assert len(first_member) == 2
    assert first_member[0].model_index == 0
    assert first_member[0].weight == pytest.approx(10.0 / 15.0)
    assert first_member[1].model_index == 1
    assert first_member[1].weight == pytest.approx(5.0 / 15.0)


def test_long_second_integration_start_spans_fractional_model_boundary():
    config = _ConfigStub()
    config.num_short_exposures_to_stack = 1
    config.num_long_exposures_to_stack = 1
    schedule = stack_schedule.build_stack_schedule_at_time(
        15.0, config, available_indices=[1, 2])

    member = schedule.long_members[0]
    assert [contribution.model_index for contribution in member] == [1, 2]
    assert [contribution.weight for contribution in member] == pytest.approx(
        [5.0 / 15.0, 10.0 / 15.0]
    )


def test_arbitrary_configured_timings_drive_cadence_and_overlap_weights():
    config = _ConfigStub()
    config.model_timestep = 7 * u.s
    config.exposure_time_short = 2 * u.s
    config.num_short_exposures_to_stack = 4
    config.exposure_time_long = 11 * u.s
    config.num_long_exposures_to_stack = 3
    config.timesteps_to_process = [2, 20, 1]

    observations = stack_schedule.build_observation_sequence(config)
    schedule = stack_schedule.build_stack_schedule_at_time(47.0, config)

    assert stack_schedule.output_stride_model_steps(config) == pytest.approx(33.0 / 7.0)
    assert [observation.start_seconds for observation in observations] == [14, 47, 80, 113]
    assert [contribution.model_index for contribution in schedule.long_members[0]] == [6, 7, 8]
    assert [contribution.weight for contribution in schedule.long_members[0]] == pytest.approx(
        [2.0 / 11.0, 7.0 / 11.0, 2.0 / 11.0]
    )


def test_output_stride_matches_observation_window():
    config = _ConfigStub()
    assert stack_schedule.output_stride_model_steps(config) == 6


def test_model_indices_to_load_for_sixty_second_window():
    config = _ConfigStub()
    assert stack_schedule.model_indices_to_load(0, config) == [0, 1, 2, 3, 4, 5]


def test_filter_false_forces_single_integration(config_filename):
    filename = config_filename
    config = config_parser.Config(filename)
    assert config.num_short_exposures_to_stack == 1
    assert config.num_long_exposures_to_stack == 1
    assert config.output_stride_model_steps == pytest.approx(1.5)


def test_filter_false_sequence_uses_instrument_cadence(config_filename):
    config = config_parser.Config(config_filename)
    observations = stack_schedule.build_observation_sequence(config)

    assert len(observations) == 241
    assert observations[0] == stack_schedule.Observation(0, 0.0)
    assert observations[1] == stack_schedule.Observation(1, 15.0)
    assert observations[-1] == stack_schedule.Observation(240, 3600.0)


def test_filter_false_pads_terminal_long_exposure_with_last_available_map(config_filename):
    config = config_parser.Config(config_filename)
    schedule = stack_schedule.build_stack_schedule(360, config, available_indices=[360])

    assert schedule.unique_model_indices == [360]
    assert len(schedule.long_members) == 1
    assert [contribution.model_index for contribution in schedule.long_members[0]] == [360, 360]
    assert [contribution.weight for contribution in schedule.long_members[0]] == pytest.approx(
        [10.0 / 15.0, 5.0 / 15.0]
    )


def test_filter_true_requires_multi_member_stack():
    filename = CONFIG_DIR / 'config_default.ini'
    config = config_parser.Config(filename)
    assert config.num_short_exposures_to_stack >= 2
    assert config.num_long_exposures_to_stack >= 2
    assert config.output_stride_model_steps == math.ceil(60 / 10)


def test_filter_true_sequence_uses_collapsed_stack_window():
    filename = CONFIG_DIR / 'config_default.ini'
    config = config_parser.Config(filename)
    observations = stack_schedule.build_observation_sequence(config)

    assert len(observations) == 61
    assert observations[0].start_seconds == 0
    assert observations[1].start_seconds == 60
    assert observations[-1].start_seconds == 3600


@pytest.fixture
def config_filename():
    return CONFIG_DIR / 'config_mars_telecom_network.ini'
