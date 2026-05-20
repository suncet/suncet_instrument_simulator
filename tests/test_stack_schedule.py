import math

import astropy.units as u
import pytest

from suncet_instrument_simulator import config_parser, stack_schedule


class _ConfigStub:
    def __init__(self):
        self.exposure_time_short = 0.0035 * u.s
        self.exposure_time_long = 15 * u.s
        self.model_timestep = 10 * u.s
        self.num_short_exposures_to_stack = 9
        self.num_long_exposures_to_stack = 4
        self.filter_out_particle_hits = True


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
    assert config.output_stride_model_steps == 1


def test_filter_true_requires_multi_member_stack():
    filename = (
        '/Users/masonjp2/Dropbox/suncet_dropbox/9000 Processing/code/suncet_instrument_simulator/'
        'suncet_instrument_simulator/config_files/config_default.ini'
    )
    config = config_parser.Config(filename)
    assert config.num_short_exposures_to_stack >= 2
    assert config.num_long_exposures_to_stack >= 2
    assert config.output_stride_model_steps == math.ceil(60 / 10)


@pytest.fixture
def config_filename():
    return (
        '/Users/masonjp2/Dropbox/suncet_dropbox/9000 Processing/code/suncet_instrument_simulator/'
        'suncet_instrument_simulator/config_files/config_mars_telecom_network.ini'
    )
