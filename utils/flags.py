from gutter.client import arguments
from gutter.client.models import Switch
from gutter.client.default import gutter
from gutter.client.models import Condition
# from gutter.client.operators.string import EqualsStripIgnoreCase
# from gutter.client.operators.identity import Truthy
from src.phase_parameters import PhaseParameters
from gutter.client.operators.comparable import Equals

class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

class PhaseParametersArguments(arguments.Container):
    COMPATIBLE_TYPE = PhaseParameters
    recompute_peak_positions = arguments.Boolean('recompute_peak_positions')
    preferred_orientation = arguments.Boolean('preferred_orientation')


pref_orient_switch = Switch('pref_or', state=Switch.states.SELECTIVE)
gutter.register(pref_orient_switch)
pref_orient_switch.conditions = [
    Condition(argument=PhaseParametersArguments,
              attribute='preferred_orientation',
              operator=Equals(value=True)),
]
gutter.update(pref_orient_switch)

recompute_peak_positions_switch = Switch('rpp',
    state=Switch.states.SELECTIVE)
gutter.register(recompute_peak_positions_switch)
recompute_peak_positions_switch.conditions = [
    Condition(argument=PhaseParametersArguments,
              attribute='recompute_peak_positions',
              operator=Equals(value=True)),
]
gutter.update(recompute_peak_positions_switch)
