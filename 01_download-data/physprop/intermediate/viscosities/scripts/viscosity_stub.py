"""Viscosity stub to allow multiprocessing to work"""

from openff.evaluator import unit
from openff.evaluator.datasets import PhysicalProperty, PropertyPhase
from openff.evaluator.datasets.thermoml import thermoml_property
from openff.evaluator.plugins import register_default_plugins
from openff.evaluator.datasets.curation.components.thermoml import ImportThermoMLData

register_default_plugins()

VISCOSITY_STR = "Viscosity, Pa*s"

def _process_archive(xml):
    @thermoml_property(VISCOSITY_STR, supported_phases=PropertyPhase.Liquid)
    class Viscosity(PhysicalProperty):
        """Stub representation of a viscosity property"""
        @classmethod
        def default_unit(cls):
            return unit.pascal * unit.seconds
        
    return ImportThermoMLData._process_archive(xml)

