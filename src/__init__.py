"""
AstroLibrary - Celestial Mechanics Python Library

A comprehensive physics simulation library for celestial mechanics,
combining classical Newtonian mechanics with General Relativity corrections.

Modules:
--------
- classical: Classical celestial mechanics and Kepler problem
- relativity: General relativity corrections 
- utils: Utility functions and constants
"""

from . import classical
from . import relativity
from . import utils

__version__ = "1.0.0"
__author__ = "AstroLibrary Contributors"
__description__ = "Celestial Mechanics Python Library"

# Re-export commonly used functions
from .classical import (
    solve_kepler,
    orbital_elements_to_state,
    orbital_period,
    vis_viva_equation
)

from .relativity import (
    pericenter_precession,
    shapiro_delay,
    mercury_precession_test
)

from .utils import (
    deg_to_rad,
    rad_to_deg,
    au_to_meters,
    solar_masses_to_kg,
    G, c, AU, solar_mass
)

__all__ = [
    'classical',
    'relativity', 
    'utils',
    'solve_kepler',
    'orbital_elements_to_state',
    'orbital_period',
    'vis_viva_equation',
    'pericenter_precession',
    'shapiro_delay',
    'mercury_precession_test',
    'deg_to_rad',
    'rad_to_deg',
    'au_to_meters',
    'solar_masses_to_kg',
    'G', 'c', 'AU', 'solar_mass'
]