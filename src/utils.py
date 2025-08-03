# src/utils.py
"""
Utility functions for celestial mechanics calculations

This module provides common utility functions, unit conversions,
and astronomical constants used throughout the library.
"""
import numpy as np

# Physical constants
G = 6.67430e-11  # Gravitational constant [m^3 kg^-1 s^-2]
c = 299792458.0  # Speed of light [m/s]
h = 6.62607015e-34  # Planck constant [J⋅s]
k_B = 1.380649e-23  # Boltzmann constant [J/K]

# Astronomical constants
AU = 1.495978707e11  # Astronomical unit [m]
pc = 3.0857e16  # Parsec [m]
ly = 9.4607e15  # Light year [m]
solar_mass = 1.9885e30  # Solar mass [kg]
earth_mass = 5.9722e24  # Earth mass [kg]
solar_radius = 6.96e8  # Solar radius [m]
earth_radius = 6.371e6  # Earth radius [m]

# Time conversions
seconds_per_day = 86400.0
seconds_per_year = 365.25 * seconds_per_day
julian_year = 365.25 * seconds_per_day


def deg_to_rad(degrees):
    """Convert degrees to radians"""
    return degrees * np.pi / 180.0


def rad_to_deg(radians):
    """Convert radians to degrees"""
    return radians * 180.0 / np.pi


def arcsec_to_rad(arcseconds):
    """Convert arcseconds to radians"""
    return arcseconds * np.pi / (180.0 * 3600.0)


def rad_to_arcsec(radians):
    """Convert radians to arcseconds"""
    return radians * 180.0 * 3600.0 / np.pi


def au_to_meters(au):
    """Convert astronomical units to meters"""
    return au * AU


def meters_to_au(meters):
    """Convert meters to astronomical units"""
    return meters / AU


def solar_masses_to_kg(solar_masses):
    """Convert solar masses to kilograms"""
    return solar_masses * solar_mass


def kg_to_solar_masses(kg):
    """Convert kilograms to solar masses"""
    return kg / solar_mass


def years_to_seconds(years):
    """Convert years to seconds"""
    return years * seconds_per_year


def seconds_to_years(seconds):
    """Convert seconds to years"""
    return seconds / seconds_per_year


def days_to_seconds(days):
    """Convert days to seconds"""
    return days * seconds_per_day


def seconds_to_days(seconds):
    """Convert seconds to days"""
    return seconds / seconds_per_day


def km_per_s_to_m_per_s(km_per_s):
    """Convert km/s to m/s"""
    return km_per_s * 1000.0


def m_per_s_to_km_per_s(m_per_s):
    """Convert m/s to km/s"""
    return m_per_s / 1000.0


def standard_gravitational_parameter(mass_kg):
    """
    Calculate standard gravitational parameter μ = GM
    
    Parameters:
    -----------
    mass_kg : float
        Mass in kilograms
    
    Returns:
    --------
    float
        Standard gravitational parameter [m^3/s^2]
    """
    return G * mass_kg


def orbital_velocity_circular(r, mu):
    """
    Calculate circular orbital velocity
    
    Parameters:
    -----------
    r : float
        Orbital radius [m]
    mu : float
        Standard gravitational parameter [m^3/s^2]
    
    Returns:
    --------
    float
        Circular orbital velocity [m/s]
    """
    return np.sqrt(mu / r)


def escape_velocity(r, mu):
    """
    Calculate escape velocity
    
    Parameters:
    -----------
    r : float
        Distance from center [m]
    mu : float
        Standard gravitational parameter [m^3/s^2]
    
    Returns:
    --------
    float
        Escape velocity [m/s]
    """
    return np.sqrt(2 * mu / r)


def sphere_of_influence(m1, m2, a):
    """
    Calculate sphere of influence radius
    
    Parameters:
    -----------
    m1 : float
        Mass of primary body [kg]
    m2 : float
        Mass of secondary body [kg]
    a : float
        Semi-major axis of secondary's orbit around primary [m]
    
    Returns:
    --------
    float
        Sphere of influence radius [m]
    """
    return a * (m2 / m1)**(2/5)


def hill_sphere(m1, m2, a):
    """
    Calculate Hill sphere radius
    
    Parameters:
    -----------
    m1 : float
        Mass of primary body [kg]
    m2 : float
        Mass of secondary body [kg]
    a : float
        Semi-major axis of secondary's orbit around primary [m]
    
    Returns:
    --------
    float
        Hill sphere radius [m]
    """
    return a * (m2 / (3 * m1))**(1/3)


def roche_limit_rigid(rho_primary, rho_satellite):
    """
    Calculate Roche limit for a rigid satellite
    
    Parameters:
    -----------
    rho_primary : float
        Density of primary body [kg/m^3]
    rho_satellite : float
        Density of satellite [kg/m^3]
    
    Returns:
    --------
    float
        Roche limit as multiple of primary radius
    """
    return 2.44 * (rho_primary / rho_satellite)**(1/3)


def roche_limit_fluid(rho_primary, rho_satellite):
    """
    Calculate Roche limit for a fluid satellite
    
    Parameters:
    -----------
    rho_primary : float
        Density of primary body [kg/m^3]
    rho_satellite : float
        Density of satellite [kg/m^3]
    
    Returns:
    --------
    float
        Roche limit as multiple of primary radius
    """
    return 2.44 * (rho_primary / rho_satellite)**(1/3)


def synodic_period(P1, P2):
    """
    Calculate synodic period of two orbiting bodies
    
    Parameters:
    -----------
    P1 : float
        Orbital period of first body [s]
    P2 : float
        Orbital period of second body [s]
    
    Returns:
    --------
    float
        Synodic period [s]
    """
    return abs(1 / (1/P1 - 1/P2))


def angular_separation(ra1, dec1, ra2, dec2):
    """
    Calculate angular separation between two celestial coordinates
    
    Parameters:
    -----------
    ra1, ra2 : float
        Right ascension [rad]
    dec1, dec2 : float
        Declination [rad]
    
    Returns:
    --------
    float
        Angular separation [rad]
    """
    cos_sep = (np.sin(dec1) * np.sin(dec2) + 
               np.cos(dec1) * np.cos(dec2) * np.cos(ra2 - ra1))
    return np.arccos(np.clip(cos_sep, -1, 1))


def magnitude_to_flux(magnitude, zero_point=0.0):
    """
    Convert magnitude to flux
    
    Parameters:
    -----------
    magnitude : float
        Apparent magnitude
    zero_point : float, optional
        Zero point magnitude
    
    Returns:
    --------
    float
        Relative flux
    """
    return 10**(-0.4 * (magnitude - zero_point))


def flux_to_magnitude(flux, zero_point=0.0):
    """
    Convert flux to magnitude
    
    Parameters:
    -----------
    flux : float
        Relative flux
    zero_point : float, optional
        Zero point magnitude
    
    Returns:
    --------
    float
        Apparent magnitude
    """
    return -2.5 * np.log10(flux) + zero_point


def luminosity_distance(absolute_magnitude, apparent_magnitude):
    """
    Calculate luminosity distance from magnitudes
    
    Parameters:
    -----------
    absolute_magnitude : float
        Absolute magnitude
    apparent_magnitude : float
        Apparent magnitude
    
    Returns:
    --------
    float
        Luminosity distance [pc]
    """
    distance_modulus = apparent_magnitude - absolute_magnitude
    return 10**(distance_modulus / 5 - 1)  # in parsecs


def parallax_to_distance(parallax_arcsec):
    """
    Convert parallax to distance
    
    Parameters:
    -----------
    parallax_arcsec : float
        Parallax [arcseconds]
    
    Returns:
    --------
    float
        Distance [pc]
    """
    return 1.0 / parallax_arcsec


def blackbody_temperature(wavelength_peak):
    """
    Calculate blackbody temperature from peak wavelength (Wien's law)
    
    Parameters:
    -----------
    wavelength_peak : float
        Peak wavelength [m]
    
    Returns:
    --------
    float
        Temperature [K]
    """
    wien_constant = 2.897771955e-3  # m⋅K
    return wien_constant / wavelength_peak


def blackbody_flux(temperature, wavelength):
    """
    Calculate blackbody flux at given wavelength (Planck function)
    
    Parameters:
    -----------
    temperature : float
        Temperature [K]
    wavelength : float
        Wavelength [m]
    
    Returns:
    --------
    float
        Spectral radiance [W⋅sr^-1⋅m^-3]
    """
    return (2 * h * c**2 / wavelength**5) / (np.exp(h * c / (wavelength * k_B * temperature)) - 1)


def create_time_array(start_time, end_time, num_points):
    """
    Create evenly spaced time array
    
    Parameters:
    -----------
    start_time : float
        Start time [s]
    end_time : float
        End time [s]
    num_points : int
        Number of time points
    
    Returns:
    --------
    numpy.ndarray
        Time array [s]
    """
    return np.linspace(start_time, end_time, num_points)


def rotation_matrix_z(angle):
    """
    Create rotation matrix around z-axis
    
    Parameters:
    -----------
    angle : float
        Rotation angle [rad]
    
    Returns:
    --------
    numpy.ndarray
        3x3 rotation matrix
    """
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    return np.array([
        [cos_a, -sin_a, 0],
        [sin_a,  cos_a, 0],
        [0,      0,     1]
    ])


def coordinate_transform_2d(x, y, angle):
    """
    Transform 2D coordinates by rotation
    
    Parameters:
    -----------
    x, y : float or array-like
        Input coordinates
    angle : float
        Rotation angle [rad]
    
    Returns:
    --------
    tuple
        Transformed (x', y') coordinates
    """
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    x_new = x * cos_a - y * sin_a
    y_new = x * sin_a + y * cos_a
    return x_new, y_new