# src/relativity.py
"""
General Relativity corrections for celestial mechanics

This module provides functions for calculating relativistic corrections
to classical orbital mechanics, including pericenter precession and time delays.
"""
import numpy as np

# Physical constants
G = 6.67430e-11  # Gravitational constant [m^3 kg^-1 s^-2]
c = 299792458.0  # Speed of light [m/s]


def pericenter_precession(a, e, M_central):
    """
    Calculate pericenter precession angle per orbit (Δω) in Schwarzschild metric
    
    This is the additional precession due to general relativistic effects,
    beyond classical perturbations.
    
    Parameters:
    -----------
    a : float
        Semi-major axis [m]
    e : float
        Eccentricity
    M_central : float
        Mass of central body [kg]
    
    Returns:
    --------
    float
        Precession angle per orbit [rad]
    """
    return 6 * np.pi * G * M_central / (c**2 * a * (1 - e**2))


def pericenter_precession_per_year(a, e, M_central, orbital_period):
    """
    Calculate pericenter precession rate per year
    
    Parameters:
    -----------
    a : float
        Semi-major axis [m]
    e : float
        Eccentricity
    M_central : float
        Mass of central body [kg]
    orbital_period : float
        Orbital period [s]
    
    Returns:
    --------
    float
        Precession rate [arcseconds/year]
    """
    precession_per_orbit = pericenter_precession(a, e, M_central)
    seconds_per_year = 365.25 * 24 * 3600
    orbits_per_year = seconds_per_year / orbital_period
    
    # Convert radians to arcseconds
    arcsec_per_rad = 206264.806
    
    return precession_per_orbit * orbits_per_year * arcsec_per_rad


def shapiro_delay(r1, r2, r_closest, M_central):
    """
    Calculate Shapiro time delay for light passing near a massive body
    
    Parameters:
    -----------
    r1 : float
        Distance from central body to transmission point [m]
    r2 : float
        Distance from central body to reception point [m]
    r_closest : float
        Closest approach distance to central body [m]
    M_central : float
        Mass of central body [kg]
    
    Returns:
    --------
    float
        Time delay [s]
    """
    rs = schwarzschild_radius(M_central)
    
    # Approximation for small deflection angles
    # More accurate formula would require numerical integration
    return rs / c * np.log((r1 + r2 + np.sqrt((r1 + r2)**2 - 4*r_closest**2)) / 
                          (2 * r_closest))


def schwarzschild_radius(M):
    """
    Calculate Schwarzschild radius
    
    Parameters:
    -----------
    M : float
        Mass [kg]
    
    Returns:
    --------
    float
        Schwarzschild radius [m]
    """
    return 2 * G * M / c**2


def gravitational_redshift(r, M_central):
    """
    Calculate gravitational redshift factor
    
    Parameters:
    -----------
    r : float
        Distance from central body [m]
    M_central : float
        Mass of central body [kg]
    
    Returns:
    --------
    float
        Redshift factor (1 + z)
    """
    rs = schwarzschild_radius(M_central)
    return 1 / np.sqrt(1 - rs / r)


def light_deflection_angle(impact_parameter, M_central):
    """
    Calculate light deflection angle in Schwarzschild spacetime
    
    Parameters:
    -----------
    impact_parameter : float
        Impact parameter (closest approach distance) [m]
    M_central : float
        Mass of central body [kg]
    
    Returns:
    --------
    float
        Deflection angle [rad]
    """
    rs = schwarzschild_radius(M_central)
    return 2 * rs / impact_parameter


def time_dilation_factor(r, M_central):
    """
    Calculate gravitational time dilation factor
    
    Parameters:
    -----------
    r : float
        Distance from central body [m]
    M_central : float
        Mass of central body [kg]
    
    Returns:
    --------
    float
        Time dilation factor (proper time / coordinate time)
    """
    rs = schwarzschild_radius(M_central)
    return np.sqrt(1 - rs / r)


def effective_potential(r, l, M_central, test_mass=1.0):
    """
    Calculate effective potential in Schwarzschild metric
    
    Parameters:
    -----------
    r : float or array-like
        Radial distance [m]
    l : float
        Angular momentum per unit mass [m^2/s]
    M_central : float
        Mass of central body [kg]
    test_mass : float, optional
        Test particle mass [kg]
    
    Returns:
    --------
    float or array-like
        Effective potential [J] (if test_mass=1, then [J/kg])
    """
    rs = schwarzschild_radius(M_central)
    
    # Classical terms
    centrifugal = l**2 / (2 * test_mass * r**2)
    gravitational = -G * M_central * test_mass / r
    
    # Relativistic correction
    relativistic = -rs * c**2 * l**2 / (2 * r**3)
    
    return centrifugal + gravitational + relativistic


def circular_orbit_radius_gr(l, M_central):
    """
    Find radius of circular orbit in Schwarzschild metric
    
    Parameters:
    -----------
    l : float
        Angular momentum per unit mass [m^2/s]
    M_central : float
        Mass of central body [kg]
    
    Returns:
    --------
    float
        Circular orbit radius [m]
    """
    rs = schwarzschild_radius(M_central)
    mu = G * M_central
    
    # Solve cubic equation for circular orbit
    # This is an approximation; exact solution requires numerical methods
    r_classical = l**2 / mu
    
    # First-order correction
    correction = 3 * rs / 2
    
    return r_classical + correction


def frame_dragging_precession(a, J, M_central, r):
    """
    Calculate frame-dragging (Lense-Thirring) precession rate
    
    This is an additional effect in the Kerr metric for rotating central bodies.
    
    Parameters:
    -----------
    a : float
        Kerr parameter (J*c/(M*G)) [m]
    J : float
        Angular momentum of central body [kg*m^2/s]
    M_central : float
        Mass of central body [kg]
    r : float
        Orbital radius [m]
    
    Returns:
    --------
    float
        Frame-dragging precession rate [rad/s]
    """
    rs = schwarzschild_radius(M_central)
    
    # Approximate formula for weak field
    return 2 * G * J / (c**2 * r**3)


def post_newtonian_correction(r, v, M_central):
    """
    Calculate first-order post-Newtonian correction to acceleration
    
    Parameters:
    -----------
    r : float
        Distance from central body [m]
    v : float
        Orbital velocity [m/s]
    M_central : float
        Mass of central body [kg]
    
    Returns:
    --------
    float
        Acceleration correction factor
    """
    rs = schwarzschild_radius(M_central)
    beta = v / c
    
    # First-order post-Newtonian correction
    pn1 = rs / r * (1 + 3*beta**2)
    
    return 1 + pn1


def mercury_precession_test(a_mercury=5.791e10, e_mercury=0.2056):
    """
    Calculate Mercury's perihelion precession as a test case
    
    Parameters:
    -----------
    a_mercury : float, optional
        Mercury's semi-major axis [m]
    e_mercury : float, optional
        Mercury's eccentricity
    
    Returns:
    --------
    dict
        Dictionary with calculated and observed precession values
    """
    M_sun = 1.9885e30  # Solar mass [kg]
    orbital_period_mercury = 87.97 * 24 * 3600  # seconds
    
    # Calculate GR precession
    gr_precession = pericenter_precession_per_year(
        a_mercury, e_mercury, M_sun, orbital_period_mercury
    )
    
    # Observed total precession is about 574.10"/century
    # GR contributes about 42.98"/century
    observed_gr_precession = 42.98  # arcsec/century
    observed_gr_per_year = observed_gr_precession / 100.0  # arcsec/year
    
    return {
        'calculated_gr_precession_per_year': gr_precession,
        'observed_gr_precession_per_year': observed_gr_per_year,
        'relative_error': abs(gr_precession - observed_gr_per_year) / observed_gr_per_year
    }