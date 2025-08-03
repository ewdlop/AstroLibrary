# src/classical.py
"""
Classical celestial mechanics module

This module provides functions for classical two-body problem calculations,
Kepler's equations, and orbital element conversions.
"""
import numpy as np
from scipy.optimize import newton

# Gravitational constant [m^3 kg^-1 s^-2]
G = 6.67430e-11


def kepler_equation(E, M, e):
    """
    Kepler's equation: E - e*sin(E) = M
    
    Parameters:
    -----------
    E : float
        Eccentric anomaly [rad]
    M : float
        Mean anomaly [rad]
    e : float
        Eccentricity
    
    Returns:
    --------
    float
        The residual of Kepler's equation
    """
    return E - e * np.sin(E) - M


def solve_kepler(M, e, tolerance=1e-12):
    """
    Solve Kepler's equation for eccentric anomaly
    
    Parameters:
    -----------
    M : float or array-like
        Mean anomaly [rad]
    e : float
        Eccentricity
    tolerance : float, optional
        Numerical tolerance for solution
    
    Returns:
    --------
    float or array-like
        Eccentric anomaly E [rad]
    """
    if np.isscalar(M):
        # Initial guess
        E0 = M if e < 0.8 else np.pi
        return newton(kepler_equation, E0, args=(M, e), tol=tolerance)
    else:
        # Handle array input
        E_array = np.zeros_like(M)
        for i, Mi in enumerate(M):
            E0 = Mi if e < 0.8 else np.pi
            E_array[i] = newton(kepler_equation, E0, args=(Mi, e), tol=tolerance)
        return E_array


def true_anomaly_from_eccentric(E, e):
    """
    Calculate true anomaly from eccentric anomaly
    
    Parameters:
    -----------
    E : float or array-like
        Eccentric anomaly [rad]
    e : float
        Eccentricity
    
    Returns:
    --------
    float or array-like
        True anomaly [rad]
    """
    return 2 * np.arctan2(
        np.sqrt(1 + e) * np.sin(E / 2),
        np.sqrt(1 - e) * np.cos(E / 2)
    )


def orbital_elements_to_state(a, e, M, mu):
    """
    Convert orbital elements to position and velocity vectors (classical mechanics)
    
    Parameters:
    -----------
    a : float
        Semi-major axis [m]
    e : float
        Eccentricity
    M : float
        Mean anomaly [rad]
    mu : float
        Standard gravitational parameter μ = G*(m1+m2) [m^3/s^2]
    
    Returns:
    --------
    tuple
        (position, velocity) where both are numpy arrays of shape (3,)
        position [m], velocity [m/s]
    """
    # Solve Kepler's equation
    E = solve_kepler(M, e)
    
    # Calculate distance and true anomaly
    r = a * (1 - e * np.cos(E))
    nu = true_anomaly_from_eccentric(E, e)
    
    # Position in orbital plane
    x = r * np.cos(nu)
    y = r * np.sin(nu)
    z = 0.0
    position = np.array([x, y, z])
    
    # Velocity in orbital plane
    h = np.sqrt(mu * a * (1 - e**2))  # Specific angular momentum
    vx = -mu / h * np.sin(nu)
    vy = mu / h * (e + np.cos(nu))
    vz = 0.0
    velocity = np.array([vx, vy, vz])
    
    return position, velocity


def orbital_period(a, mu):
    """
    Calculate orbital period using Kepler's third law
    
    Parameters:
    -----------
    a : float
        Semi-major axis [m]
    mu : float
        Standard gravitational parameter [m^3/s^2]
    
    Returns:
    --------
    float
        Orbital period [s]
    """
    return 2 * np.pi * np.sqrt(a**3 / mu)


def mean_motion(a, mu):
    """
    Calculate mean motion
    
    Parameters:
    -----------
    a : float
        Semi-major axis [m]
    mu : float
        Standard gravitational parameter [m^3/s^2]
    
    Returns:
    --------
    float
        Mean motion [rad/s]
    """
    return np.sqrt(mu / a**3)


def vis_viva_equation(r, a, mu):
    """
    Calculate orbital velocity using vis-viva equation
    
    Parameters:
    -----------
    r : float
        Distance from central body [m]
    a : float
        Semi-major axis [m] (positive for ellipse, negative for hyperbola)
    mu : float
        Standard gravitational parameter [m^3/s^2]
    
    Returns:
    --------
    float
        Orbital speed [m/s]
    """
    return np.sqrt(mu * (2/r - 1/a))


def specific_energy(a, mu):
    """
    Calculate specific orbital energy
    
    Parameters:
    -----------
    a : float
        Semi-major axis [m]
    mu : float
        Standard gravitational parameter [m^3/s^2]
    
    Returns:
    --------
    float
        Specific energy [J/kg]
    """
    return -mu / (2 * a)


def escape_velocity(r, mu):
    """
    Calculate escape velocity at given distance
    
    Parameters:
    -----------
    r : float
        Distance from central body [m]
    mu : float
        Standard gravitational parameter [m^3/s^2]
    
    Returns:
    --------
    float
        Escape velocity [m/s]
    """
    return np.sqrt(2 * mu / r)


def orbit_from_state_vectors(position, velocity, mu):
    """
    Calculate orbital elements from state vectors
    
    Parameters:
    -----------
    position : array-like
        Position vector [m]
    velocity : array-like
        Velocity vector [m/s]
    mu : float
        Standard gravitational parameter [m^3/s^2]
    
    Returns:
    --------
    dict
        Orbital elements: {'a': semi-major axis, 'e': eccentricity, 'M': mean anomaly}
    """
    r_vec = np.array(position)
    v_vec = np.array(velocity)
    
    r = np.linalg.norm(r_vec)
    v = np.linalg.norm(v_vec)
    
    # Specific energy
    epsilon = v**2 / 2 - mu / r
    
    # Semi-major axis
    a = -mu / (2 * epsilon)
    
    # Eccentricity vector
    e_vec = ((v**2 - mu/r) * r_vec - np.dot(r_vec, v_vec) * v_vec) / mu
    e = np.linalg.norm(e_vec)
    
    # Calculate anomalies (simplified for circular/elliptical orbits)
    if e < 1e-8:  # Nearly circular
        E = 0.0
        M = 0.0
    else:
        cos_E = (1 - r/a) / e
        E = np.arccos(np.clip(cos_E, -1, 1))
        M = E - e * np.sin(E)
    
    return {
        'a': a,
        'e': e,
        'E': E,
        'M': M
    }