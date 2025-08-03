"""
Unit tests for classical celestial mechanics module
"""
import unittest
import numpy as np
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from classical import (
    kepler_equation,
    solve_kepler,
    true_anomaly_from_eccentric,
    orbital_elements_to_state,
    orbital_period,
    mean_motion,
    vis_viva_equation,
    specific_energy,
    escape_velocity
)
from utils import G, AU, solar_mass


class TestClassicalMechanics(unittest.TestCase):
    
    def setUp(self):
        """Set up test parameters"""
        self.tolerance = 1e-10
        self.mu_sun = G * solar_mass
        
    def test_kepler_equation_zero_eccentricity(self):
        """Test Kepler equation with circular orbit (e=0)"""
        M = np.pi / 2  # 90 degrees
        e = 0.0
        E = solve_kepler(M, e)
        self.assertAlmostEqual(E, M, places=10)
        self.assertAlmostEqual(kepler_equation(E, M, e), 0.0, places=10)
    
    def test_kepler_equation_moderate_eccentricity(self):
        """Test Kepler equation with moderate eccentricity"""
        M = np.pi / 3  # 60 degrees
        e = 0.5
        E = solve_kepler(M, e)
        residual = kepler_equation(E, M, e)
        self.assertLess(abs(residual), self.tolerance)
    
    def test_kepler_equation_high_eccentricity(self):
        """Test Kepler equation with high eccentricity"""
        M = np.pi / 4  # 45 degrees
        e = 0.9
        E = solve_kepler(M, e)
        residual = kepler_equation(E, M, e)
        self.assertLess(abs(residual), self.tolerance)
    
    def test_kepler_equation_array_input(self):
        """Test Kepler equation with array input"""
        M_array = np.array([0, np.pi/4, np.pi/2, np.pi])
        e = 0.3
        E_array = solve_kepler(M_array, e)
        
        for M, E in zip(M_array, E_array):
            residual = kepler_equation(E, M, e)
            self.assertLess(abs(residual), self.tolerance)
    
    def test_true_anomaly_circular(self):
        """Test true anomaly for circular orbit"""
        E = np.pi / 2
        e = 0.0
        nu = true_anomaly_from_eccentric(E, e)
        self.assertAlmostEqual(nu, E, places=10)
    
    def test_true_anomaly_eccentric(self):
        """Test true anomaly for eccentric orbit"""
        E = np.pi / 2
        e = 0.5
        nu = true_anomaly_from_eccentric(E, e)
        # For e=0.5 and E=π/2, ν should be approximately 2.094 radians
        self.assertAlmostEqual(nu, 2.0943951023931953, places=6)
    
    def test_orbital_elements_to_state_circular(self):
        """Test state vector conversion for circular orbit"""
        a = 1.0 * AU
        e = 0.0
        M = 0.0  # At periapsis
        
        pos, vel = orbital_elements_to_state(a, e, M, self.mu_sun)
        
        # Should be at distance a along x-axis
        self.assertAlmostEqual(np.linalg.norm(pos), a, places=6)
        self.assertAlmostEqual(pos[0], a, places=6)
        self.assertAlmostEqual(pos[1], 0.0, places=10)
        
        # Velocity should be purely in y-direction
        expected_v = np.sqrt(self.mu_sun / a)
        self.assertAlmostEqual(np.linalg.norm(vel), expected_v, places=6)
        self.assertAlmostEqual(vel[0], 0.0, places=10)
    
    def test_orbital_elements_to_state_eccentric(self):
        """Test state vector conversion for eccentric orbit"""
        a = 1.0 * AU
        e = 0.5
        M = 0.0  # At periapsis
        
        pos, vel = orbital_elements_to_state(a, e, M, self.mu_sun)
        
        # At periapsis, distance should be a(1-e)
        expected_r = a * (1 - e)
        self.assertAlmostEqual(np.linalg.norm(pos), expected_r, places=6)
    
    def test_orbital_period_earth(self):
        """Test orbital period calculation for Earth"""
        a_earth = 1.0 * AU
        T = orbital_period(a_earth, self.mu_sun)
        T_years = T / (365.25 * 24 * 3600)
        
        # Should be approximately 1 year
        self.assertAlmostEqual(T_years, 1.0, places=3)
    
    def test_mean_motion(self):
        """Test mean motion calculation"""
        a = 1.0 * AU
        n = mean_motion(a, self.mu_sun)
        T = orbital_period(a, self.mu_sun)
        
        # n = 2π/T
        expected_n = 2 * np.pi / T
        self.assertAlmostEqual(n, expected_n, places=10)
    
    def test_vis_viva_equation_circular(self):
        """Test vis-viva equation for circular orbit"""
        r = 1.0 * AU
        a = 1.0 * AU  # circular orbit
        
        v = vis_viva_equation(r, a, self.mu_sun)
        expected_v = np.sqrt(self.mu_sun / r)
        
        self.assertAlmostEqual(v, expected_v, places=10)
    
    def test_vis_viva_equation_periapsis(self):
        """Test vis-viva equation at periapsis"""
        a = 1.0 * AU
        e = 0.5
        r_peri = a * (1 - e)
        
        v_peri = vis_viva_equation(r_peri, a, self.mu_sun)
        
        # At periapsis: v = sqrt(μ(1+e)/(a(1-e)))
        expected_v = np.sqrt(self.mu_sun * (1 + e) / (a * (1 - e)))
        self.assertAlmostEqual(v_peri, expected_v, places=10)
    
    def test_specific_energy(self):
        """Test specific energy calculation"""
        a = 1.0 * AU
        epsilon = specific_energy(a, self.mu_sun)
        
        # ε = -μ/(2a)
        expected_epsilon = -self.mu_sun / (2 * a)
        self.assertAlmostEqual(epsilon, expected_epsilon, places=10)
    
    def test_escape_velocity_earth_surface(self):
        """Test escape velocity from Earth's surface"""
        R_earth = 6.371e6  # meters
        M_earth = 5.972e24  # kg
        mu_earth = G * M_earth
        
        v_esc = escape_velocity(R_earth, mu_earth)
        
        # Should be approximately 11.2 km/s
        expected_v_esc = 11186  # m/s
        self.assertAlmostEqual(v_esc, expected_v_esc, places=0)
    
    def test_energy_conservation(self):
        """Test energy conservation in orbit"""
        a = 1.5 * AU
        e = 0.3
        
        # Calculate specific energy
        epsilon = specific_energy(a, self.mu_sun)
        
        # Test at different positions
        for M in [0, np.pi/2, np.pi]:
            pos, vel = orbital_elements_to_state(a, e, M, self.mu_sun)
            r = np.linalg.norm(pos)
            v = np.linalg.norm(vel)
            
            # Calculate energy: ε = v²/2 - μ/r
            energy = v**2 / 2 - self.mu_sun / r
            
            self.assertAlmostEqual(energy, epsilon, places=6)
    
    def test_angular_momentum_conservation(self):
        """Test angular momentum conservation in orbit"""
        a = 1.0 * AU
        e = 0.4
        
        # Specific angular momentum: h = sqrt(μa(1-e²))
        h_expected = np.sqrt(self.mu_sun * a * (1 - e**2))
        
        # Test at different positions
        for M in [0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi]:
            pos, vel = orbital_elements_to_state(a, e, M, self.mu_sun)
            
            # Calculate angular momentum: h = |r × v|
            h_calculated = np.linalg.norm(np.cross(pos, vel))
            
            # Test relative error is small
            relative_error = abs(h_calculated - h_expected) / h_expected
            self.assertLess(relative_error, 1e-10)


class TestKepler3rdLaw(unittest.TestCase):
    """Test Kepler's third law with known planetary data"""
    
    def setUp(self):
        self.mu_sun = G * solar_mass
        self.tolerance_percent = 1.0  # 1% tolerance for real planetary data
    
    def test_mercury(self):
        """Test Kepler's 3rd law for Mercury"""
        a_mercury = 0.387 * AU
        T_expected = 87.97 * 24 * 3600  # seconds
        
        T_calculated = orbital_period(a_mercury, self.mu_sun)
        error_percent = abs(T_calculated - T_expected) / T_expected * 100
        
        self.assertLess(error_percent, self.tolerance_percent)
    
    def test_earth(self):
        """Test Kepler's 3rd law for Earth"""
        a_earth = 1.0 * AU
        T_expected = 365.25 * 24 * 3600  # seconds
        
        T_calculated = orbital_period(a_earth, self.mu_sun)
        error_percent = abs(T_calculated - T_expected) / T_expected * 100
        
        self.assertLess(error_percent, self.tolerance_percent)
    
    def test_mars(self):
        """Test Kepler's 3rd law for Mars"""
        a_mars = 1.524 * AU
        T_expected = 687 * 24 * 3600  # seconds
        
        T_calculated = orbital_period(a_mars, self.mu_sun)
        error_percent = abs(T_calculated - T_expected) / T_expected * 100
        
        self.assertLess(error_percent, 2.0)  # Slightly higher tolerance for Mars


if __name__ == '__main__':
    unittest.main()