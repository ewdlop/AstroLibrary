"""
Unit tests for utilities module
"""
import unittest
import numpy as np
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from utils import *


class TestUnitConversions(unittest.TestCase):
    """Test unit conversion functions"""
    
    def test_angle_conversions(self):
        """Test angle conversion functions"""
        # Test degrees to radians
        self.assertAlmostEqual(deg_to_rad(180), np.pi, places=10)
        self.assertAlmostEqual(deg_to_rad(90), np.pi/2, places=10)
        self.assertAlmostEqual(deg_to_rad(0), 0, places=10)
        
        # Test radians to degrees
        self.assertAlmostEqual(rad_to_deg(np.pi), 180, places=10)
        self.assertAlmostEqual(rad_to_deg(np.pi/2), 90, places=10)
        self.assertAlmostEqual(rad_to_deg(0), 0, places=10)
        
        # Test roundtrip conversions
        angles_deg = [0, 30, 45, 90, 180, 270, 360]
        for angle in angles_deg:
            self.assertAlmostEqual(rad_to_deg(deg_to_rad(angle)), angle, places=10)
    
    def test_arcsec_conversions(self):
        """Test arcsecond conversion functions"""
        # 1 degree = 3600 arcseconds
        self.assertAlmostEqual(arcsec_to_rad(3600), deg_to_rad(1), places=10)
        self.assertAlmostEqual(rad_to_arcsec(deg_to_rad(1)), 3600, places=10)
        
        # Test roundtrip
        arcsec_values = [1, 10, 100, 3600]
        for arcsec in arcsec_values:
            self.assertAlmostEqual(rad_to_arcsec(arcsec_to_rad(arcsec)), arcsec, places=10)
    
    def test_distance_conversions(self):
        """Test distance conversion functions"""
        # Test AU conversions
        self.assertAlmostEqual(au_to_meters(1), AU, places=10)
        self.assertAlmostEqual(meters_to_au(AU), 1, places=10)
        
        # Test roundtrip
        au_values = [0.1, 1.0, 5.2, 30.1]
        for au_val in au_values:
            self.assertAlmostEqual(meters_to_au(au_to_meters(au_val)), au_val, places=10)
    
    def test_mass_conversions(self):
        """Test mass conversion functions"""
        # Test solar mass conversions
        self.assertAlmostEqual(solar_masses_to_kg(1), solar_mass, places=10)
        self.assertAlmostEqual(kg_to_solar_masses(solar_mass), 1, places=10)
        
        # Test roundtrip
        mass_values = [0.1, 1.0, 10.0, 100.0]
        for mass in mass_values:
            self.assertAlmostEqual(kg_to_solar_masses(solar_masses_to_kg(mass)), mass, places=10)
    
    def test_time_conversions(self):
        """Test time conversion functions"""
        # Test year conversions
        self.assertAlmostEqual(years_to_seconds(1), seconds_per_year, places=10)
        self.assertAlmostEqual(seconds_to_years(seconds_per_year), 1, places=10)
        
        # Test day conversions
        self.assertAlmostEqual(days_to_seconds(1), seconds_per_day, places=10)
        self.assertAlmostEqual(seconds_to_days(seconds_per_day), 1, places=10)
        
        # Test roundtrip
        time_values = [1, 10, 365.25, 1000]
        for time_val in time_values:
            self.assertAlmostEqual(seconds_to_years(years_to_seconds(time_val)), time_val, places=10)
            self.assertAlmostEqual(seconds_to_days(days_to_seconds(time_val)), time_val, places=10)
    
    def test_velocity_conversions(self):
        """Test velocity conversion functions"""
        # Test km/s to m/s
        self.assertAlmostEqual(km_per_s_to_m_per_s(1), 1000, places=10)
        self.assertAlmostEqual(m_per_s_to_km_per_s(1000), 1, places=10)
        
        # Test roundtrip
        vel_values = [0.1, 1.0, 30.0, 300.0]
        for vel in vel_values:
            self.assertAlmostEqual(m_per_s_to_km_per_s(km_per_s_to_m_per_s(vel)), vel, places=10)


class TestPhysicalConstants(unittest.TestCase):
    """Test that physical constants have reasonable values"""
    
    def test_fundamental_constants(self):
        """Test fundamental physical constants"""
        # Gravitational constant (approximate value)
        self.assertAlmostEqual(G / 1e-11, 6.674, places=1)
        
        # Speed of light (exact by definition)
        self.assertEqual(c, 299792458.0)
        
        # Planck constant (SI 2019 definition)
        self.assertEqual(h, 6.62607015e-34)
        
        # Boltzmann constant (SI 2019 definition)
        self.assertEqual(k_B, 1.380649e-23)
    
    def test_astronomical_constants(self):
        """Test astronomical constants"""
        # AU should be approximately 1.496e11 meters
        self.assertAlmostEqual(AU / 1e11, 1.496, places=2)
        
        # Solar mass should be approximately 1.989e30 kg
        self.assertAlmostEqual(solar_mass / 1e30, 1.989, places=2)
        
        # Earth mass should be approximately 5.972e24 kg
        self.assertAlmostEqual(earth_mass / 1e24, 5.972, places=2)


class TestOrbitalCalculations(unittest.TestCase):
    """Test orbital calculation utility functions"""
    
    def setUp(self):
        self.mu_earth = standard_gravitational_parameter(earth_mass)
        self.R_earth = earth_radius
    
    def test_standard_gravitational_parameter(self):
        """Test standard gravitational parameter calculation"""
        mu = standard_gravitational_parameter(earth_mass)
        expected_mu = G * earth_mass
        self.assertAlmostEqual(mu, expected_mu, places=10)
    
    def test_circular_orbital_velocity(self):
        """Test circular orbital velocity calculation"""
        # Low Earth orbit (~400 km altitude)
        r = self.R_earth + 400e3  # 400 km altitude
        v = orbital_velocity_circular(r, self.mu_earth)
        
        # Should be approximately 7.7 km/s
        self.assertAlmostEqual(v / 1000, 7.67, places=1)
    
    def test_escape_velocity_earth(self):
        """Test escape velocity from Earth's surface"""
        v_esc = escape_velocity(self.R_earth, self.mu_earth)
        
        # Should be approximately 11.2 km/s
        self.assertAlmostEqual(v_esc / 1000, 11.18, places=1)
    
    def test_sphere_of_influence(self):
        """Test sphere of influence calculation"""
        # Earth-Moon system
        m_earth = earth_mass
        m_moon = 7.342e22  # kg
        a_moon = 384400e3  # meters
        
        soi = sphere_of_influence(m_earth, m_moon, a_moon)
        
        # Should be approximately 66,000 km
        self.assertAlmostEqual(soi / 1000, 66100, places=-3)  # Within 1000 km
    
    def test_hill_sphere(self):
        """Test Hill sphere calculation"""
        # Earth around Sun
        m_sun = solar_mass
        m_earth = earth_mass
        a_earth = AU
        
        hill_r = hill_sphere(m_sun, m_earth, a_earth)
        
        # Should be approximately 1.5 million km
        self.assertAlmostEqual(hill_r / 1000, 1.5e6, places=-4)  # Within 10,000 km


class TestAstronomicalFunctions(unittest.TestCase):
    """Test astronomical calculation functions"""
    
    def test_angular_separation(self):
        """Test angular separation calculation"""
        # Same point
        ra1, dec1 = 0, 0
        ra2, dec2 = 0, 0
        sep = angular_separation(ra1, dec1, ra2, dec2)
        self.assertAlmostEqual(sep, 0, places=10)
        
        # 90 degree separation
        ra1, dec1 = 0, 0
        ra2, dec2 = np.pi/2, 0
        sep = angular_separation(ra1, dec1, ra2, dec2)
        self.assertAlmostEqual(sep, np.pi/2, places=10)
        
        # Antipodal points
        ra1, dec1 = 0, 0
        ra2, dec2 = np.pi, 0
        sep = angular_separation(ra1, dec1, ra2, dec2)
        self.assertAlmostEqual(sep, np.pi, places=10)
    
    def test_magnitude_flux_conversions(self):
        """Test magnitude-flux conversion functions"""
        # Test roundtrip conversions
        magnitudes = [0, 5, 10, 15, 20]
        for mag in magnitudes:
            flux = magnitude_to_flux(mag)
            mag_back = flux_to_magnitude(flux)
            self.assertAlmostEqual(mag_back, mag, places=10)
        
        # Test magnitude difference of 5 corresponds to flux ratio of 100
        mag1, mag2 = 10, 15
        flux1 = magnitude_to_flux(mag1)
        flux2 = magnitude_to_flux(mag2)
        flux_ratio = flux1 / flux2
        self.assertAlmostEqual(flux_ratio, 100, places=6)
    
    def test_parallax_distance(self):
        """Test parallax to distance conversion"""
        # 1 arcsecond parallax = 1 parsec distance
        distance = parallax_to_distance(1.0)
        self.assertAlmostEqual(distance, 1.0, places=10)
        
        # Proxima Centauri: ~0.768 arcsec parallax
        distance_proxima = parallax_to_distance(0.768)
        self.assertAlmostEqual(distance_proxima, 1.302, places=2)  # ~1.3 pc
    
    def test_blackbody_temperature(self):
        """Test Wien's law for blackbody temperature"""
        # Sun's peak wavelength ~500 nm
        wavelength_peak = 500e-9  # meters
        T = blackbody_temperature(wavelength_peak)
        
        # Should be approximately 5800 K
        self.assertAlmostEqual(T, 5795, places=-1)


class TestMathematicalUtilities(unittest.TestCase):
    """Test mathematical utility functions"""
    
    def test_rotation_matrix_z(self):
        """Test z-axis rotation matrix"""
        # 90 degree rotation
        angle = np.pi / 2
        R = rotation_matrix_z(angle)
        
        # Test that it rotates (1,0,0) to (0,1,0)
        v_in = np.array([1, 0, 0])
        v_out = R @ v_in
        expected = np.array([0, 1, 0])
        
        np.testing.assert_array_almost_equal(v_out, expected, decimal=10)
        
        # Test that it's a proper rotation matrix (det = 1, orthogonal)
        self.assertAlmostEqual(np.linalg.det(R), 1.0, places=10)
        np.testing.assert_array_almost_equal(R @ R.T, np.eye(3), decimal=10)
    
    def test_coordinate_transform_2d(self):
        """Test 2D coordinate transformation"""
        # 90 degree rotation
        x, y = 1, 0
        angle = np.pi / 2
        x_new, y_new = coordinate_transform_2d(x, y, angle)
        
        self.assertAlmostEqual(x_new, 0, places=10)
        self.assertAlmostEqual(y_new, 1, places=10)
        
        # Test with arrays
        x_arr = np.array([1, 0, -1, 0])
        y_arr = np.array([0, 1, 0, -1])
        x_new_arr, y_new_arr = coordinate_transform_2d(x_arr, y_arr, np.pi/2)
        
        expected_x = np.array([0, -1, 0, 1])
        expected_y = np.array([1, 0, -1, 0])
        
        np.testing.assert_array_almost_equal(x_new_arr, expected_x, decimal=10)
        np.testing.assert_array_almost_equal(y_new_arr, expected_y, decimal=10)
    
    def test_create_time_array(self):
        """Test time array creation"""
        start, end, num = 0.0, 10.0, 11
        t_array = create_time_array(start, end, num)
        
        self.assertEqual(len(t_array), num)
        self.assertAlmostEqual(t_array[0], start, places=10)
        self.assertAlmostEqual(t_array[-1], end, places=10)
        
        # Check spacing
        expected_spacing = (end - start) / (num - 1)
        actual_spacing = t_array[1] - t_array[0]
        self.assertAlmostEqual(actual_spacing, expected_spacing, places=10)


class TestPhysicalConsistency(unittest.TestCase):
    """Test physical consistency of utility functions"""
    
    def test_escape_velocity_greater_than_circular(self):
        """Test that escape velocity > circular velocity"""
        r = earth_radius
        mu = standard_gravitational_parameter(earth_mass)
        
        v_circular = orbital_velocity_circular(r, mu)
        v_escape = escape_velocity(r, mu)
        
        self.assertGreater(v_escape, v_circular)
        # Escape velocity should be sqrt(2) times circular velocity
        self.assertAlmostEqual(v_escape / v_circular, np.sqrt(2), places=6)
    
    def test_synodic_period_consistency(self):
        """Test synodic period calculation"""
        # Earth and Mars
        P_earth = 1.0 * julian_year  # seconds
        P_mars = 1.88 * julian_year  # seconds
        
        P_synodic = synodic_period(P_earth, P_mars)
        
        # Synodic period should be longer than both orbital periods
        self.assertGreater(P_synodic, P_earth)
        self.assertGreater(P_synodic, P_mars)
        
        # Should be approximately 2.14 years for Earth-Mars
        P_synodic_years = P_synodic / julian_year
        self.assertAlmostEqual(P_synodic_years, 2.14, places=1)


if __name__ == '__main__':
    unittest.main()