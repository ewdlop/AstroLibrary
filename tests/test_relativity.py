"""
Unit tests for general relativity module
"""
import unittest
import numpy as np
import sys
import os

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from relativity import (
    pericenter_precession,
    pericenter_precession_per_year,
    mercury_precession_test,
    shapiro_delay,
    schwarzschild_radius,
    gravitational_redshift,
    light_deflection_angle,
    time_dilation_factor,
    effective_potential,
    post_newtonian_correction
)
from utils import G, c, AU, solar_mass, solar_radius


class TestRelativisticEffects(unittest.TestCase):
    
    def setUp(self):
        """Set up test parameters"""
        self.M_sun = solar_mass
        self.tolerance_percent = 10.0  # 10% tolerance for theoretical calculations
    
    def test_schwarzschild_radius_sun(self):
        """Test Schwarzschild radius calculation for the Sun"""
        rs = schwarzschild_radius(self.M_sun)
        
        # Expected value: ~2.95 km
        expected_rs = 2 * G * self.M_sun / c**2
        self.assertAlmostEqual(rs, expected_rs, places=10)
        
        # Check approximate value
        self.assertAlmostEqual(rs / 1000, 2.95, places=1)  # ~2.95 km
    
    def test_schwarzschild_radius_earth(self):
        """Test Schwarzschild radius calculation for Earth"""
        M_earth = 5.972e24  # kg
        rs = schwarzschild_radius(M_earth)
        
        # Expected value: ~8.87 mm
        self.assertAlmostEqual(rs * 1000, 8.87, places=1)  # ~8.87 mm
    
    def test_mercury_precession_order_of_magnitude(self):
        """Test Mercury's precession is correct order of magnitude"""
        a_mercury = 5.791e10  # meters
        e_mercury = 0.2056
        
        precession_per_orbit = pericenter_precession(a_mercury, e_mercury, self.M_sun)
        
        # Should be on the order of 10^-7 radians per orbit
        self.assertGreater(precession_per_orbit, 1e-8)
        self.assertLess(precession_per_orbit, 1e-6)
    
    def test_mercury_precession_known_value(self):
        """Test Mercury's precession against known observational value"""
        result = mercury_precession_test()
        
        # Should be within reasonable agreement with observations
        self.assertLess(result['relative_error'], 0.5)  # Less than 50% error
        
        # Check that calculated value is in right ballpark (30-50 arcsec/century)
        calculated_per_century = result['calculated_gr_precession_per_year'] * 100
        self.assertGreater(calculated_per_century, 30)
        self.assertLess(calculated_per_century, 50)
    
    def test_light_deflection_solar_limb(self):
        """Test light deflection at solar limb"""
        deflection = light_deflection_angle(solar_radius, self.M_sun)
        deflection_arcsec = deflection * 180 * 3600 / np.pi
        
        # Should be approximately 1.75 arcseconds
        self.assertAlmostEqual(deflection_arcsec, 1.75, places=1)
    
    def test_light_deflection_scaling(self):
        """Test light deflection scales correctly with impact parameter"""
        b1 = solar_radius
        b2 = 2 * solar_radius
        
        deflection1 = light_deflection_angle(b1, self.M_sun)
        deflection2 = light_deflection_angle(b2, self.M_sun)
        
        # Should scale as 1/b (inverse relationship)
        ratio = deflection1 / deflection2
        self.assertAlmostEqual(ratio, 2.0, places=6)
    
    def test_shapiro_delay_positive(self):
        """Test that Shapiro delay is always positive"""
        r1 = 1.0 * AU
        r2 = 1.5 * AU
        r_closest = solar_radius
        
        delay = shapiro_delay(r1, r2, r_closest, self.M_sun)
        
        self.assertGreater(delay, 0)
        self.assertLess(delay, 1e-3)  # Should be less than 1 ms for typical values
    
    def test_gravitational_redshift_weak_field(self):
        """Test gravitational redshift in weak field limit"""
        r = 10 * solar_radius  # Far from the Sun
        
        redshift_factor = gravitational_redshift(r, self.M_sun)
        
        # In weak field: (1+z) ≈ 1 + GM/(rc²)
        weak_field_approx = 1 + G * self.M_sun / (r * c**2)
        
        self.assertAlmostEqual(redshift_factor, weak_field_approx, places=6)
    
    def test_time_dilation_factor_weak_field(self):
        """Test time dilation factor in weak field"""
        r = 10 * solar_radius
        
        factor = time_dilation_factor(r, self.M_sun)
        
        # Should be slightly less than 1
        self.assertLess(factor, 1.0)
        self.assertGreater(factor, 0.99)  # But close to 1 for weak field
    
    def test_effective_potential_basic(self):
        """Test basic effective potential calculation"""
        l = 1e13  # Angular momentum per unit mass [m²/s]
        r = 10 * solar_radius
        
        potential = effective_potential(r, l, self.M_sun)
        
        # Should be a finite number
        self.assertTrue(np.isfinite(potential))
        self.assertIsInstance(potential, (int, float, np.number))
    
    def test_post_newtonian_correction_weak_field(self):
        """Test post-Newtonian correction in weak field"""
        r = 1.0 * AU
        v = 30000  # m/s (typical orbital velocity)
        
        correction = post_newtonian_correction(r, v, self.M_sun)
        
        # Should be a small correction (1 + small number)
        self.assertGreater(correction, 1.0)
        self.assertLess(correction, 1.001)  # Less than 0.1% correction


class TestRelativisticScaling(unittest.TestCase):
    """Test scaling relationships for relativistic effects"""
    
    def setUp(self):
        self.M_sun = solar_mass
    
    def test_precession_mass_scaling(self):
        """Test that precession scales linearly with mass"""
        a = 1.0 * AU
        e = 0.1
        
        M1 = self.M_sun
        M2 = 2 * self.M_sun
        
        prec1 = pericenter_precession(a, e, M1)
        prec2 = pericenter_precession(a, e, M2)
        
        ratio = prec2 / prec1
        self.assertAlmostEqual(ratio, 2.0, places=6)
    
    def test_precession_semimajor_axis_scaling(self):
        """Test that precession scales as 1/a"""
        e = 0.1
        M = self.M_sun
        
        a1 = 1.0 * AU
        a2 = 2.0 * AU
        
        prec1 = pericenter_precession(a1, e, M)
        prec2 = pericenter_precession(a2, e, M)
        
        ratio = prec1 / prec2
        self.assertAlmostEqual(ratio, 2.0, places=6)
    
    def test_precession_eccentricity_scaling(self):
        """Test that precession scales as 1/(1-e²)"""
        a = 1.0 * AU
        M = self.M_sun
        
        e1 = 0.0
        e2 = 0.5
        
        prec1 = pericenter_precession(a, e1, M)
        prec2 = pericenter_precession(a, e2, M)
        
        # For e1=0: (1-e²) = 1
        # For e2=0.5: (1-e²) = 0.75
        expected_ratio = 1.0 / 0.75
        actual_ratio = prec2 / prec1
        
        self.assertAlmostEqual(actual_ratio, expected_ratio, places=6)
    
    def test_schwarzschild_radius_mass_scaling(self):
        """Test that Schwarzschild radius scales linearly with mass"""
        M1 = self.M_sun
        M2 = 3 * self.M_sun
        
        rs1 = schwarzschild_radius(M1)
        rs2 = schwarzschild_radius(M2)
        
        ratio = rs2 / rs1
        self.assertAlmostEqual(ratio, 3.0, places=10)


class TestRelativisticConstants(unittest.TestCase):
    """Test that relativistic calculations use correct physical constants"""
    
    def test_speed_of_light_value(self):
        """Test that speed of light has correct value"""
        # c should be exactly 299,792,458 m/s by definition
        self.assertEqual(c, 299792458.0)
    
    def test_gravitational_constant_reasonable(self):
        """Test that gravitational constant is reasonable"""
        # G should be approximately 6.674e-11
        self.assertAlmostEqual(G / 1e-11, 6.674, places=1)
    
    def test_solar_mass_reasonable(self):
        """Test that solar mass is reasonable"""
        # Solar mass should be approximately 1.989e30 kg
        self.assertAlmostEqual(solar_mass / 1e30, 1.989, places=2)


class TestPhysicalConsistency(unittest.TestCase):
    """Test physical consistency of relativistic calculations"""
    
    def setUp(self):
        self.M_sun = solar_mass
    
    def test_precession_always_positive(self):
        """Test that precession is always positive"""
        a_values = [0.3*AU, 1.0*AU, 5.0*AU]
        e_values = [0.0, 0.2, 0.5, 0.8]
        
        for a in a_values:
            for e in e_values:
                prec = pericenter_precession(a, e, self.M_sun)
                self.assertGreater(prec, 0)
    
    def test_time_dilation_less_than_one(self):
        """Test that gravitational time dilation factor is ≤ 1"""
        r_values = np.logspace(0, 3, 20) * solar_radius  # From 1 to 1000 solar radii
        
        for r in r_values:
            if r > schwarzschild_radius(self.M_sun):  # Outside event horizon
                factor = time_dilation_factor(r, self.M_sun)
                self.assertLessEqual(factor, 1.0)
                self.assertGreater(factor, 0.0)
    
    def test_redshift_factor_greater_than_one(self):
        """Test that gravitational redshift factor is ≥ 1"""
        r_values = np.logspace(0, 3, 20) * solar_radius
        
        for r in r_values:
            if r > schwarzschild_radius(self.M_sun):
                factor = gravitational_redshift(r, self.M_sun)
                self.assertGreaterEqual(factor, 1.0)


if __name__ == '__main__':
    unittest.main()