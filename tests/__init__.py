"""
Test suite for AstroLibrary celestial mechanics package
"""

import unittest
import sys
import os

# Add src directory to path for testing
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from .test_classical import *
from .test_relativity import *
from .test_utils import *


def run_all_tests():
    """Run all tests in the test suite"""
    # Discover and run all tests
    loader = unittest.TestLoader()
    suite = loader.discover(os.path.dirname(__file__), pattern='test_*.py')
    
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result


if __name__ == '__main__':
    run_all_tests()