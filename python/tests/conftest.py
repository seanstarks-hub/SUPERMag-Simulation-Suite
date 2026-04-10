"""Shared pytest fixtures for SUPERMag test suite."""

import os
import sys

# Add python/ directory to path so 'supermag' package is importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

import pytest
import numpy as np


@pytest.fixture
def nb_fe_params():
    """Standard Nb/Fe bilayer parameters."""
    return {
        "Tc0": 9.2,       # K
        "d_S": 50.0,      # nm
        "xi_S": 38.0,     # nm
        "xi_F": 0.7,      # nm
        "E_ex": 256.0,    # meV
    }


@pytest.fixture
def d_F_array():
    """Standard d_F sweep array (nm)."""
    return np.linspace(0.5, 20.0, 50)
