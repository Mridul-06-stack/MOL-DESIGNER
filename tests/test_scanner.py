"""Tests for moldesigner.scanner and red-team docking integration."""

import pytest
from rdkit import Chem

from moldesigner.scanner import MockScanner
from moldesigner.docking import RDKitScoreDocker

ERLOTINIB = "C=Cc1cccc(Nc2ncnc3cc(OCCOC)c(OCCOC)cc23)c1"
ASPIRIN = "CC(=O)Oc1ccccc1C(=O)O"

def test_mock_scanner_flags_good_affinity():
    scanner = MockScanner()
    # Baseline energy good (-9.0)
    variant = scanner.scan(ERLOTINIB, -9.0)
    assert variant is not None
    assert variant.startswith("RESIST_SIM_")
    assert ERLOTINIB in variant

def test_mock_scanner_ignores_bad_affinity():
    scanner = MockScanner()
    # Baseline energy bad (-3.0)
    variant = scanner.scan(ERLOTINIB, -3.0)
    assert variant is None

def test_redteam_variant_docking_penalty():
    """Verify that RDKitScoreDocker properly penalizes SIM variants."""
    scanner = MockScanner()
    variant = scanner.scan(ERLOTINIB, -9.0)
    
    docker = RDKitScoreDocker(targets=["WT", variant])
    
    results = docker.dock_many([ERLOTINIB, ASPIRIN])
    assert len(results) == 2
    
    erlotinib_res, aspirin_res = results
    
    # Erlotinib should be heavily penalized by the specific variant (energy > -4)
    assert erlotinib_res["WT"] < -6.0  # good on WT
    assert erlotinib_res[variant] > -4.0  # bad on escape mutation
    
    # Aspirin is structurally distinct, so it shouldn't trigger the red-team penalty logic
    # Its baseline might be poor because of RDKit properties, but WT and Variant 
    # should be roughly identical on the fallback calculation.
    import math
    assert math.isclose(aspirin_res["WT"], aspirin_res[variant], abs_tol=1.5)
