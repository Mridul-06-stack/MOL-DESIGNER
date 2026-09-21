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


def test_clinical_scanner_progression():
    from moldesigner.scanner import ClinicalEGFRScanner, CLINICAL_EGFR_MUTATIONS
    scanner = ClinicalEGFRScanner(energy_threshold=-6.0)

    # 1. First scan yields Gatekeeper T790M
    mut1 = scanner.scan(ERLOTINIB, -8.0, active_targets=["WT", "L858R"])
    assert mut1 == "T790M"
    info1 = scanner.get_mutation_info("T790M")
    assert info1 is not None
    assert "Gatekeeper" in info1.name
    assert "Steric Clash" in info1.mechanism

    # 2. Second scan yields C797S
    mut2 = scanner.scan(ERLOTINIB, -8.0, active_targets=["WT", "L858R", "T790M"])
    assert mut2 == "C797S"

    # 3. Third scan yields L718Q
    mut3 = scanner.scan(ERLOTINIB, -8.0, active_targets=["WT", "L858R", "T790M", "C797S"])
    assert mut3 == "L718Q"

    # 4. Weak molecule (> threshold) returns None
    assert scanner.scan(ERLOTINIB, -3.0, active_targets=["WT", "L858R"]) is None


def test_clinical_docking_mechanisms():
    """Verify that clinical mutations correctly apply biophysical mechanisms."""
    docker = RDKitScoreDocker(targets=["WT", "L858R", "T790M", "C797S"])

    osimertinib = "C=CC(=O)Nc1cc(Nc2nccc(-c3cn(C)c4ccccc34)n2)c(OC)cc1N(C)CCN(C)C"
    compact_breaker = "CNS(=O)(=O)c1cccc(-c2ccncc2N2CCOCC2)c1"

    results = docker.dock_many([osimertinib, compact_breaker])
    osim_res, breaker_res = results

    # Osimertinib has covalent acrylamide: binds WT and L858R, but suffers on C797S
    assert osim_res["WT"] < -6.0
    assert osim_res["C797S"] > osim_res["WT"] + 2.0  # significant drop on C797S

    # Compact breaker has sulfonamide/morpholine and no covalent reliance:
    # binds stably across both WT and T790M
    assert breaker_res["WT"] < -7.5
    assert breaker_res["T790M"] < -7.5
