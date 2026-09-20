"""Tests for moldesigner.validation module."""

import pytest

from moldesigner.docking import RDKitScoreDocker
from moldesigner.validation import redock, enrichment_summary


class TestValidation:
    def setup_method(self):
        self.docker = RDKitScoreDocker(targets=["WT"])
        self.smiles = [
            "CCO",
            "c1ccccc1",
            "c1ccc2c(c1)cc(=O)oc2"
        ]
        
    def test_redock_returns_reports(self):
        # First "original" docking
        orig_scores = self.docker.dock_many(self.smiles)
        
        reports = redock(self.docker, self.smiles, orig_scores)
        
        assert len(reports) == 3
        assert reports[0].smiles == "CCO"
        assert "WT" in reports[0].original_scores
        assert "WT" in reports[0].redock_scores
        assert reports[0].consistent is True

    def test_redock_consistency_threshold(self):
        # Fake original scores
        orig_scores = [
            {"WT": -5.0}, # CCO
            {"WT": -5.0}, # c1ccccc1
            {"WT": -5.0}  # Coumarin
        ]
        
        # Redocking should use RDKitScoreDocker which gives scores around -5 to -10
        # which might be more than 2.0 away from our fake scores
        reports = redock(self.docker, self.smiles, orig_scores, consistency_threshold=0.1)
        
        for r in reports:
            diff = abs(r.original_scores["WT"] - r.redock_scores["WT"])
            if diff > 0.1:
                assert r.consistent is False
            else:
                assert r.consistent is True
        
    def test_enrichment_summary_empty(self):
        summary = enrichment_summary([])
        assert summary["total"] == 0
        
    def test_enrichment_summary(self):
        orig_scores = self.docker.dock_many(self.smiles)
        reports = redock(self.docker, self.smiles, orig_scores)
        
        summary = enrichment_summary(reports)
        
        assert summary["total"] == 3
        assert summary["consistent"] == 3
        assert summary["consistency_rate"] == 1.0
        assert summary["mean_deviation"] == 0.0 # RDKitScoreDocker is 100% deterministic
