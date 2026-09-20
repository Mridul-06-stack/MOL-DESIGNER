"""Tests for moldesigner.docking module."""

import pytest

from moldesigner.docking import RDKitScoreDocker, VinaDocker, load_docker
from moldesigner.scoring import Docker


class TestRDKitScoreDocker:
    def test_dock_many_returns_correct_shape(self):
        docker = RDKitScoreDocker(targets=["WT", "Mut1"])
        smiles = ["CCO", "c1ccccc1", "not_a_smiles"]
        results = docker.dock_many(smiles)
        
        assert len(results) == 3
        
        # Valid SMILES 1
        assert results[0] is not None
        assert "WT" in results[0]
        assert "Mut1" in results[0]
        
        # Valid SMILES 2
        assert results[1] is not None
        assert "WT" in results[1]
        assert "Mut1" in results[1]
        
        # Invalid SMILES 3
        assert results[2] is None

    def test_dock_many_energies_negative(self):
        docker = RDKitScoreDocker(targets=["target1"])
        results = docker.dock_many(["c1ccc2c(c1)cc(=O)oc2"])
        
        assert results[0] is not None
        # Our pseudo-energy should be negative
        assert results[0]["target1"] < 0

    def test_dock_many_deterministic(self):
        docker1 = RDKitScoreDocker(targets=["tgt1", "tgt2"])
        docker2 = RDKitScoreDocker(targets=["tgt1", "tgt2"])
        
        smiles = ["CCO", "c1ccc(-c2ccncc2)cc1"]
        res1 = docker1.dock_many(smiles)
        res2 = docker2.dock_many(smiles)
        
        assert res1 == res2

    def test_docker_protocol(self):
        """Verify RDKitScoreDocker satisfies the Docker protocol."""
        docker = RDKitScoreDocker()
        assert isinstance(docker, Docker)


class TestVinaDockerImport:
    def test_vina_import_error(self):
        """Should raise ImportError on instantiation if vina is not installed."""
        try:
            import vina
            pytest.skip("vina is installed, skipping import error test")
        except ImportError:
            with pytest.raises(ImportError, match="AutoDock Vina python bindings are required"):
                VinaDocker(receptor_pdbqt_paths={"WT": "fake.pdbqt"}, center=(0,0,0))


class TestLoadDocker:
    def test_load_rdkit_score(self):
        docker = load_docker({
            "backend": "rdkit_score",
            "targets": ["A", "B"]
        })
        assert isinstance(docker, RDKitScoreDocker)
        assert docker.targets == ["A", "B"]

    def test_load_unknown_backend(self):
        with pytest.raises(ValueError, match="Unknown docking backend"):
            load_docker({"backend": "magic"})
            
    def test_load_vina_missing_config(self):
        # We expect a ValueError about receptor_paths before the ImportError triggers
        with pytest.raises(ValueError, match="VinaDocker requires receptor_paths"):
             load_docker({"backend": "vina"})
