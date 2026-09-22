"""Tests for FastAPI endpoints and CORS security configuration."""

import pytest
from fastapi.testclient import TestClient
from server.api import app

client = TestClient(app)


def test_cors_allowed_origin():
    """Allowed origin should receive Access-Control-Allow-Origin."""
    response = client.options(
        "/api/targets",
        headers={
            "Origin": "https://alchemist-ai-8lhb.onrender.com",
            "Access-Control-Request-Method": "GET",
        },
    )
    assert response.headers.get("access-control-allow-origin") == "https://alchemist-ai-8lhb.onrender.com"
    # Ensure credentials are not allowed for public endpoints (RepoSentinel Finding 01)
    assert response.headers.get("access-control-allow-credentials") != "true"


def test_cors_disallowed_origin():
    """Arbitrary untrusted origin should not receive Access-Control-Allow-Origin."""
    response = client.options(
        "/api/targets",
        headers={
            "Origin": "https://evil-attacker.example.com",
            "Access-Control-Request-Method": "GET",
        },
    )
    assert "access-control-allow-origin" not in response.headers


def test_health_check():
    """GET /api/health should return status ok."""
    response = client.get("/api/health")
    assert response.status_code == 200
    data = response.json()
    assert data["status"] == "ok"
    assert data["version"] == "0.1.0"


def test_get_molblock():
    """GET /api/molblock should return valid 3D SDF block."""
    response = client.get("/api/molblock?smiles=CC(=O)Oc1ccccc1C(=O)O")
    assert response.status_code == 200
    assert "M  END" in response.text
