import pytest
from flask import url_for
import json
from blueprints.gwas import GWASTaskManager, GWASTaskStatus
from datetime import datetime


@pytest.fixture
def mock_gwas_manager(mocker):
    manager = GWASTaskManager()
    # Mock the variant_level_assoc function
    mocker.patch("blueprints.gwas.variant_level_assoc", return_value=None)
    return manager


def test_run_task_endpoint(client, mock_gwas_manager):
    """Test the run-task endpoint starts a GWAS task."""
    response = client.post("/run-task/250.2")
    assert response.status_code == 202
    data = json.loads(response.data)
    assert data["status"] == "Task started"


def test_task_result_endpoint_not_found(client, mock_gwas_manager):
    """Test getting results for non-existent task."""
    response = client.get("/task-result/999.9")
    assert response.status_code == 404
    data = json.loads(response.data)
    assert data["status"] == "not_found"


def test_task_result_endpoint_running(client, mock_gwas_manager):
    """Test getting results for running task."""
    # Start a task
    mock_gwas_manager.start_task("250.2")

    response = client.get("/task-result/250.2")
    assert response.status_code == 200
    data = json.loads(response.data)
    assert data["status"] == "running"
    assert data["result"] is None


def test_task_result_endpoint_completed(client, mock_gwas_manager):
    """Test getting results for completed task."""
    # Manually set a completed task
    mock_gwas_manager._tasks["250.2"] = GWASTaskStatus(
        status="completed",
        result="GWAS identified 10 missense variants in 5 unique genes",
        started_at=datetime.now(),
        completed_at=datetime.now(),
        associations=[{"gene_id": "GENE1", "P": 0.01}],
    )

    response = client.get("/task-result/250.2")
    assert response.status_code == 200
    data = json.loads(response.data)
    assert data["status"] == "completed"
    assert "GWAS identified" in data["result"]
    assert len(data["associations"]) == 1
