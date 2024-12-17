import pytest

# from blueprints.gwas import GWASTaskManager, GWASTaskStatus
import pandas as pd


@pytest.fixture
def mock_variant_level_assoc(mocker):
    return mocker.patch("blueprints.gwas.variant_level_assoc")


# def test_gwas_task_manager_initialization():
#     """Test GWASTaskManager initializes correctly."""
#     manager = GWASTaskManager()
#     assert manager._tasks == {}


# def test_start_task(mock_variant_level_assoc):
#     """Test starting a new GWAS task."""
#     manager = GWASTaskManager()
#     manager.start_task("250.2")

#     task = manager.get_task_status("250.2")
#     assert task is not None
#     assert task.status == "running"
#     assert task.started_at is not None


# def test_task_completion(mock_variant_level_assoc):
#     """Test task completion with successful GWAS."""
#     # Create mock GWAS results
#     mock_results = pd.DataFrame(
#         {"P": [0.01, 0.02], "gene_id": ["GENE1", "GENE2"], "RSID": ["rs1", "rs2"]}
#     )
#     mock_variant_level_assoc.return_value = mock_results

#     manager = GWASTaskManager()
#     manager.start_task("250.2")

#     # Wait for background task to complete
#     import time

#     time.sleep(0.1)

#     task = manager.get_task_status("250.2")
#     assert task.status == "completed"
#     assert task.completed_at is not None
#     assert len(task.associations) == 2


# def test_task_failure(mock_variant_level_assoc):
#     """Test task handling when GWAS fails."""
#     mock_variant_level_assoc.side_effect = Exception("GWAS failed")

#     manager = GWASTaskManager()
#     manager.start_task("250.2")

#     # Wait for background task to complete
#     import time

#     time.sleep(0.1)

#     task = manager.get_task_status("250.2")
#     assert task.status == "failed"
#     assert "GWAS failed" in task.result


# def test_format_associations():
#     """Test the _format_associations static method."""
#     manager = GWASTaskManager()

#     # Test with empty DataFrame
#     assert manager._format_associations(pd.DataFrame()) == []

#     # Test with valid DataFrame
#     df = pd.DataFrame(
#         {
#             "P": [0.01, 0.06],  # Only first row should be included (P < 0.05)
#             "RSID": ["rs1", "rs2"],
#         }
#     )
#     result = manager._format_associations(df)
#     assert len(result) == 1
#     assert "href" in result[0]["RSID"]
