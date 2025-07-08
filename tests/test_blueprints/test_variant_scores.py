from flask import template_rendered
from contextlib import contextmanager
import pytest

from tests.data_utils import production_data_available

pytestmark = pytest.mark.skipif(
    not production_data_available(),
    reason="requires production data services",
)


@contextmanager
def captured_templates(integration_app):
    recorded = []

    def record(sender, template, context, **extra):
        recorded.append((template, context))

    template_rendered.connect(record, integration_app)
    try:
        yield recorded
    finally:
        template_rendered.disconnect(record, integration_app)


def test_variant_scores_template(auth_integration_app_client):
    with captured_templates(auth_integration_app_client.application) as templates:
        response = auth_integration_app_client.get("/variant_scores/290.11/GO:0016861")
        assert response.status_code == 200

        # Check that the correct template was rendered
        template, context = templates[0]
        assert template.name == "variant_scores.html"

        # Check that required context variables are present
        assert "data" in context
        assert "disease_code" in context
        assert "term" in context
