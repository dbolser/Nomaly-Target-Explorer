import pytest
from flask import template_rendered
from contextlib import contextmanager

@contextmanager
def captured_templates(app):
    recorded = []
    
    def record(sender, template, context, **extra):
        recorded.append((template, context))
        
    template_rendered.connect(record, app)
    try:
        yield recorded
    finally:
        template_rendered.disconnect(record, app)

def test_variant_scores_template(client):
    with captured_templates(client.application) as templates:
        response = client.get('/variant_scores/290.11/GO:0016861')
        assert response.status_code == 200
        
        # Check that the correct template was rendered
        template, context = templates[0]
        assert template.name == 'variant_scores.html'
        
        # Check that required context variables are present
        assert 'data' in context
        assert 'disease_code' in context
        assert 'term' in context
        
        # Check that initial stats values are set
        data = context['data']
        assert 'metric1_pvalue' in data
        assert 'metric1_tpr' in data
        assert 'metric1_fpr' in data
        assert 'metric1_lrp' in data
        assert 'metric1_tp' in data
