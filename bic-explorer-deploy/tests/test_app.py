from pathlib import Path
from streamlit.testing.v1 import AppTest
APP = Path(__file__).resolve().parents[1] / 'app.py'


def test_default_and_cutoff():
    app=AppTest.from_file(APP).run(timeout=30)
    assert not app.exception
    assert len(app.get('plotly_chart')) == 4
    freq=next(x for x in app.number_input if x.label=='Assumed resonance f₀ (GHz)')
    freq.set_value(1.).run()
    assert not app.exception
    assert any('below cutoff' in e.value for e in app.error)


def test_odd_lossless_and_invalid_geometry():
    app=AppTest.from_file(APP).run(timeout=30)
    next(x for x in app.number_input if x.label=='Cosine angular number m').set_value(1)
    next(x for x in app.number_input if x.label=='Internal decay γi (s⁻¹)').set_value(0.)
    app.run()
    assert not app.exception
    next(x for x in app.number_input if x.label=='Upper core radius a (mm)').set_value(4.)
    app.run()
    assert not app.exception
    assert any('a < b' in e.value for e in app.error)
