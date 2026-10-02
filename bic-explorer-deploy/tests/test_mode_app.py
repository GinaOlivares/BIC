from pathlib import Path
from streamlit.testing.v1 import AppTest

APP=Path(__file__).resolve().parents[1]/'pages/1_Waveguide_Modes.py'


def test_mode_controls_and_planes():
    app=AppTest.from_file(APP).run(timeout=30)
    assert not app.exception
    assert len(app.get('plotly_chart'))==2
    app.radio[0].set_value('TM').run()
    assert not app.exception
    assert any('Below cutoff' in warning.value for warning in app.warning)
    next(x for x in app.number_input if x.label=='Frequency (GHz)').set_value(40.)
    for plane in ['x–y','x–z','y–z']:
        next(x for x in app.selectbox if x.label=='Viewing plane').set_value(plane)
        next(x for x in app.selectbox if x.label=='Field representation').set_value('Instantaneous Re(F e⁻ⁱωᵗ)')
        app.run(timeout=30)
        assert not app.exception
    next(x for x in app.selectbox if 'Mode preset' in x.label).set_value('Custom indices').run()
    next(x for x in app.number_input if x.label=='Index m').set_value(0).run()
    assert not app.exception
    assert any('TM requires' in err.value for err in app.error)


def test_every_preset_renders():
    from physics.modes import PRESETS
    app=AppTest.from_file(APP).run(timeout=30)
    for family,indices in PRESETS.items():
        app.radio[0].set_value(family).run()
        for m,n in indices:
            next(x for x in app.selectbox if 'Mode preset' in x.label).set_value(f'{family}{m}{n}').run(timeout=30)
            assert not app.exception
            assert len(app.get('plotly_chart'))==2
