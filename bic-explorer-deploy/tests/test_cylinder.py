from pathlib import Path
import numpy as np
import pytest
from scipy.constants import epsilon_0,mu_0
from physics.cylinder import fields,CARTESIAN
from streamlit.testing.v1 import AppTest

@pytest.mark.parametrize('region',['Lower disk','Upper annulus'])
@pytest.mark.parametrize('m,angular',[(0,'cos'),(1,'cos'),(1,'sin'),(3,'cos')])
def test_curls(region,m,angular):
    args=(10e9,10.,.001,.003,.01,m,.4,region,angular)
    p=np.array([.0016,.0004,.002 if region=='Lower disk' else .007])
    def evaluate(p):
        values,_=fields(*args,*p)
        return np.array([values[k] for k in CARTESIAN])
    values=evaluate(p);grad=[];delta=1e-8
    for axis in range(3):
        off=np.eye(3)[axis]*delta
        grad.append((evaluate(p+off)-evaluate(p-off))/(2*delta))
    g=np.array(grad)
    ce=np.array([g[1,2]-g[2,1],g[2,0]-g[0,2],g[0,1]-g[1,0]])
    ch=np.array([g[1,5]-g[2,4],g[2,3]-g[0,5],g[0,4]-g[1,3]])
    np.testing.assert_allclose(ce,1j*2*np.pi*10e9*mu_0*values[3:],rtol=1e-5,atol=1e-6)
    np.testing.assert_allclose(ch,-1j*2*np.pi*10e9*epsilon_0*10*values[:3],rtol=1e-5,atol=1e-6)


def test_wall_axis_and_mask():
    for m in range(4):
        args=(10e9,10.,.001,.003,.01,m,.4)
        values,_=fields(*args,'Upper annulus','cos',.001,0,.007)
        assert abs(values['Ez'])<1e-12
        assert abs(values['Eφ'])<1e-12
        values,_=fields(*args,'Lower disk','cos',0,0,.002)
        assert all(np.isfinite(v) for v in values.values())
        near,_=fields(*args,'Lower disk','cos',1e-12,0,.002)
        for k in CARTESIAN: assert values[k]==pytest.approx(near[k],abs=1e-8)
    for region,z,x in [('Lower disk',.007,.002),('Upper annulus',.007,0),('Lower disk',.005,0)]:
        values,mask=fields(10e9,10,.001,.003,.01,1,.4,region,'cos',x,0,z)
        assert not mask and np.isnan(values['Ez'])


def test_ui():
    app=AppTest.from_file(Path(__file__).resolve().parents[1]/'pages/2_Cylindrical_Inclusion.py').run(timeout=30)
    assert not app.exception
    assert len(app.get('plotly_chart'))==4
    next(v for v in app.selectbox if v.label=='Dielectric region').set_value('Upper annulus').run()
    assert not app.exception
    next(v for v in app.radio if v.label=='Components').set_value('Cylindrical')
    next(v for v in app.selectbox if v.label=='Viewing plane').set_value('x–z').run()
    assert not app.exception
