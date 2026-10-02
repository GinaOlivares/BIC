import numpy as np
import pytest
from scipy.constants import mu_0, epsilon_0
from physics.modes import PRESETS, COMPONENTS, fields, propagation, reference_scales, plane_grid

W,H=.02286,.01016
ALL=[(fam,m,n) for fam, modes in PRESETS.items() for m,n in modes]

@pytest.mark.parametrize('fam,m,n',ALL)
def test_all_presets_pec_and_longitudinal(fam,m,n):
    fc,_,_=propagation(fam,m,n,1e10,W,H)
    f=1.3*fc
    for y in [-W/2,W/2]:
        a=fields(fam,m,n,f,W,H,.01,y,np.linspace(0,H,23))
        assert np.max(np.abs(a['Ex']))<1e-10
        assert np.max(np.abs(a['Ez']))<1e-10
    for z in [0,H]:
        a=fields(fam,m,n,f,W,H,.01,np.linspace(-W/2,W/2,23),z)
        assert np.max(np.abs(a['Ex']))<1e-10
        assert np.max(np.abs(a['Ey']))<1e-10
    assert np.all(a['Ex' if fam=='TE' else 'Hx']==0)

@pytest.mark.parametrize('fam,m,n',[('TE',1,0),('TE',0,2),('TE',2,1),('TM',2,1)])
@pytest.mark.parametrize('ratio',[.7,1.,1.4])
def test_maxwell_curls(fam,m,n,ratio):
    fc,_,_=propagation(fam,m,n,1e10,W,H)
    f=ratio*fc; omega=2*np.pi*f
    p=np.array([.003,.002,.004]); step=1e-8
    def evaluate(p):
        a=fields(fam,m,n,f,W,H,*p)
        return np.array([a[k] for k in COMPONENTS])
    a=evaluate(p)
    grad=[]
    for j in range(3):
        offset=np.eye(3)[j]*step
        grad.append((evaluate(p+offset)-evaluate(p-offset))/(2*step))
    grad=np.array(grad)
    curlE=np.array([grad[1,2]-grad[2,1],grad[2,0]-grad[0,2],grad[0,1]-grad[1,0]])
    curlH=np.array([grad[1,5]-grad[2,4],grad[2,3]-grad[0,5],grad[0,4]-grad[1,3]])
    np.testing.assert_allclose(curlE,1j*omega*mu_0*a[3:],rtol=2e-6,atol=1e-6)
    np.testing.assert_allclose(curlH,-1j*omega*epsilon_0*a[:3],rtol=2e-6,atol=1e-6)

@pytest.mark.parametrize('fam,m,n',ALL)
def test_cutoff_evanescence_and_scales(fam,m,n):
    fc,_,_=propagation(fam,m,n,1e10,W,H)
    assert propagation(fam,m,n,fc,W,H)[1]==0
    _,beta,state=propagation(fam,m,n,.8*fc,W,H)
    assert beta.imag>0 and state=='EVANESCENT'
    base=fields(fam,m,n,.8*fc,W,H,0,.001,.003)
    later=fields(fam,m,n,.8*fc,W,H,.003,.001,.003)
    for k in COMPONENTS:
        assert later[k]==pytest.approx(base[k]*np.exp(-beta.imag*.003))
    for ratio in [.8,1.,1.2]:
        _,_,xyz,_=plane_grid('y–z',W,H,.1,0,151)
        a=fields(fam,m,n,ratio*fc,W,H,**xyz)
        es,hs=reference_scales(fam,m,n,ratio*fc,W,H)
        for k in COMPONENTS:
            assert np.max(np.abs(a[k])) <= (es if k[0]=='E' else hs)*(1+1e-12)


def test_te10_known_profile_and_invalid_modes():
    f=1e10; y=np.linspace(-W/2,W/2,21)
    a=fields('TE',1,0,f,W,H,0,y,0)
    np.testing.assert_allclose(a['Ez'],1j*2*np.pi*f*mu_0*W/np.pi*np.cos(np.pi*y/W),atol=1e-12)
    for mode in [('TE',0,0),('TM',0,1),('TM',1,0),('XX',1,1)]:
        with pytest.raises(ValueError): fields(*mode,f,W,H,0,0,0)
    for plane in ['y–z','x–y','x–z']:
        u,v,xyz,axes=plane_grid(plane,W,H,.1,.3,31)
        assert all(val.shape==(31,31) for val in xyz.values())
        assert np.all(xyz[axes[2]]==axes[3])
