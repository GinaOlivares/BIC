import numpy as np
from physics.correspondence import to_cylindrical,boundary_samples


def test_rotation_and_norm():
    phi=np.array([0,np.pi/2,np.pi,1.234])
    source={k:np.full(4,v,dtype=complex) for k,v in {'Ex':1+2j,'Ey':3-1j,'Ez':4,'Hx':2j,'Hy':7,'Hz':0}.items()}
    result=to_cylindrical(source,phi)
    for g in ['E','H']:
        np.testing.assert_allclose(result[g+'ρ']*np.cos(phi)-result[g+'φ']*np.sin(phi),source[g+'x'],atol=1e-14)
        np.testing.assert_allclose(np.abs(result[g+'ρ'])**2+np.abs(result[g+'φ'])**2,np.abs(source[g+'x'])**2+np.abs(source[g+'y'])**2)
    assert result['Eρ'][0]==source['Ex'][0]
    np.testing.assert_allclose(result['Eφ'][1],-source['Ex'][1])


def test_te10_boundary_correspondence():
    args=(10e9,10.,.001,.003,.01016,1,.4,'Lower disk','cos')
    local,guide,cyl,z=boundary_samples(args,.02286,1.,.06,'TE',1,0,.002,1.,np.linspace(0,2*np.pi,361))
    assert all(np.all(np.isfinite(v)) for v in local.values())
    assert np.all(cyl['Eφ']==0) and np.all(cyl['Hz']==0)
    np.testing.assert_allclose(cyl['Ez'],guide['Ez'])
