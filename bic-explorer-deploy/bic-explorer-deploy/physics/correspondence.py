"""Coordinate transforms and local boundary samples, not a scattering solver."""
import numpy as np
from . import modes, cylinder


def to_cylindrical(values,phi):
    c,s=np.cos(phi),np.sin(phi)
    result={}
    for group in ['E','H']:
        result[group+'ρ']=values[group+'x']*c+values[group+'y']*s
        result[group+'φ']=-values[group+'x']*s+values[group+'y']*c
        result[group+'z']=values[group+'z']
    return result


def guide_reference(family,p,q,f,w,h,eps,length,x,y,z,amplitude):
    """Uniform guide basis, origin at the left end x=-length/2."""
    return modes.fields(family,p,q,f,w,h,np.asarray(x)+length/2,y,z,eps,amplitude)


def boundary_samples(cylinder_args,w,guide_eps,length,family,p,q,guide_amplitude,cylinder_amplitude,phi):
    f,eps,a,b,h,m,axial,region,angular=cylinder_args
    if w<=2*b or length<=2*b: raise ValueError('Guide width and shown length must exceed the cylinder diameter.')
    z=h*(.25 if region=='Lower disk' else .75)
    # One-sided interior trace avoids floating-point points just outside rho=b.
    r=b*(1-1e-12)
    local,_=cylinder.fields(*cylinder_args,r*np.cos(phi),r*np.sin(phi),z,cylinder_amplitude)
    guide=guide_reference(family,p,q,f,w,h,guide_eps,length,b*np.cos(phi),b*np.sin(phi),z,guide_amplitude)
    return local,guide,to_cylindrical(guide,phi),z
