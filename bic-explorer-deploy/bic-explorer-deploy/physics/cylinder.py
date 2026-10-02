"""Local TM-z basis only; Calculation_II.pdf §§12–20, pp.24–29.

No matching between regions or to exterior ports is solved.
"""
import numpy as np
from scipy.constants import epsilon_0, mu_0
from scipy.special import jv, yv, jvp, yvp
from .validation import positive, integer

CARTESIAN=('Ex','Ey','Ez','Hx','Hy','Hz')
CYLINDRICAL=('Eρ','Eφ','Ez','Hρ','Hφ','Hz')


def parameters(frequency,eps_r,a,b,h,m,axial_fraction):
    for name,val in [('frequency',frequency),('eps_r',eps_r),('a',a),('b',b),('h',h)]: positive(name,val)
    integer('m',m)
    if a>=b: raise ValueError('Require 0 < a < b.')
    if not 0<=axial_fraction<1: raise ValueError('Require 0 ≤ kz/kd < 1; kρ=0 is excluded.')
    kd=2*np.pi*frequency*np.sqrt(mu_0*epsilon_0*eps_r)
    kz=axial_fraction*kd
    return kz,np.sqrt(kd*kd-kz*kz)


def radial(r,m,kr,a,b,region):
    if region=='Lower disk':
        return jv(m,kr*r),kr*jvp(m,kr*r)
    if region!='Upper annulus': raise ValueError('Unknown region.')
    # Scale the Bessel combination without changing the PEC zero at rho=a.
    ca,cb=yv(m,kr*a),jv(m,kr*a)
    scale=max(abs(ca),abs(cb))
    ca,cb=ca/scale,cb/scale
    return ca*jv(m,kr*r)-cb*yv(m,kr*r),kr*(ca*jvp(m,kr*r)-cb*yvp(m,kr*r))


def fields(frequency,eps_r,a,b,h,m,axial_fraction,region,angular,x,y,z,amplitude=1.):
    kz,kr=parameters(frequency,eps_r,a,b,h,m,axial_fraction)
    positive('amplitude',amplitude)
    if angular not in ('cos','sin') or (angular=='sin' and m==0):
        raise ValueError('Use cosine for m=0; sine is identically zero.')
    x,y,z=np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float),np.asarray(z,float))
    if not all(np.all(np.isfinite(v)) for v in [x,y,z]): raise ValueError('Coordinates must be finite.')
    rho=np.hypot(x,y); phi=np.arctan2(y,x)
    if region=='Lower disk': mask=(rho<=b)&(z>=0)&(z<h/2)
    elif region=='Upper annulus': mask=(rho>=a)&(rho<=b)&(z>h/2)&(z<=h)
    else: raise ValueError('Unknown region.')
    safe=np.clip(rho,0 if region=='Lower disk' else a,b)
    R,dR=radial(safe,m,kr,a,b,region)
    rr=np.linspace(0 if region=='Lower disk' else a,b,1001)
    norm=np.max(np.abs(radial(rr,m,kr,a,b,region)[0]))
    if not np.isfinite(norm) or norm<=0: raise ValueError('Basis normalization is unresolved for these parameters.')
    R,dR=R/norm,dR/norm
    angle=np.cos(m*phi) if angular=='cos' else np.sin(m*phi)
    dangle=-m*np.sin(m*phi) if angular=='cos' else m*np.cos(m*phi)
    ratio=np.divide(R,safe,out=np.zeros_like(R),where=safe!=0)
    if region=='Lower disk' and m==1: ratio=np.where(safe==0,kr/(2*norm),ratio)
    p=amplitude*np.exp(1j*kz*(z-(0 if region=='Lower disk' else h/2)))
    ez=R*angle*p; dr=dR*angle*p; dp=ratio*dangle*p
    er=1j*kz/kr**2*dr; ep=1j*kz/kr**2*dp
    hr=-1j*2*np.pi*frequency*epsilon_0*eps_r/kr**2*dp
    hp=1j*2*np.pi*frequency*epsilon_0*eps_r/kr**2*dr
    c,s=np.cos(phi),np.sin(phi)
    result={'Ex':er*c-ep*s,'Ey':er*s+ep*c,'Ez':ez,
            'Hx':hr*c-hp*s,'Hy':hr*s+hp*c,'Hz':np.zeros_like(ez),
            'Eρ':er,'Eφ':ep,'Hρ':hr,'Hφ':hp}
    return {key:np.where(mask,value,np.nan+0j) for key,value in result.items()},mask


def slice_grid(plane,b,h,fraction,resolution=121):
    if not 0<=fraction<=1: raise ValueError('Slice fraction must be in [0,1].')
    axes={'x–y':('x','y','z'),'x–z':('x','z','y'),'y–z':('y','z','x')}
    u,v,fixed=axes[plane]
    coords={'x':np.linspace(-b,b,resolution),'y':np.linspace(-b,b,resolution),'z':np.linspace(0,h,resolution)}
    U,V=np.meshgrid(coords[u],coords[v]); pos=coords[fixed][0]+fraction*(coords[fixed][-1]-coords[fixed][0])
    return coords[u],coords[v],{u:U,v:V,fixed:np.full_like(U,pos)},(u,v,fixed,pos)
