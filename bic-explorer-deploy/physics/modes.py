"""Uniform rectangular PEC guide, Calculation_II.pdf pp. 17–21.

Coordinates x (propagation), -w/2 <= y <= w/2, 0 <= z <= h.
Phasors use exp(i beta x), physical fields Re[F exp(-i omega t)].
No resonator, termination or modal superposition is included.
"""
import numpy as np
from scipy.constants import mu_0, epsilon_0
from .validation import positive, integer
from .waveguide import cutoff

PRESETS = {
    'TE': [(1,0),(2,0),(0,1),(1,1),(3,0),(2,1),(0,2),(1,2),(3,1),(2,2)],
    'TM': [(1,1),(2,1),(1,2),(2,2),(3,1),(1,3),(3,2),(2,3),(3,3),(4,1)],
}
COMPONENTS = ('Ex','Ey','Ez','Hx','Hy','Hz')


def validate_mode(family, m, n):
    if family not in PRESETS:
        raise ValueError('Mode family must be TE or TM.')
    integer('m', m); integer('n', n)
    if (family == 'TE' and m == n == 0) or (family == 'TM' and (m == 0 or n == 0)):
        raise ValueError('TE requires (m,n) ≠ (0,0); TM requires m ≥ 1 and n ≥ 1.')


def propagation(family, m, n, frequency, w, h, eps_r=1.):
    validate_mode(family,m,n); positive('frequency',frequency)
    fc = float(cutoff(w,h,eps_r,m,n))
    kc = np.hypot(m*np.pi/w,n*np.pi/h)
    ratio = frequency/fc
    beta = complex(kc*np.lib.scimath.sqrt((ratio-1)*(ratio+1)))
    state = 'PROPAGATING' if frequency > fc else 'EVANESCENT' if frequency < fc else 'AT CUTOFF'
    return fc, beta, state


def fields(family, m, n, frequency, w, h, x, y, z, eps_r=1., amplitude=1.):
    """Return six complex SI phasors. TE amplitude is H0 [A/m]; TM is E0 [V/m].

    x >= 0 selects the passive, decaying evanescent solution, beta=+i alpha.
    """
    _, beta, _ = propagation(family,m,n,frequency,w,h,eps_r)
    positive('amplitude',amplitude)
    x,y,z = np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float),np.asarray(z,float))
    if not all(np.all(np.isfinite(v)) for v in (x,y,z)):
        raise ValueError('Coordinates must be finite.')
    if np.any(x<0) or np.any(np.abs(y)>w/2+1e-14) or np.any(z<0) or np.any(z>h+1e-14):
        raise ValueError('Coordinates require x ≥ 0, −w/2 ≤ y ≤ w/2, 0 ≤ z ≤ h.')
    ky,kz=m*np.pi/w,n*np.pi/h
    kc2=ky*ky+kz*kz
    omega=2*np.pi*frequency
    cy,sy=np.cos(ky*(y+w/2)),np.sin(ky*(y+w/2))
    cz,sz=np.cos(kz*z),np.sin(kz*z)
    p=amplitude*np.exp(1j*beta*x)
    zero=np.zeros(x.shape,dtype=complex)
    if family=='TE':
        longitudinal=cy*cz*p
        dy=-ky*sy*cz*p; dz=-kz*cy*sz*p
        return dict(zip(COMPONENTS,(zero,1j*omega*mu_0/kc2*dz,
             -1j*omega*mu_0/kc2*dy,longitudinal,1j*beta/kc2*dy,1j*beta/kc2*dz)))
    longitudinal=sy*sz*p
    dy=ky*cy*sz*p; dz=kz*sy*cz*p
    return dict(zip(COMPONENTS,(longitudinal,1j*beta/kc2*dy,1j*beta/kc2*dz,
         zero,-1j*omega*epsilon_0*eps_r/kc2*dz,1j*omega*epsilon_0*eps_r/kc2*dy)))


def reference_scales(family,m,n,frequency,w,h,eps_r=1.,amplitude=1.):
    """Shared E/H scales: largest component phasor peak over full section at x=0.

    These analytic maxima are independent of plane, time phase and grid sampling.
    """
    _,beta,_=propagation(family,m,n,frequency,w,h,eps_r)
    positive('amplitude',amplitude)
    ky,kz=m*np.pi/w,n*np.pi/h
    kc2=ky*ky+kz*kz; omega=2*np.pi*frequency
    if family=='TE':
        return (amplitude*omega*mu_0*max(ky,kz)/kc2,
                amplitude*max(1,abs(beta)*ky/kc2,abs(beta)*kz/kc2))
    return (amplitude*max(1,abs(beta)*ky/kc2,abs(beta)*kz/kc2),
            amplitude*omega*epsilon_0*eps_r*max(ky,kz)/kc2)


def plane_grid(plane,w,h,length,slice_fraction,resolution=121):
    for name,val in [('w',w),('h',h),('length',length)]: positive(name,val)
    if not 0<=slice_fraction<=1: raise ValueError('Slice fraction must be in [0,1].')
    integer('resolution',resolution)
    if resolution<3: raise ValueError('Resolution must be at least 3.')
    coords={'x':np.linspace(0,length,resolution), 'y':np.linspace(-w/2,w/2,resolution),
            'z':np.linspace(0,h,resolution)}
    axes={'y–z':('y','z','x'), 'x–y':('x','y','z'), 'x–z':('x','z','y')}
    if plane not in axes: raise ValueError('Unknown plane.')
    u,v,fixed=axes[plane]; U,V=np.meshgrid(coords[u],coords[v])
    fixed_value=coords[fixed][0]+slice_fraction*(coords[fixed][-1]-coords[fixed][0])
    xyz={u:U,v:V,fixed:np.full_like(U,fixed_value)}
    return coords[u],coords[v],xyz,(u,v,fixed,fixed_value)


def displayed(field,representation,phase_degrees=0.):
    if representation=='Magnitude |F|': return np.abs(field)
    if representation=='Instantaneous Re(F e⁻ⁱωᵗ)':
        return np.real(field*np.exp(-1j*np.deg2rad(phase_degrees)))
    if representation=='Imaginary Im(F)': return np.imag(field)
    raise ValueError('Unknown field representation.')
