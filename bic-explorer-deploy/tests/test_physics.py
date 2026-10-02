import numpy as np
import pytest
from physics.waveguide import cutoff, beta10, guide_wavelength, single_te10
from physics.coupling import radiation, normalized_coupling
from physics.bic import position, quadratic_rate, pa_positions, log_quality
from physics.cmt import spectra, regime
from physics.validation import geometry
from physics.units import ghz_to_hz, mm_to_m

B = 123.4


def test_cutoff():
    fc = cutoff(.02286, .01016)
    assert beta10(fc,.02286,.01016) == 0
    assert 6.55e9 < fc < 6.56e9
    assert beta10(fc*.9,.02286,.01016).imag > 0
    with pytest.raises(ValueError): guide_wavelength(0)

@pytest.mark.parametrize('m', range(6))
@pytest.mark.parametrize('q', range(4))
def test_bic(m,q):
    lb=position(B,m,q)
    assert radiation(B,lb,m) < 1e-26
    assert abs(normalized_coupling(B,lb,m)) < 1e-13


def test_quadratic():
    lb=position(B,1,1); delta=1e-6/B
    a=radiation(B,lb+delta,1)
    assert radiation(B,lb+2*delta,1)/a == pytest.approx(4,rel=1e-8)
    assert a == pytest.approx(quadratic_rate(B,lb+delta,lb),rel=1e-8)


def test_critical_and_passive():
    assert spectra(0,2,2) == pytest.approx((0,1))
    for e in [0,1e-8,1,100]:
        for i in [0,1e-8,1,100]:
            R,A=spectra(np.linspace(-100,100,1001),e,i)
            assert np.all((R>=-1e-14)&(R<=1+1e-14))
            assert np.all((A>=-1e-14)&(A<=1+1e-14))
    assert spectra(0,0,0) == pytest.approx((1,0))
    assert 'CRITICAL' not in regime(0,0)

@pytest.mark.parametrize('m',[0,1,2,3])
def test_exact_pa(m):
    for L in pa_positions(B,position(B,m,1),.3,exact=True):
        assert radiation(B,L,m) == pytest.approx(.3)
    assert pa_positions(B,1,5) == []
    assert pa_positions(B,1,0) == []


def test_safeguards():
    assert not single_te10(10e9,.02286,.030)
    assert single_te10(10e9,.02286,.01016)
    for args in [(-1,1,1), (1,0,1), (1,1,-1)]:
        with pytest.raises(ValueError): cutoff(*args)
    with pytest.raises(ValueError): geometry(1,1,.3,.2)
    with pytest.raises(ValueError): spectra(0,1,-1)
    with pytest.raises(ValueError): normalized_coupling(1j,1,1)
    logq,mask=log_quality(1e10,[0,1e-300,1])
    assert np.all(np.isfinite(logq)) and mask[0]
    assert ghz_to_hz(1)==1e9 and mm_to_m(1)==.001
