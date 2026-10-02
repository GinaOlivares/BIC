"""Phase 1 presentation; all physical relations live in physics/."""
import numpy as np
import streamlit as st
from physics.units import mm_to_m, ghz_to_hz, hz_to_ghz, m_to_mm, hz_to_omega
from physics.validation import geometry
from physics.waveguide import cutoff, beta10, guide_wavelength, single_te10
from physics.coupling import normalized_coupling, radiation
from physics.bic import position, quadratic_rate, pa_positions, log_quality, BIC_TOL
from physics.cmt import spectra, regime
from visualization.curves import curves, markers
from visualization.geometry import schematic

st.set_page_config(page_title='BIC Explorer',layout='wide')
st.title('BIC Explorer')
st.page_link('pages/1_Waveguide_Modes.py', label='Start here: TE / TM waveguide field visualizer', icon='🔬')
st.caption('Phase 1 · Analytical TE10 / coupled-mode model · Calculation_II.pdf')
with st.sidebar:
    st.header('Model controls')
    w=mm_to_m(st.number_input('Width w (mm)',min_value=1.,max_value=1000.,value=22.86))
    h=mm_to_m(st.number_input('Height h (mm)',min_value=.1,max_value=1000.,value=10.16))
    eps=st.number_input('Guide relative permittivity εr',min_value=.01,max_value=1000.,value=1.)
    f=ghz_to_hz(st.number_input('Assumed resonance f₀ (GHz)',min_value=.001,max_value=1000.,value=10.,format='%.6f'))
    m=st.number_input('Cosine angular number m',min_value=0,max_value=30,value=2,step=1)
    q=st.number_input('BIC order q',min_value=0,max_value=20,value=1,step=1)
    L=mm_to_m(st.number_input('Wall distance L (mm)',min_value=0.,max_value=1e6,value=20.,format='%.6f',key='length'))
    gi=st.number_input('Internal decay γi (s⁻¹)',min_value=0.,max_value=1e12,value=.1,format='%.6g')
    gp=st.number_input('Reference decay γ+ (s⁻¹)',min_value=1e-9,max_value=1e12,value=1.,format='%.6g')
    st.caption('Default γ+ = 1 s⁻¹ is arbitrary and uncalibrated. |W̃|² remains dimensionless.')
    with st.expander('Schematic dimensions'):
        b=mm_to_m(st.number_input('Outer radius b (mm)',min_value=.001,value=3.))
        a=mm_to_m(st.number_input('Upper core radius a (mm)',min_value=.0001,value=1.))
try:
    geometry(w,h,a,b)
except ValueError as e:
    st.error(str(e)); st.stop()
fc=cutoff(w,h,eps); fc20=cutoff(w,h,eps,2,0); fc01=cutoff(w,h,eps,0,1)
c1,c2,c3=st.columns(3)
c1.metric('TE10 cutoff',f'{hz_to_ghz(fc):.6f} GHz')
c2.metric('TE20 cutoff',f'{hz_to_ghz(fc20):.6f} GHz')
c3.metric('TE01 cutoff',f'{hz_to_ghz(fc01):.6f} GHz')
st.plotly_chart(schematic(w,h,a,b,L),width='stretch')
st.caption('Schematic — not a full-wave field solution. Outer dielectric: 0 ≤ z ≤ h; gray PEC core: h/2 ≤ z ≤ h. L is measured from the cylinder center to x = −L.')
if L <= b:
    st.error('Wall intersects or touches the resonator: require L > b. Mathematical zero positions may lie outside this physical geometry.')
if f <= fc:
    if f < fc:
        st.error('TE10 is below cutoff. No propagating TE10 channel exists.')
        st.write(f'β10 = i {float(beta10(f,w,h,eps).imag):.6g} rad/m (evanescent).')
    else:
        st.warning('TE10 is at cutoff: β10 = 0. No propagating channel; λg is undefined.')
    st.info('Wall-controlled BIC, Q and one-port propagation plots require f₀ > fc10.')
    st.stop()
beta=float(np.real(beta10(f,w,h,eps))); lg=guide_wavelength(beta)
lb=float(position(beta,m,q)); omega=float(hz_to_omega(f))
we=normalized_coupling(beta,L,m); norm=float(abs(we)**2); ge=float(radiation(beta,L,m,gp))
logq,zero=log_quality(omega,ge,gp)
single=single_te10(f,w,h,eps)
if single: st.success('Single TE10 propagation regime — other rectangular TE/TM channels are closed.')
else: st.warning('Additional rectangular channels may propagate or be at cutoff. TE10 cancellation alone cannot establish a BIC; spectra below are a conditional one-port model.')
if m==0: st.warning('m = 0 generally couples to the TEM port. A TE10 radiation zero is not a true BIC of the complete structure.')
st.caption('A true BIC additionally requires the localized eigenfrequency condition and zero coupling to every open port. f₀ is assumed, not solved. Pure m ≥ 1 suppresses ideal TEM coupling only under the reference symmetry assumptions.')
if lg > 10: st.warning('Near cutoff: guide wavelength exceeds 10 m. Wall phase and positions become extremely sensitive to frequency.')
c1,c2,c3=st.columns(3)
c1.metric('β10',f'{beta:.6g} rad/m'); c2.metric('λg',f'{m_to_mm(lg):.6g} mm'); c3.metric('Selected BIC candidate LB',f'{m_to_mm(lb):.6g} mm')
state=('IDEAL BIC CANDIDATE (lossless model)' if gi==0 else 'RADIATION ZERO / internally lossy') if zero else ('NEAR BIC / quasi-BIC candidate' if norm<1e-3 else 'RADIATIVE RESONANCE')
st.subheader('Physics status')
st.code(f'TE10: OPEN | TE20: {"OPEN" if f>fc20 else "AT CUTOFF" if f==fc20 else "CLOSED"}\nm = {m} | parity = {"even" if m%2==0 else "odd"} | px = {(-1)**m:+d}\n2βL = {2*beta*L:.9g} rad\nNormalized |W̃eff|² = {norm:.9g}\nγe = {ge:.9g} s⁻¹ | γi = {gi:.9g} s⁻¹\n{state}\n{regime(ge,gi)}')
st.write('Qrad: BIC threshold reached; ideal limit diverges.' if zero else f'log₁₀ Qrad = {float(logq):.6g} (display cap 16).')
orders=list(range(max(5,q+2)))
zeros=[float(position(beta,m,k)) for k in orders]
st.dataframe({'q':orders,'LB (mm)':[float(m_to_mm(x)) for x in zeros],'Clears cylinder':[x>b for x in zeros]},hide_index=True)
end=max(2.,lb/lg+.5,float(L/lg)+.2)
u=np.unique(np.r_[np.linspace(0,end,1801),np.asarray(zeros)/lg,lb/lg])
lengths=u*lg; rates=radiation(beta,lengths,m,gp)
fig=curves(u,{'Exact normalized coupling':rates/gp},'L / λg','|W̃eff|² = γe / γ+')
near=np.abs(lengths-lb)*beta <= .35
fig.add_scatter(x=u[near],y=quadratic_rate(beta,lengths[near],lb,gp)/gp,name='Near-BIC quadratic',mode='lines',line_dash='dash')
markers(fig,[x/lg for x in zeros if x/lg<=end]); markers(fig,[x/lg for x in pa_positions(beta,lb,gi,gp,exact=True) if x>b], 'Critical coupling', '#2a9d8f'); fig.add_vline(x=float(L/lg),line_color='#e07a3f')
st.plotly_chart(fig,width='stretch')
st.caption('Even/odd zeros shift by λg/4; each family repeats every λg/2. Orange line: selected L. Quadratic curve shown only for |β(L−LB)| ≤ 0.35.')
logs,mask=log_quality(omega,rates,gp)
fig=curves(u,{'log₁₀ Qrad (capped)':logs},'L / λg','log₁₀ Qrad')
markers(fig,[x/lg for x in zeros if x/lg<=end]); st.plotly_chart(fig,width='stretch')
st.caption(f'Normalized |W̃|² ≤ {BIC_TOL:g} marks the numerical BIC threshold. Display cap: log₁₀ Q = 16. Exact ideal Q diverges; Q here includes TE10 radiation only.')
st.subheader('Critical coupling and one-port spectrum')
exact=pa_positions(beta,lb,gi,gp,exact=True); approx=pa_positions(beta,lb,gi,gp)
if exact:
    st.write('Exact positive-loss critical positions near selected LB (mm): '+', '.join(f'{m_to_mm(x):.6g}'+(' (wall overlap)' if x<=b else '') for x in exact))
    st.write('Near-BIC approximation (mm): '+', '.join(f'{m_to_mm(x):.6g}'+(' (wall overlap)' if x<=b else '') for x in approx))
    if gi/gp>.1: st.caption('γi/γ+ is not small; the near-BIC approximation may be inaccurate. Use exact positions.')
else: st.info('No positive-absorption critical position: require 0 < γi ≤ 4γ+.')
span=max(ge+gi,gp)*5
detuning=np.linspace(-span,span,1001)
R,A=spectra(detuning,ge,gi)
fig=curves(detuning/(2*np.pi),{'Reflectance R':R,'Absorption A':A},'Frequency detuning f − f₀ (Hz)','Power fraction')
fig.add_vline(x=0,line_dash='dot'); fig.update_yaxes(range=[-.02,1.02]); st.plotly_chart(fig,width='stretch')
r0,a0=spectra(0,ge,gi)
st.write(f'At resonance: R = {float(r0):.6g}, A = {float(a0):.6g}. {regime(ge,gi)}.')
st.caption('Narrowband CMT sweep: rates and β fixed at f₀. Detuning avoids precision loss with the arbitrary 1 s⁻¹ reference scale. γ/(2π) converts rate to linewidth in Hz.')
if span/omega>.01: st.warning('Sweep span is large relative to ω₀; the frozen-rate narrowband approximation may not be valid.')
with st.expander('Equations used',expanded=False):
    st.latex(r'f_{c,mn}=\frac{\sqrt{(m/w)^2+(n/h)^2}}{2\sqrt{\mu_0\epsilon_0\epsilon_r}},\quad\beta_{10}=\sqrt{\omega_0^2\mu_0\epsilon_0\epsilon_r-(\pi/w)^2},\quad\lambda_g=2\pi/\beta_{10}')
    st.latex(r'\widetilde W=1-(-1)^m e^{2i\beta L},\quad\gamma_e=\gamma_+|\widetilde W|^2,\quad Q_{rad}=\omega_0/\gamma_e')
    st.latex(r'L_B=\begin{cases}q\lambda_g/2&m\;\mathrm{even}\\(2q+1)\lambda_g/4&m\;\mathrm{odd}\end{cases},\quad\gamma_e\simeq4\gamma_+\beta^2(L-L_B)^2')
    st.latex(r'r=\frac{i\Delta+(\gamma_e-\gamma_i)/2}{-i\Delta+(\gamma_e+\gamma_i)/2},\quad R=|r|^2,\quad A=1-R')
    st.latex(r'L_{PA}^{\pm}\simeq L_B\pm\frac{\sqrt{\gamma_i/\gamma_+}}{2\beta},\quad L_{PA,exact}^{\pm}=L_B\pm\frac{\arcsin(\sqrt{\gamma_i/\gamma_+}/2)}{\beta}')
    st.caption('Calculation_II.pdf: pp. 19–20, 68–80, 100–113. Geometry is a schematic only; it does not calculate resonator fields.')
with st.expander('Explain the physics',expanded=True):
    st.write(f'For m = {m}, px = {(-1)**m:+d}. The direct and PEC-reflected amplitudes combine as W̃ = 1 − px exp(2iβL). At LB = {m_to_mm(lb):.6g} mm, their phases cancel the TE10 radiation. At your selected L, |W̃|² = {norm:.6g}. Positive equal external and internal rates give perfect absorption at resonance; zero external coupling cannot absorb incident power.')
    st.write('Unresolved: actual eigenfrequency and wall-induced shift, absolute coupling calibration, full mode matching, angular mixing and additional-port losses. Radii affect the schematic only. This is an analytical/CMT explorer, not a full-wave Maxwell solver.')
