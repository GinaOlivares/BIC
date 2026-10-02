import numpy as np
import streamlit as st
import plotly.graph_objects as go
from physics.cylinder import fields,parameters,slice_grid,CARTESIAN,CYLINDRICAL
from physics.modes import displayed
from visualization.mode_fields import field_panels
from visualization.correspondence_app import render as render_correspondence


def render():
    st.markdown('''<style>@media(min-width:800px){[data-testid="stSidebar"]{width:240px!important;min-width:240px!important;max-width:240px!important;}[data-testid="stMainBlockContainer"]{padding:2.5rem 2rem;}}h1{font-size:2.1rem!important;}</style>''',unsafe_allow_html=True)
    st.title('Cylindrical Inclusion Fields')
    st.caption('Local TM-z basis · Calculation_II.pdf §§12–20, pp. 24–29 · SI calculations')
    st.warning('Local basis visualization, not a solved resonant field. Lower and upper regions are explored separately; interface matching, exterior waveguide and TEM-port matching are not solved.')
    with st.sidebar:
        st.header('Inclusion controls')
        region=st.selectbox('Dielectric region',['Lower disk','Upper annulus'])
        m=st.number_input('Angular order m',0,8,1)
        angular=st.selectbox('Angular family',['cos'] if m==0 else ['cos','sin'])
        f=st.number_input('Frequency (GHz)',.1,100.,10.)*1e9
        eps=st.number_input('Cylinder εr',1.,100.,10.)
        b=st.number_input('Outer radius b (mm)',.1,100.,3.)*.001
        a=st.number_input('Core radius a (mm)',.01,99.,1.)*.001
        h=st.number_input('Height h (mm)',.1,100.,10.16)*.001
        axial=st.slider('Axial fraction kz / kd',0.,.95,.4,.01)
        amplitude=st.number_input('Longitudinal reference E₀ (V/m)',.001,1000.,1.)
    try: kz,kr=parameters(f,eps,a,b,h,m,axial)
    except ValueError as e: st.error(str(e));return
    st.caption(f'kz = {kz:.5g} rad/m · kρ = {kr:.5g} rad/m · Hz = 0 for this local TM basis. kz is a chosen basis parameter, not an eigenvalue.')
    c1,c2,c3=st.columns([1,1,1.5])
    plane=c1.selectbox('Viewing plane',['x–y','x–z','y–z'])
    fraction=c2.slider('Slice position (%)',0.,100.,25. if region=='Lower disk' else 75.,key='cylinder_slice_'+region)/100
    representation=c3.selectbox('Field representation',['Magnitude |F|','Instantaneous Re(F e⁻ⁱωᵗ)','Imaginary Im(F)'])
    c1,c2,c3=st.columns(3)
    coordinates=c1.radio('Components',['Cartesian','Cylindrical'],horizontal=True)
    phase=c2.slider('Time phase ωt (degrees)',0,360,45,5,disabled=representation!='Instantaneous Re(F e⁻ⁱωᵗ)')
    normalized=c3.toggle('Normalize E and H separately',True)
    names=CARTESIAN if coordinates=='Cartesian' else CYLINDRICAL
    u,v,xyz,axes=slice_grid(plane,b,h,fraction)
    args=(f,eps,a,b,h,m,axial,region,angular)
    values,mask=fields(*args,**xyz,amplitude=amplitude)
    # Reference sampled in full selected cross-section, not on the user's slice.
    radius=np.linspace(0 if region=='Lower disk' else a,b,201)
    phi=np.linspace(0,2*np.pi,193)
    rr,pp=np.meshgrid(radius,phi)
    reference,_=fields(*args,rr*np.cos(pp),rr*np.sin(pp),np.full_like(rr,h*.25 if region=='Lower disk' else h*.75),amplitude=amplitude)
    scales=tuple(max(float(np.nanmax(np.abs(reference[k]))) for k in names if k.startswith(group)) for group in ['E','H'])
    scales=tuple(max(s,1e-30) for s in scales)
    st.markdown(f'**{plane} slice · {axes[2]} = {axes[3]*1000:.4g} mm**')
    st.caption('Blank regions are outside the selected dielectric region, inside PEC, or on the unresolved z=h/2 matching plane. They are not computed zero fields.')
    if not np.any(mask): st.info('This slice does not intersect the selected dielectric region. Move the slice or change region.')
    panels=field_panels(u,v,values,axes,scales,representation,phase,normalized,names)
    for i in range(1,7):
        suffix='' if i==1 else str(i)
        panels.layout['yaxis'+suffix].update(scaleanchor='x'+suffix,scaleratio=1)
        panels.layout['xaxis'+suffix].update(constrain='domain')
    st.plotly_chart(panels,width='stretch')
    st.subheader('Fields inside the cylindrical geometry')
    component=st.selectbox('3D field component',names,index=2)
    fig=go.Figure()
    theta=np.linspace(0,2*np.pi,81)
    for radius,z0,z1,color in [(b,0,h,'#90bed1'),(a,h/2,h,'#666666')]:
        t,z=np.meshgrid(theta,[z0,z1])
        fig.add_trace(go.Surface(x=radius*np.cos(t)*1000,y=radius*np.sin(t)*1000,z=z*1000,
            surfacecolor=np.zeros_like(t),colorscale=[[0,color],[1,color]],opacity=.14 if radius==b else .7,showscale=False,hoverinfo='skip'))
        for level in [z0,z1]:
            fig.add_trace(go.Scatter3d(x=radius*np.cos(theta)*1000,y=radius*np.sin(theta)*1000,z=np.full_like(theta,level*1000),mode='lines',line=dict(color=color,width=3),showlegend=False,hoverinfo='skip'))
    # PEC lower face: a disk, part of the upper half-height solid core.
    t,r=np.meshgrid(theta,np.linspace(0,a,20))
    for level in [h/2,h]:
        fig.add_trace(go.Surface(x=r*np.cos(t)*1000,y=r*np.sin(t)*1000,z=np.full_like(t,level*1000),surfacecolor=np.zeros_like(t),colorscale=[[0,'#777'],[1,'#777']],opacity=.65,showscale=False,hoverinfo='skip'))
    idx=0 if component.startswith('E') else 1
    scale=scales[idx] if normalized else 1.
    limit=1. if normalized else scales[idx]
    data=displayed(values[component],representation,phase)/scale
    unit=('normalized '+component[0]) if normalized else ('V/m' if idx==0 else 'A/m')
    fig.add_trace(go.Surface(x=np.where(mask,xyz['x']*1000,np.nan),y=xyz['y']*1000,z=xyz['z']*1000,
        surfacecolor=data,cmin=0 if representation=='Magnitude |F|' else -limit,cmax=limit,
        colorscale='Viridis' if representation=='Magnitude |F|' else 'RdBu_r',
        colorbar=dict(title=unit,len=.7),lighting=dict(ambient=1,diffuse=0,specular=0),
        hovertemplate=f'x=%{{x:.3g}} mm<br>y=%{{y:.3g}} mm<br>z=%{{z:.3g}} mm<br>{component}=%{{surfacecolor:.4g}}<extra></extra>'))
    fig.update_layout(height=540,margin=dict(l=0,r=50,t=15,b=20),uirevision='cylinder',scene=dict(aspectmode='data',xaxis_title='x (mm)',yaxis_title='y (mm)',zaxis_title='z (mm)'))
    st.plotly_chart(fig,width='stretch')
    st.caption('Blue outline: dielectric outer cylinder over 0≤z≤h. Gray solid core: ρ≤a, h/2≤z≤h. The other dielectric region is uncomputed, not absent. Rotate the geometry to inspect the field slice.')
    render_correspondence(args,amplitude,fig,representation,phase)
    with st.expander('Equations, normalization and unresolved matching',expanded=True):
        st.latex(r'E_z=E_0\widehat R_m(\rho)\Phi_m(\phi)e^{ik_z(z-z_{ref})},\quad H_z=0,\quad k_\rho^2=\omega^2\mu_0\epsilon_0\epsilon_r-k_z^2')
        st.latex(r'R_m^{I}=J_m(k_\rho\rho),\quad R_m^{II}=Y_m(k_\rho a)J_m(k_\rho\rho)-J_m(k_\rho a)Y_m(k_\rho\rho)')
        st.latex(r'E_\rho=\frac{ik_z}{k_\rho^2}\partial_\rho E_z,\quad E_\phi=\frac{ik_z}{k_\rho^2\rho}\partial_\phi E_z,\quad H_\rho=-\frac{i\omega\epsilon}{k_\rho^2\rho}\partial_\phi E_z,\quad H_\phi=\frac{i\omega\epsilon}{k_\rho^2}\partial_\rho E_z')
        st.write('Φ is cos(mφ) or sin(mφ). R̂ has unit sampled peak in the selected radial interval; E₀ is an arbitrary field amplitude, not a power calibration. Optional E/H display scales are shared by their three components and calculated across the selected region. Time convention: exp(−iωt).')
        st.write('The lower basis is regular at the axis. The upper basis enforces Ez=Eφ=0 on the vertical PEC core wall. A single traveling basis term does not enforce the horizontal PEC faces or guide plates. The full field requires a superposition with matched coefficients at z=h/2, ρ=b, and the upper port. No artificial PEC boundary is imposed at the outer dielectric surface.')
        st.caption('Cylindrical components at ρ=0 depend on the angular basis direction; Cartesian fields use their finite analytic limits. This is a local homogeneous Maxwell basis, not the globally matched inclusion or BIC field.')
