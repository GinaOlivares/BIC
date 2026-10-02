import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
from physics.correspondence import boundary_samples,guide_reference,to_cylindrical
from physics.modes import propagation,displayed


def render(cylinder_args,amplitude,base_figure,representation,phase):
    f,eps,a,b,h,m,axial,region,angular=cylinder_args
    st.header('Waveguide–cylinder component correspondence')
    st.write('The cylinder axis is z; waveguide propagation is x. Cylinder TM-z means Hz = 0, whereas guide TE-x means Ex = 0. These describe different axes, so a TM-like cylinder can couple to a TE waveguide mode.')
    c1,c2,c3=st.columns(3)
    w=c1.number_input('Surrounding guide width (mm)',min_value=.1,max_value=1000.,value=22.86)*.001
    length=c2.number_input('Shown guide length (mm)',min_value=.1,max_value=2000.,value=60.)*.001
    guide_eps=c3.number_input('Guide background εr',min_value=.01,max_value=100.,value=1.)
    c1,c2,c3,c4=st.columns(4)
    family=c1.selectbox('Reference guide family',['TE','TM'])
    p=c2.number_input('Guide index p',min_value=0,max_value=8,value=1)
    q=c3.number_input('Guide index q',min_value=0,max_value=8,value=0 if family=='TE' else 1,key='guide_q_'+family)
    ga=c4.number_input('Guide H₀ (A/m)' if family=='TE' else 'Guide E₀ (V/m)',min_value=.000001,max_value=1000.,value=.002 if family=='TE' else 1.,format='%.6g',key='guide_amp_'+family)
    phi=np.linspace(0,2*np.pi,361)
    try:
        fc,beta,status=propagation(family,p,q,f,w,h,guide_eps)
        local,guide,gc,z=boundary_samples(cylinder_args,w,guide_eps,length,family,p,q,ga,amplitude,phi)
    except ValueError as exc: st.error(str(exc));return
    st.caption(f'{family}({p},{q}): {status.lower()} · cutoff {fc/1e9:.5g} GHz · same height h = {h*1000:.4g} mm as the inclusion. Reference guide phase/amplitude is specified at x = −{length*500:.4g} mm.')
    st.info('The guide field below is an unperturbed reference mode. The cylinder field is an independent local basis. Their amplitudes and phases are not matched; no transmission, reflection, or numerical coupling coefficient is claimed.')
    st.markdown('**Exact coordinate correspondence at angle φ**')
    st.table({'Cartesian component':['Ex','Ey','Ez','Hx','Hy','Hz'],
              'Cylinder expression':['Eρ cosφ − Eφ sinφ','Eρ sinφ + Eφ cosφ','Ez','Hρ cosφ − Hφ sinφ','Hρ sinφ + Hφ cosφ','Hz']})
    st.latex(r'F_\rho=F_x\cos\phi+F_y\sin\phi,\qquad F_\phi=-F_x\sin\phi+F_y\cos\phi\quad(F=E,H)')
    angle=st.slider('Probe angle φ (degrees)',0,360,45)
    idx=angle
    st.caption(f'Probe on ρ=b at z={z*1000:.4g} mm: x={b*np.cos(np.deg2rad(angle))*1000:.4g} mm, y={b*np.sin(np.deg2rad(angle))*1000:.4g} mm. Values below are complex SI phasors, not normalized display values.')
    def fmt(value): return f'{value.real:.4g} {value.imag:+.4g}i'
    keys=['Eφ','Ez','Hφ','Hz']
    st.dataframe({'Tangential component':keys,'Unit':['V/m','V/m','A/m','A/m'],
        'Cylinder basis':[fmt(local[k][idx]) for k in keys],
        'Unperturbed guide reference':[fmt(gc[k][idx]) for k in keys]},hide_index=True,width='stretch')
    st.write('At the actual dielectric boundary, Eφ, Ez, Hφ and Hz must be continuous. Normal components obey εin Eρ,in = εout Eρ,out and μin Hρ,in = μout Hρ,out in the absence of free surface charge/current. It is the total exterior field (incident plus scattered) that must match—not the reference mode alone.')
    shape=st.toggle('Compare boundary shapes with independent E/H normalization',True)
    fig=make_subplots(rows=2,cols=2,subplot_titles=keys,horizontal_spacing=.14,vertical_spacing=.22)
    for source,values,color in [('Cylinder basis',local,'#b75c35'),('Guide reference',gc,'#287b9b')]:
        for i,k in enumerate(keys):
            scale=max(max(float(np.max(np.abs(values[c]))) for c in [k[0]+'φ',k[0]+'z']),1e-30) if shape else 1.
            fig.add_trace(go.Scatter(x=np.rad2deg(phi),y=displayed(values[k],representation,phase)/scale,name=source,
                legendgroup=source,showlegend=i==0,line=dict(color=color)),row=i//2+1,col=i%2+1)
            fig.update_xaxes(title_text='φ (degrees)',range=[0,360],row=i//2+1,col=i%2+1)
            fig.update_yaxes(title_text='Shape (own E/H scale)' if shape else ('V/m' if i<2 else 'A/m'),row=i//2+1,col=i%2+1)
    fig.update_layout(height=580,template='plotly_white',legend=dict(orientation='h'),margin=dict(l=45,r=30,t=45,b=40))
    st.plotly_chart(fig,width='stretch')
    st.caption('Boundary traces use the same field representation and time phase selected above. Independent normalization compares shape only; disable it to compare the arbitrary SI amplitudes. A visual mismatch is not a computed scattering residual.')
    st.subheader('Inclusion inside the waveguide')
    layer=st.radio('Field shown in shared geometry',['Guide reference outside cylinder','Cylinder basis slice'],horizontal=True)
    component=st.selectbox('Guide reference component',['Ex','Ey','Ez','Hx','Hy','Hz'],index=2,disabled=layer=='Cylinder basis slice')
    fig=go.Figure(base_figure)
    # Keep cylinder geometry; remove its existing field slice for a guide-only reference view.
    if layer=='Guide reference outside cylinder':
        fig.data=fig.data[:-1]
        x=np.linspace(-length/2,length/2,181);y=np.linspace(-w/2,w/2,121)
        xx,yy=np.meshgrid(x,y);zz=np.full_like(xx,z)
        ref=guide_reference(family,p,q,f,w,h,guide_eps,length,xx,yy,zz,ga)
        outside=np.hypot(xx,yy)>b
        data=displayed(ref[component],representation,phase)
        group=component[0]
        scale=max(max(float(np.max(np.abs(ref[k]))) for k in [group+'x',group+'y',group+'z']),1e-30)
        fig.add_trace(go.Surface(x=np.where(outside,xx*1000,np.nan),y=yy*1000,z=zz*1000,
            surfacecolor=data/scale,cmin=0 if representation=='Magnitude |F|' else -1,cmax=1,
            colorscale='Viridis' if representation=='Magnitude |F|' else 'RdBu_r',colorbar=dict(title=component+' / ref',len=.7),
            lighting=dict(ambient=1,diffuse=0,specular=0),hovertemplate='x=%{x:.3g} mm<br>y=%{y:.3g} mm<br>field=%{surfacecolor:.4g}<extra></extra>'))
    else:
        st.caption('Cylinder layer uses the component and slice chosen in the cylinder plot above. The guide component selector is disabled for the cylinder layer.')
    bounds=[(-length/2,length/2),(-w/2,w/2),(0,h)]
    for axis in range(3):
        other=[i for i in range(3) if i!=axis]
        for v1 in bounds[other[0]]:
            for v2 in bounds[other[1]]:
                points=np.zeros((2,3));points[:,axis]=bounds[axis];points[:,other[0]]=v1;points[:,other[1]]=v2
                fig.add_trace(go.Scatter3d(x=points[:,0]*1000,y=points[:,1]*1000,z=points[:,2]*1000,mode='lines',line=dict(color='#426179',width=3),showlegend=False,hoverinfo='skip'))
    fig.update_layout(height=580,margin=dict(l=30,r=85,t=35,b=40),uirevision='shared-guide',scene=dict(aspectmode='data',camera=dict(eye=dict(x=1.8,y=1.8,z=1.6))))
    st.plotly_chart(fig,width='stretch')
    st.caption('PEC guide outline: y=±w/2 and z=0,h; the x ends are reference cross-sections, not metallic terminations. Field layers are separate, never added or presented as a continuous matched field.')
    with st.expander('How these components enter the coupling',expanded=True):
        st.latex(r'W_{n\alpha}\propto\int_{\Sigma_\alpha}(\mathbf E_n\times\mathbf H_\alpha^*+\mathbf E_\alpha^*\times\mathbf H_n)\cdot d\mathbf S')
        st.write('For a waveguide port plane normal to x, the overlap contains Ey Hαz* − Ez Hαy* + Eαy* Hz − Eαz* Hy. This requires the resonant field on that port, which the local cylinder basis does not provide.')
        st.latex(r'W_{n,10}\propto-\int [E_{n,z}H_{10,y}^*+E_{10,z}^*H_{n,y}]\,dy\,dz,\quad H_{n,y}=H_{n,\rho}\sin\phi+H_{n,\phi}\cos\phi')
        st.write('For TE10 the key pairs are cylinder Ez with guide Hy, and cylinder Hy with guide Ez. Coupling is a field-overlap integral, not a one-to-one identification of mode names.')
        st.write('Under the centered mirror-symmetric reference assumptions, sine angular families are forbidden from coupling to TE10 by y-parity; cosine families are allowed but can still have zero overlap. Moving the cylinder or mixing angular families can change that selection rule.')
        st.caption('Calculation_II.pdf: coordinate unit vectors pp.23–24; dielectric matching p.29; channel overlaps pp.30–32. To obtain actual coupling, solve the interior/exterior and upper-port matching, then normalize resonance energy and channel power. No calibrated W is inferred from these independent basis plots.')
