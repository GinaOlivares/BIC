import io
import numpy as np
import streamlit as st
from physics.modes import PRESETS, COMPONENTS, fields, propagation, reference_scales, plane_grid
from physics.units import mm_to_m, ghz_to_hz
from visualization.mode_fields import field_panels, locator


def render():
    st.markdown("""<style>
    @media (min-width: 800px) {
      [data-testid="stSidebar"] {width:240px !important; min-width:240px !important; max-width:240px !important;}
      [data-testid="stMainBlockContainer"] {padding:2.4rem 2rem 2rem;}
    }
    [data-testid="stSidebarUserContent"] {padding:1rem 1rem 1.5rem;}
    [data-testid="stSidebar"] [data-testid="stVerticalBlock"] {gap:.75rem;}
    [data-testid="stMetricValue"] {font-size:1.65rem;}
    h1 {font-size:2.1rem !important; padding-bottom:.3rem !important;}
    h3 {font-size:1.35rem !important;}
    </style>""", unsafe_allow_html=True)
    st.title('Waveguide Mode Explorer')
    st.caption('Analytical rectangular PEC guide · Calculation_II.pdf, pp. 17–21')
    st.caption('Propagation: x · Width: y · Height: z · Six analytical electric and magnetic field components')
    with st.sidebar:
        st.header('Waveguide modes')
        family=st.radio('Mode family',['TE','TM'],horizontal=True,key='mode_family')
        options=[f'{family}{m}{n}' for m,n in PRESETS[family]]+['Custom indices']
        preset=st.selectbox('Mode preset — 10 per family',options,key=f'preset_{family}')
        if preset=='Custom indices':
            m=st.number_input('Index m',min_value=0,max_value=10,value=1,step=1,key=f'index_m_{family}')
            n=st.number_input('Index n',min_value=0,max_value=10,value=0 if family=='TE' else 1,step=1,key=f'index_n_{family}')
        else:
            m,n=PRESETS[family][options.index(preset)]
            st.caption(f'm = {m} · n = {n}. Select Custom indices to edit either index.')
        frequency=ghz_to_hz(st.number_input('Frequency (GHz)',min_value=.001,max_value=1000.,value=10.,format='%.6f',key='mode_frequency'))
        w=mm_to_m(st.number_input('Width w (mm)',min_value=1.,max_value=1000.,value=22.86,key='mode_width'))
        h=mm_to_m(st.number_input('Height h (mm)',min_value=.1,max_value=1000.,value=10.16,key='mode_height'))
        eps=st.number_input('Relative permittivity εr',min_value=.01,max_value=1000.,value=1.,key='mode_eps')
        length=mm_to_m(st.number_input('Shown x length (mm)',min_value=.1,max_value=10000.,value=60.,key='mode_length'))
        amplitude=st.number_input('H₀ (A/m)' if family=='TE' else 'E₀ (V/m)',min_value=.000001,max_value=1e6,value=1.,format='%.6g',key=f'amp_{family}')
    try:
        fc,beta,state=propagation(family,m,n,float(frequency),float(w),float(h),eps)
    except ValueError as exc:
        st.error(str(exc)); return
    st.subheader(f'{family}({m}, {n}) · {state.lower()}')
    columns=st.columns(3)
    columns[0].metric('Cutoff frequency',f'{fc/1e9:.6f} GHz')
    columns[1].metric('β (rad/m)' if beta.imag==0 else 'Decay α (m⁻¹)',f'{beta.real if beta.imag==0 else beta.imag:.6g}')
    columns[2].metric('Guide wavelength (mm)' if beta.real>0 else 'Decay length (mm)' if beta.imag>0 else 'Guide wavelength',
        f'{2*np.pi/beta.real*1e3:.6g}' if beta.real>0 else f'{1/beta.imag*1e3:.6g}' if beta.imag>0 else 'Undefined at cutoff')
    if state=='EVANESCENT':
        st.warning('Below cutoff: no propagating power in this single decaying mode. β = iα and fields decay as exp(−αx). The reference amplitude is specified at x = 0.')
    elif state=='AT CUTOFF':
        st.info('At cutoff β = 0. The limiting field profile is shown; it carries no axial power and has no finite guide wavelength.')
    else:
        st.caption('Propagating toward +x · Phasor exp(iβx) · Time convention exp(−iωt)')
    st.caption('These m,n are rectangular-guide indices, distinct from the cylindrical angular m used in BIC Explorer. TE: Ex ≡ 0. TM: Hx ≡ 0. Other components can also vanish for particular indices or slices.')
    c1,c2,c3=st.columns([1,1.6,1.2])
    plane=c1.selectbox('Viewing plane',['y–z','x–y','x–z'])
    slice_fraction=c2.slider('Slice position (%)',0.,100.,50.,step=1.)/100
    representation=c3.selectbox('Field representation',['Magnitude |F|','Instantaneous Re(F e⁻ⁱωᵗ)','Imaginary Im(F)'])
    c1,c2,c3=st.columns([1.5,1,1])
    phase=c1.slider('Time phase ωt (degrees)',0,360,45,step=5,disabled=representation!='Instantaneous Re(F e⁻ⁱωᵗ)')
    normalized=c2.toggle('Normalize E and H separately',value=True)
    resolution=c3.select_slider('Grid points per axis',options=[61,121,181,241],value=121)
    u,v,xyz,axes=plane_grid(plane,w,h,length,slice_fraction,resolution)
    values=fields(family,m,n,float(frequency),float(w),float(h),**xyz,eps_r=eps,amplitude=amplitude)
    scales=reference_scales(family,m,n,float(frequency),float(w),float(h),eps,amplitude)
    st.markdown(f'**{plane} plane · {axes[2]} = {axes[3]*1e3:.4g} mm**')
    if plane!='y–z' and beta.real*length/(2*np.pi)>(resolution-1)/12:
        st.warning('This longitudinal span has fewer than 12 samples per guide wavelength. Shorten the shown x length or raise grid resolution to avoid aliasing.')
    st.plotly_chart(field_panels(u,v,values,axes,scales,representation,phase,normalized),width='stretch',key='six_components')
    st.caption('Analytical mode solution of the uniform, homogeneous, lossless PEC waveguide. No cylinder, metallic termination, scattering or BIC fields are included. H₀/E₀ are chosen amplitudes, not a fixed input-power normalization.')
    st.subheader('Field component inside the 3D guide')
    component=st.selectbox('3D field component',COMPONENTS,index=2,key='guide_component')
    st.caption('Rotate or zoom the guide. The colored plane uses the same slice, phase and normalization as the maps above; move it with Slice position. Color represents the selected component, not displacement of the guide.')
    st.plotly_chart(locator(w,h,length,plane,axes[3],xyz,values[component],component,
        scales[0 if component[0]=='E' else 1],representation,phase,normalized),
        width='stretch',key='plane_locator')
    if np.max(np.abs(values[component])) < 1e-12*scales[0 if component[0]=='E' else 1]:
        st.caption(f'{component} is zero on this slice for the selected mode. Choose another component or move the slice.')
    with st.expander('Reading the field plots'):
        st.write('Top row: Ex, Ey, Ez. Bottom row: Hx, Hy, Hz. Hover for values; drag to zoom.')
        st.caption('Magnitude shows the phasor envelope. Use the instantaneous view and time-phase slider to see signed lobes and longitudinal oscillation. A zero snapshot can reflect temporal phase or a nodal slice, rather than a missing mode.')
        st.caption(f'Eref = {scales[0]:.6g} V/m; Href = {scales[1]:.6g} A/m. Each is the largest component amplitude over the entire cross-section at x = 0. One shared scale per E/H row preserves component ratios and evanescent decay; no per-panel auto-normalization.')
    with st.expander('All 20 presets and cutoff frequencies'):
        rows=[]
        for fam,modes in PRESETS.items():
            for im,jn in modes:
                cut,_,status=propagation(fam,im,jn,float(frequency),float(w),float(h),eps)
                rows.append({'Mode':f'{fam}({im},{jn})','m':im,'n':jn,'Cutoff (GHz)':cut/1e9,'At current frequency':status})
        st.dataframe(sorted(rows,key=lambda row:row['Cutoff (GHz)']),hide_index=True,width='stretch')
    with st.expander('Equations and conventions'):
        st.latex(r'Y=y+w/2,\quad k_y=m\pi/w,\quad k_z=n\pi/h,\quad k_c^2=k_y^2+k_z^2,\quad \beta=\sqrt{\omega^2\mu\epsilon-k_c^2}')
        st.latex(r'\mu=\mu_0,\quad\epsilon=\epsilon_0\epsilon_r,\quad F_{physical}=\mathrm{Re}\{F(x,y,z)e^{-i\omega t}\}')
        st.markdown('**TE modes** — m,n ≥ 0, excluding (0,0).')
        st.latex(r'H_x=H_0\cos(k_yY)\cos(k_zz)e^{i\beta x},\quad E_x=0')
        st.latex(r'E_y=\frac{i\omega\mu}{k_c^2}\partial_zH_x,\quad E_z=-\frac{i\omega\mu}{k_c^2}\partial_yH_x,\quad H_y=\frac{i\beta}{k_c^2}\partial_yH_x,\quad H_z=\frac{i\beta}{k_c^2}\partial_zH_x')
        st.markdown('**TM modes** — m,n ≥ 1.')
        st.latex(r'E_x=E_0\sin(k_yY)\sin(k_zz)e^{i\beta x},\quad H_x=0')
        st.latex(r'E_y=\frac{i\beta}{k_c^2}\partial_yE_x,\quad E_z=\frac{i\beta}{k_c^2}\partial_zE_x,\quad H_y=-\frac{i\omega\epsilon}{k_c^2}\partial_zE_x,\quad H_z=\frac{i\omega\epsilon}{k_c^2}\partial_yE_x')
        st.caption('Source: Calculation_II.pdf §§5–9, pp. 17–21. At y walls Ex=Ez=0; at z walls Ex=Ey=0. These PEC conditions and both Maxwell curl equations are tested.')
    with st.expander('Download this slice'):
        buffer=io.BytesIO()
        np.savez_compressed(buffer,**values,**xyz,frequency_hz=float(frequency),w_m=float(w),h_m=float(h),
                            eps_r=eps,family=family,m=m,n=n,amplitude=amplitude,plane=plane,
                            convention='exp(i beta x - i omega t); E in V/m; H in A/m')
        st.download_button('Download complex SI fields (.npz)',buffer.getvalue(),file_name=f'{family}_{m}_{n}_slice.npz',mime='application/octet-stream')
