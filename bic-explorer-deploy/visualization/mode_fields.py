"""Plotting only; input fields and scales come from physics.modes."""
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from physics.modes import COMPONENTS, displayed


def field_panels(u,v,values,axes,scales,representation,phase,normalized=True,components=COMPONENTS):
    magnitude=representation=='Magnitude |F|'
    fig=make_subplots(rows=2,cols=3,subplot_titles=list(components),
                      horizontal_spacing=.065,vertical_spacing=.22)
    for i,k in enumerate(components):
        scale=scales[0 if k[0]=='E' else 1]
        data=displayed(values[k],representation,phase)
        limit=1. if normalized else scale
        unit=('E / Eref' if k[0]=='E' else 'H / Href') if normalized else ('V/m' if k[0]=='E' else 'A/m')
        if normalized: data=data/scale
        row,col=i//3+1,i%3+1
        fig.add_trace(go.Heatmap(x=u*1e3,y=v*1e3,z=data,
            coloraxis='coloraxis' if row==1 else 'coloraxis2',
            hovertemplate=f'{axes[0]}=%{{x:.3f}} mm<br>{axes[1]}=%{{y:.3f}} mm<br>{k}=%{{z:.5g}} {unit}<extra></extra>'),row=row,col=col)
        fig.update_xaxes(title_text=f'{axes[0]} (mm)',title_standoff=8,nticks=5,automargin=True,row=row,col=col)
        fig.update_yaxes(title_text=f'{axes[1]} (mm)' if col==1 else None,title_standoff=8,nticks=5,automargin=True,row=row,col=col)
    for index,kind in enumerate(['E','H']):
        limit=1. if normalized else scales[index]
        unit=f'{kind} / {kind}ref' if normalized else ('V/m' if kind=='E' else 'A/m')
        fig.update_layout(**{'coloraxis' if index==0 else 'coloraxis2':dict(
            cmin=0 if magnitude else -limit,cmax=limit,
            colorscale='Viridis' if magnitude else 'RdBu_r',
            colorbar=dict(title=dict(text=unit,side='top'),len=.38,thickness=14,
                          x=1.025,y=.805 if index==0 else .195,tickformat='.2g'))})
    fig.update_annotations(font=dict(size=18,color='#263648'))
    fig.update_layout(template='plotly_white',height=640,font=dict(size=13,color='#334155'),
                      margin=dict(l=55,r=95,t=40,b=45))
    return fig


def locator(w,h,length,plane,position,xyz,field,component,scale,representation,phase,normalized=True):
    fig=go.Figure()
    for axis in range(3):
        other=[i for i in range(3) if i!=axis]
        bounds=[(0,length),(-w/2,w/2),(0,h)]
        for a in bounds[other[0]]:
            for b in bounds[other[1]]:
                pts=np.zeros((2,3)); pts[:,axis]=bounds[axis];pts[:,other[0]]=a;pts[:,other[1]]=b
                fig.add_trace(go.Scatter3d(x=pts[:,0]*1e3,y=pts[:,1]*1e3,z=pts[:,2]*1e3,
                    mode='lines',line=dict(color='#607080',width=3),showlegend=False,hoverinfo='skip'))
    data=displayed(field,representation,phase)
    if normalized:
        data=data/scale
    limit=1. if normalized else scale
    unit=f'{component[0]} / {component[0]}ref' if normalized else ('V/m' if component[0]=='E' else 'A/m')
    magnitude=representation=='Magnitude |F|'
    fig.add_trace(go.Surface(x=xyz['x']*1e3,y=xyz['y']*1e3,z=xyz['z']*1e3,
        surfacecolor=data,cmin=0 if magnitude else -limit,cmax=limit,
        colorscale='Viridis' if magnitude else 'RdBu_r',opacity=1,
        colorbar=dict(title=dict(text=unit,side='top'),thickness=16,len=.7),
        lighting=dict(ambient=1,diffuse=0,specular=0,fresnel=0),
        hovertemplate=f'x=%{{x:.3f}} mm<br>y=%{{y:.3f}} mm<br>z=%{{z:.3f}} mm<br>{component}=%{{surfacecolor:.5g}} {unit}<extra></extra>'))
    fig.update_layout(height=520,margin=dict(l=20,r=65,t=25,b=25),
        title=dict(text=f'{component} · {plane} slice',font=dict(size=17)),
        uirevision='guide-camera',
        scene=dict(xaxis_title='x (mm) →',yaxis_title='y (mm)',zaxis_title='z (mm)',
                   aspectmode='data',camera=dict(eye=dict(x=1.65,y=2.,z=1.6))))
    return fig
