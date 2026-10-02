import plotly.graph_objects as go
from physics.units import m_to_mm


def schematic(w,h,a,b,L):
    w,h,a,b,L = [float(m_to_mm(x)) for x in [w,h,a,b,L]]
    right=max(w,2*b); left=-L
    fig=go.Figure()
    for y in [-w/2,w/2]:
        fig.add_shape(type='line',x0=left,x1=right,y0=y,y1=y,line=dict(color='#444',width=3))
    fig.add_shape(type='line',x0=left,x1=left,y0=-w/2,y1=w/2,line=dict(color='#222',width=6))
    for radius,color in [(b,'#90bed1'),(a,'#757575')]:
        fig.add_shape(type='circle',x0=-radius,x1=radius,y0=-radius,y1=radius,
                      fillcolor=color,line_color='#334')
    fig.add_annotation(x=left,y=w*.65,text='PEC wall',showarrow=False)
    fig.add_annotation(x=right,y=0,text='Open TE10 port',showarrow=False,xanchor='right')
    fig.add_annotation(x=-L/2,y=-w*.35,text=f'L = {L:.3f} mm',showarrow=False)
    fig.add_annotation(x=right*.65,y=w*.35,ax=right*.25,ay=w*.35,axref='x',ayref='y',text='+x',showarrow=True)
    fig.update_layout(template='plotly_white',height=280,margin=dict(l=25,r=20,t=20,b=30),
        xaxis=dict(title='x (mm)',range=[left-w*.15,right+w*.1]),
        yaxis=dict(title='y (mm)',scaleanchor='x',range=[-w*.7,w*.8]))
    return fig
