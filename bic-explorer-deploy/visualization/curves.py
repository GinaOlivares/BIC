import plotly.graph_objects as go


def curves(x, series, xlabel, ylabel):
    fig = go.Figure()
    for name, y in series.items():
        fig.add_trace(go.Scatter(x=x, y=y, name=name, mode='lines'))
    fig.update_layout(template='plotly_white', height=360, margin=dict(l=30,r=20,t=30,b=40),
                      xaxis_title=xlabel, yaxis_title=ylabel, legend=dict(orientation='h'))
    return fig


def markers(fig, positions, label='BIC candidate', color='#577590'):
    for x in positions:
        fig.add_vline(x=float(x), line_dash='dot', line_color=color)
    if positions:
        fig.add_annotation(x=float(positions[0]), y=1, yref='paper', text=label, showarrow=False)
    return fig
