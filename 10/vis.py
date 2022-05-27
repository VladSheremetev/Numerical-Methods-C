import re
import sys

import plotly.graph_objects as go
import matplotlib.pyplot as plt


def main():
    filename = sys.argv[1]
    l = []
    with open(filename) as f:
        for line in f:
            l.append(list(map(float, line.split())))
    fig = go.Figure(data=[go.Surface(z=l[2:], y=l[0], x=l[1])])
    name = re.sub('.txt', '', filename)
    fig.update_layout(title=f'{name}', autosize=True,
                      scene = dict(xaxis_title='t',
                                   yaxis_title='x',
                                  ),
                      )
    fig.update_xaxes(title='t')
    fig.show()


main()
