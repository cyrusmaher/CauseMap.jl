import numpy as np
from bokeh.plotting import *
import pandas as pd
import bokeh.objects as bobj
import seaborn as sns
import pdb
from matplotlib.colors import rgb2hex

tools = "pan,wheel_zoom,box_zoom,reset,previewsave,select"
width = 500
height = 500

N = 80

x = np.linspace(0, 4*np.pi, N)
y = np.sin(x)
data = pd.read_table('vr.raw_fixed.txt', header=None)
data.columns = ['Day', 'x', 'y']

data['x_lag'] = np.nan
data.ix[1:, 'x_lag'] = data['x'].values[:-1]

data['y_lag'] = np.nan
data.ix[1:, 'y_lag'] = data['y'].values[:-1]

source = bobj.ColumnDataSource(data=dict(data))

output_file("ParaDidiExample.html", title="Para-Didi example")

paracolor, didicolor, fullcolor = [rgb2hex(xx) for xx in sns.color_palette("Set2", 3)]


lineproperties = {
                            'alpha': .3,
                            'line_width': 2.5
                            } 
scatterproperties = {
                                'tools': tools,
                                'width': width,
                                'height': height,
                                'source': source,
                                'size': 6,
                                'line_width': 0
                            }
# Create Manifold panel
maniftitle = 'Full Manifold'
hold(True)
manif = scatter('x', 'y', color=fullcolor, title=maniftitle, **scatterproperties)
manif = line(data['x'], data['y'], color=fullcolor, **lineproperties)
xaxis().axis_label = 'Paramecium abundance'
yaxis().axis_label = 'Didinium abundance'
hold(False)

# Create time series panel

##TS1
TS1title = 'Time Series (Paramecium)'
ts1 = scatter('Day', 'x', color=paracolor, title=TS1title, **scatterproperties)
hold(True)
ts1 = line(data['Day'], data['x'], color=paracolor,  **lineproperties)
xaxis().axis_label = 'Day'
yaxis().axis_label = 'Abundance'
hold(False)

## TS2
TS2title = 'Time Series (Didinium)'
ts2 = scatter('Day', 'y', color=didicolor, title=TS2title, **scatterproperties)
hold(True)
ts2 = line(data['Day'], data['y'], color=didicolor,  **lineproperties)
xaxis().axis_label = 'Day'
yaxis().axis_label = 'Abundance'
hold(False)

## x manif
manifxtitle = 'Reconstructed Manifold (Paramecium)'
manifx = scatter('x_lag', 'x', color=paracolor, title=manifxtitle, **scatterproperties)
hold(True)
manifx = line(data['x_lag'], data['x'], color=paracolor,  **lineproperties)
hold(False)

## y manif
manifytitle = 'Reconstructed Manifold (Didinium)'
manify = scatter('y_lag', 'y', color=didicolor, title=manifytitle, **scatterproperties)
hold(True)
manify = line(data['y_lag'], data['y'], color=didicolor,  **lineproperties)
hold(False)

gridplot([[ts1, ts2], [manifx, manify], [manif, manif]])

snippet= curplot().create_html_snippet()
print snippet

show()