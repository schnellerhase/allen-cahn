import pyvista as pv
import argparse

pv.start_xvfb(wait=0.1)
parser = argparse.ArgumentParser()
parser.add_argument('file', type=str)
file = parser.parse_args().file
    
pv.global_theme.cmap = 'jet'
pv.global_theme.transparent_background = True

mesh = pv.read(file)
pl = pv.Plotter(window_size=[720, 720], off_screen=True)
pl.add_mesh(mesh, scalar_bar_args=dict(width=.9, position_x=(1-.9)/2))
pl.view_isometric()
pl.set_position([.5, .5, 1])
pl.set_viewup([0, 1, 0])
pl.show(screenshot='plot.png')
