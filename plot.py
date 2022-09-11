import pyvista as pv
import argparse
import sys
# must change font!

pv.start_xvfb(wait=0.1)
parser = argparse.ArgumentParser()
parser.add_argument("file", type=str)
file = parser.parse_args().file

pv.global_theme.cmap = "jet"
pv.global_theme.transparent_background = True
pv.global_theme.font.size = 10

mesh = pv.read(file)
pl = pv.Plotter(window_size=[1920, 1920], off_screen=True)
mesh = pl.add_mesh(
    mesh, scalar_bar_args=dict(width=0.9, position_x=(1 - 0.9) / 2), clim=[-1, 1]
)

pl.view_isometric()
pl.set_position([0.5, 0.5, 1])
pl.set_viewup([0, 1, 0])
if sys.version_info[0] > 2:
    pl.camera.zoom("tight")
pl.save_graphic(f"{file}.svg")
pl.show(screenshot=f"{file}.png")
