from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
import matplotlib.patches as patches
from matplotlib.transforms import Affine2D


def make_text_elements(text, x=0, y=0, width=1, height=1, font = FontProperties(family='monospace')):
    tp1 = TextPath((0.0, 0.0), text, size=1, prop=font)
    bbox = tp1.get_extents()
    bwidth = bbox.x1 - bbox.x0
    bheight = bbox.y1 - bbox.y0
    trafo = Affine2D()
    trafo.translate(-bbox.x0, -bbox.y0)
    trafo.scale(1 / bwidth * width, 1 / bheight * height)
    trafo.translate(x,y)
    tp1 = tp1.transformed(trafo)
    return patches.PathPatch(tp1)