from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
import matplotlib.patches as patches
from matplotlib.transforms import Affine2D


def make_text_elements(text, x=0.0, y=0.0, width=1.0, height=1.0, font = FontProperties(family='monospace')):
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


def make_bar_plot(axes, texts, heights, width=0.8):
    '''
    Makes a bar plot but each bar is not just a rectangle but an element from the texts list
    :param axes: the axes that is modified
    :param texts: a list of strings, where each element is plotted as a "bar"
    :param heights: a list of the height of each texts element
    :param width: the width of the bar. Default: 0.8
    :return: None
    '''
    texts = list(texts)
    heights = list(heights)
    n_elem = len(texts)
    if n_elem != len(heights):
        raise ValueError("Texts and heights must be of the same length")

    axes.set_ylim(min(0,min(heights)), max(0,max(heights)))
    axes.set_xlim(0, n_elem)
    for idx, (text, height) in enumerate(zip(texts, heights)):
        text_shape = make_text_elements(text, x=idx+(1-width)/2, y=0, width=width, height=height)
        axes.add_patch(text_shape)


def make_stacked_bar_plot(axes, texts, heights, width=0.8):
    '''
    Makes a stackedbar plot but each bar is not just a rectangle but an element from the texts list
    :param axes: the axes that is modified
    :param texts: a list of list of strings, where each element is plotted as a "bar"
    :param heights: a list of lists of the height of each texts element
    :param width: the width of the bar. Default: 0.8
    :return: None
    '''
    n_elem = len(texts)
    if n_elem != len(heights):
        raise ValueError("Texts and heights must be of the same length")
    for idx, (text, height) in enumerate(zip(texts, heights)):
        y_stack = 0
        for jdx, (t, h) in enumerate(zip(text, height)):
            text_shape = make_text_elements(t, x=idx+(1-width)/2, y=y_stack, width=width, height=h)
            axes.add_patch(text_shape)
            y_stack += h

    axes.autoscale()