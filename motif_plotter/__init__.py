from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
import matplotlib.patches as patches
from matplotlib.transforms import Affine2D
from motif_plotter.information_helper import *


def make_text_elements(text, x=0.0, y=0.0, width=1.0, height=1.0, color='blue', edgecolor="black",
                       font = FontProperties(family='monospace')):
    tp = TextPath((0.0, 0.0), text, size=1, prop=font)
    bbox = tp.get_extents()
    bwidth = bbox.x1 - bbox.x0
    bheight = bbox.y1 - bbox.y0
    trafo = Affine2D()
    trafo.translate(-bbox.x0, -bbox.y0)
    trafo.scale(1 / bwidth * width, 1 / bheight * height)
    trafo.translate(x,y)
    tp = tp.transformed(trafo)
    return patches.PathPatch(tp, facecolor=color, edgecolor=edgecolor)


def make_bar_plot(axes, texts, heights, width=0.8, colors=None, edgecolor="black"):
    """
    Makes a bar plot but each bar is not just a rectangle but an element from the texts list
    :param axes: the axes that is modified
    :param texts: a list of strings, where each element is plotted as a "bar"
    :param heights: a list of the height of each texts element
    :param width: the width of the bar. Default: 0.8
    :param colors: A list of colors, a list with a single entry or None. Default: None, which is plotted as blue
    :return: None
    """
    texts = list(texts)
    heights = list(heights)
    n_elem = len(texts)
    if n_elem != len(heights):
        raise ValueError("Texts and heights must be of the same length")
    if colors is None:
        colors = ['blue'] * n_elem
    elif len(colors) == 1:
        colors *= n_elem

    axes.set_ylim(min(0,min(heights)), max(0,max(heights)))
    axes.set_xlim(0, n_elem)
    for idx, (text, height, color) in enumerate(zip(texts, heights, colors)):
        text_shape = make_text_elements(text, x=idx+(1-width)/2, y=0, width=width, height=height,
                                        color=color, edgecolor=edgecolor)
        axes.add_patch(text_shape)


def make_stacked_bar_plot(axes, texts, heights, width=0.8, colors=None, edgecolor="black"):
    """
    Makes a stackedbar plot but each bar is not just a rectangle but an element from the texts list
    :param axes: the axes that is modified
    :param texts: a list of list of strings, where each element is plotted as a "bar"
    :param heights: a list of lists of the height of each texts element
    :param width: the width of the bar. Default: 0.8
    :param colors:
    :return: None
    """
    if colors is None:
        colors = [['blue'] * len(text) for text in texts]
    elif len(colors) == 1:
        colors = [colors * len(text) for text in texts]

    if len(texts) != len(heights):
        raise ValueError("Texts and heights must be of the same length")
    for idx, (text, height, color) in enumerate(zip(texts, heights, colors)):
        y_stack = 0
        for jdx, (t, h, c) in enumerate(zip(text, height, color)):
            text_shape = make_text_elements(t, x=idx+(1-width)/2, y=y_stack, width=width, height=h,
                                            color=c, edgecolor=edgecolor)
            axes.add_patch(text_shape)
            y_stack += h

    axes.autoscale()
    axes.set_xlim(0, len(texts))


class ConsensusMotifPlotter:

    def __init__(self, elements, weights, colors=None):
        self.n_elem = len(elements)
        self.colors = colors
        self.elements = elements
        self.weights = weights

    @classmethod
    def from_importance_scoring(cls, value):
        return cls([['A', 'C', 'T', 'G']] * len(value.Sequence), value.Scores,
                   [['#008000', '#0000cc', '#cc0000','#ffb300']]* len(value.Sequence))

    @classmethod
    def from_bio_motif(cls, motif, scale_info_content=True):
        n_elem = len(motif)
        colors_scheme = {'G': '#ffb300', 'A': '#008000', 'C': '#0000cc', 'T': '#cc0000'}
        bases = ['A', 'T', 'G', 'C']
        if scale_info_content:
            rel_info = calc_relative_information(motif)
        else:
            rel_info = motif.counts

        basess = []
        scoress = []
        colorss = []
        for i in range(0, n_elem):
            scores = [(b, rel_info[b][i], colors_scheme[b]) for b in bases]
            scores.sort(key=lambda t: t[1])
            basess += [[x[0] for x in scores]]
            scoress += [[x[1] for x in scores]]
            colorss += [[x[2] for x in scores]]
        return cls(basess, scoress, colorss)


    def plot(self, axes):
        """
        Add the motif to an axes
        :return: modifies the axes object with all the necessary characters
        """
        make_stacked_bar_plot(axes, self.elements, self.weights, width=1, colors=self.colors, edgecolor="none")