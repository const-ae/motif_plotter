'''
This code is mostly taken from this discussion: https://github.com/biopython/biopython/issues/850
But as the plotting didn't seem to be scalable, I will probably not use it
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtrans
from matplotlib import transforms
import seaborn
from motif_plotter import *

from Bio import motifs
with open("examples/transfac_motif.txt") as handle:
    m = motifs.parse(handle, "transfac")


def approximate_error(motif):
    """Calculate approximate error"""
    pwm = motif.pwm
    bases = list(pwm.keys())
    n = sum(motif.counts[bases[0]])
    approx_error = (len(bases)-1)/(2 * np.log(2) * n)
    return approx_error


def exact_error(motif):
    """Calculate exact error, using multinomial(na,nc,ng,nt)"""
    ## Super Slow. O(n^3)
    pwm = motif.pwm
    bases = pwm.keys()
    na = sum(motif.counts['A'])
    n = na
    nc = 0
    ng = 0
    nt = 0
    done = False
    exact_error = 0
    while not done:
        print (na,nc,ng,nt)
        exact_error += sum([-p*np.log2(p) for p in [na/n, nc/n, ng/n, nt/n]])
        if nt<=0:
            ## iterate inner loop
            if ng > 0:
                ## g => t
                ng = ng - 1
                nt = nt + 1
            elif nc > 0:
                ## c -> g
                nc = nc - 1;
                ng = ng + 1;
            else:
                ## a->c
                na = na - 1
                nc = nc + 1
        else:
            if ng > 0:
                ## g => t
                ng = ng - 1
                nt = nt + 1
            elif nc>0:
                ## c => g; all t -> g
                nc = nc - 1
                ng = nt + 1
                nt = 0
            elif na>0:
                ## a => c; all g,t -> c
                nc = nt + 1
                na = na - 1
                nt = 0
            else:
                done = True
    return exact_error


def calc_info_matrix(motif, correction_type='approx'):
    """Calculate information matrix with small sample correction"""
    pwm = motif.pwm
    bases = pwm.keys()
    if correction_type=='approx':
        error = approximate_error(motif)
    else:
        error = exact_error(motif)
    info_matrix = [2-error+sum([pwm[b][l]*np.nan_to_num(np.log2(pwm[b][l])) for b in bases]) for l in range(0, len(motif))]
    return info_matrix

def calc_relative_information(motif, correction_type='approx'):
    """Calculate relative information matrix"""
    pwm = motif.pwm
    bases = pwm.keys()
    if correction_type=='approx':
        info_matrix = calc_info_matrix(motif)
    else:
        info_matrix = calc_info_matrix(motif, 'exact')
    relative_info = {base: [prob*info for prob,info in zip(pwm[base], info_matrix)]  for base in bases}
    return relative_info


motif = m[1]
rel_info = calc_relative_information(motif)
bases = ['A', 'T', 'G', 'C']
colors_scheme = {'G': '#ffb300', 'T': '#008000', 'C': '#0000cc', 'A': '#cc0000'}
#ax.set_xticks(range(0.1,0.1*len(motif),0.1))
#ax.set_xticklabels(range(1,len(motif)+1))

fig = plt.figure()
fig.set_size_inches(len(motif)*1.3,3)
ax = fig.add_subplot(111)

xshift = 0
trans_offset = transforms.offset_copy(ax.transAxes,
                                  fig=fig,
                                  x=0,
                                  y=0,
                                  units='dots')

for i in range(0, len(motif)):
    scores = [(b,rel_info[b][i]) for b in bases]
    scores.sort(key=lambda t: t[1])
    yshift = 0
    for base, score in scores:
        txt = ax.text(0,
                      0,
                      base,
                      transform=trans_offset,
                      fontsize=82,
                      color=colors_scheme[base],
                      family='monospace'
                      )
        txt.set_path_effects([Scale(1, score)])
        fig.canvas.draw()
        pos_bounding_box = txt.get_window_extent(txt._renderer)
        yshift = pos_bounding_box.height * score
        trans_offset = transforms.offset_copy(txt._transform, fig=fig, y=yshift, units='dots')
    xshift+=pos_bounding_box.width
    trans_offset = mtrans.offset_copy(ax.transAxes, fig=fig, x=xshift, units='dots')

ax.set_yticks(range(0,3))
#ax.set_xticks(np.arange(window_ext.width/100,len(motif)*window_ext.width/100+0.5,window_ext.width/100))

#ax.set_xticklabels(range(1, len(motif)+1))
ax.set_yticklabels(np.arange(0,3,1))

plt.show()
plt.clf()