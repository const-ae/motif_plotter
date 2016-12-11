# Motif plotter

This library provides the tools to create plots of DNA motifs. It 
builds upon the matplotlib framework and provides flexible parts to
construct individual plots with letters that represent data.

# Usage

Highlevel API example:

Example of a single sequence:
```python
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import motif_plotter

# This is a function in cross_aligner/parse_helper.py
def parse_importance_df(df, col_names):
    # Iterate over every entry
    parsed_cols = []
    for name in col_names:
        col = df[name].as_matrix()
        parsed_col = np.apply_along_axis(lambda e: np.array([float(x) for x in e[0][1:-1].split(",")]), 1, col.reshape(len(col),1))
        parsed_cols.append(parsed_col)
    return np.stack(parsed_cols, 2)

importance_file = "examples/sequence_only-importances.tsv"
importance_df = pd.read_csv(importance_file, sep="\t")
values = parse_importance_df(importance_df, ["A", "T", "C", "G"])
scores = parse_importance_df(importance_df, ["A_Scores", "T_Scores", "C_Scores", "G_Scores"])

fig=plt.figure()
ax=fig.add_subplot(111)

motif_plotter.make_single_sequence_spectrum(ax, values[0], scores[0])

plt.show()
```


Example of a ConsensusMotif:
```python
from Bio import motifs

with open("examples/transfac_motif.txt") as handle:
    m = motifs.parse(handle, "transfac")


fig=plt.figure()
ax=fig.add_subplot(111)

cbp = ConsensusMotifPlotter.from_bio_motif(m[0], scale_info_content=True)
cbp.plot(ax)

plt.show()
```


Low level API:

```python
import matplotlib.pylab as plt
from motif_plotter import *

fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.set_ylim(-0.5,1.5)
ax1.set_xlim(-0.5,1.5)

sentence_shape = make_text_elements('Hello', x=-0.5,y=0.25,width=0.5,height=0.5)
a_shape = make_text_elements('A', x=0,   y=0.2,  width=0.5)
b_shape = make_text_elements('B', x=0.5, y=-0.2, width=0.5)
c_shape = make_text_elements('C', x=1,   y=0.6,  width=0.5, height=0.9)
i_shape = make_text_elements('I', x=1, y=-0.5, width=0.5, height=1.1)
ax1.add_patch(sentence_shape)
ax1.add_patch(a_shape)
ax1.add_patch(b_shape)
ax1.add_patch(c_shape)
ax1.add_patch(i_shape)
plt.show()
```