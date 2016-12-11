import matplotlib.pyplot as plt
from motif_plotter import *
import numpy as np

fig=plt.figure()
ax1=fig.add_subplot(111)
ax1.set_ylim(-0.5,1.5)
ax1.set_xlim(-0.5,1.5)

sentence_shape = make_text_elements('Hello', x=-0.5,y=0.25,width=0.5,height=0.5)
a_shape = make_text_elements('A', x=0,   y=0.2,  width=0.5)
b_shape = make_text_elements('B', x=0.5, y=-0.2, width=0.5)
c_shape = make_text_elements('C', x=1,   y=0.6,  width=0.5)
ax1.add_patch(sentence_shape)
ax1.add_patch(a_shape)
ax1.add_patch(b_shape)
ax1.add_patch(c_shape)
plt.show()


fig=plt.figure()
ax=fig.add_subplot(111)
make_stacked_bar_plot(ax, [list("abcd"), list("defg")], np.array([[1,-2,3,4], [4,0,-2,0]]))
plt.show()





from Bio import motifs

with open("examples/transfac_motif-negative.txt") as handle:
    m = motifs.parse(handle, "transfac")


fig=plt.figure()
ax=fig.add_subplot(111)

cbp = ConsensusMotifPlotter.from_bio_motif(m[0], scale_info_content=False)
cbp.plot(ax)

plt.show()




import pandas as pd
from copy import deepcopy
import motif_plotter

location = []
sequence = []
scores = []
with open("examples/examples_scores.tsv", "rb") as file:
    for line in file:
        elem = line.decode("UTF-8").strip().split("\t")
        location.append(elem[0])
        sequence.append(elem[1])
        scores.append(np.array(elem[2:]).reshape((int((len(elem)-2)/4),4)).astype(np.float))

data_frame = pd.DataFrame({"Location": location, "Sequence": sequence, "Scores": scores})

fig=plt.figure()
ax=fig.add_subplot(111)

value = deepcopy(data_frame.iloc[1])
value.Scores = value.Scores[50:100]

value.Scores[1, :] = [1, -2, -0.4, 0.5]

value.Sequence = value.Sequence[50:100]
cbp = motif_plotter.ConsensusMotifPlotter.from_importance_scoring(value)
cbp.plot(ax)

plt.show()



import pandas as pd
import numpy as np
import motif_plotter

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

fig=plt.figure(figsize=(20, 10))
for i in range(0, 10):
    ax=fig.add_subplot(10, 1, i+1)
    motif_plotter.make_single_sequence_spectrum(ax, values[100+i], scores[100+i])

plt.show()