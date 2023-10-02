import numpy as np
import pandas as pd
import tmap as tm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from tqdm import tqdm
import math
from faerun import Faerun
import matplotlib
from rxnfp.transformer_fingerprints import (
    RXNBERTFingerprintGenerator, get_default_model_and_tokenizer, generate_fingerprints
)



colors = ['DarkOrchid', 'BlueViolet', 'DarkViolet', 'RoyalBlue', 'DeepSkyBlue',
          'DarkTurquoise', 'MediumTurquoise',
          'MediumSpringGreen', 'SpringGreen', 'LimeGreen', 'Chartreuse',
          'GreenYellow', 'Gold', 'Orange', 'DarkOrange', 'OrangeRed', 'Red', 'gray']
cmap = matplotlib.colors.ListedColormap(colors)
VDR_labels = [
    (1, 'A: -14~-13'),
    (2, 'A: -13~-12'),
    (3, 'A: -12~--11'),
    (4, 'A: -11~-10'),
    (5, 'A: -10~-9'),
    (6, 'A: -9~-8'),
    (7, 'A: -8~-7'),
    (8, 'A: -7~-6'),
    (9, 'A: -6~-5'),
    (10, 'A: -5~-4'),
    (11, 'A: -4~-3'),
    (12, 'A: -3~-2'),
    (13, 'A: -2~-1'),
    (14, 'A: -1~0'),
    (15, 'A: 0~1'),
    (16, 'A: 1~2'),
    (17, 'A: 2~3'),
    (18, 'no docking_score data'),
]


#df_tmap = pd.read_csv('input_data/OECD_4000/step3_PFASOECD_class_0801.csv')
df_tmap = pd.read_csv('/home/pfas2.0/0921/apply_PFAS/VDR/VDR_OECD_0812_review.csv')
list_main_smi = df_tmap['RDKIT_SMILES'].tolist()
list_main_num = df_tmap['range of docking_score'].tolist()


labels = []
for index, i in df_tmap.iterrows():
    if not math.isnan(i['docking_score']):
        label = (
                i['SMILES'] +
                "__<h1>" + i['Second_Class'] + '<br />' +
                'docking_score:' + '<br />' +
                str(round(i['docking_score'], 3)) + "</h1>"
        )
    else:
        label = (
                i['SMILES'] +
                "__<h1>" + i['Second_Class'] + '<br />' +
                'no data' + '<br />' + "</h1>"
        )
    labels.append(label)


# compute reaction fingerprint
model, tokenizer = get_default_model_and_tokenizer('best_model')
rxnfp_generator = RXNBERTFingerprintGenerator(model, tokenizer)
pfas_fp = generate_fingerprints(list_main_smi, rxnfp_generator, batch_size=16)
np.savez_compressed('pfas', fps=pfas_fp)
pfas_fp = np.load('pfas.npz')['fps']


dims = 256
enc = tm.Minhash(dims)

lf = tm.LSHForest(dims, 128)
cfg = tm.LayoutConfiguration()
cfg.k = 50
cfg.kc = 50
cfg.sl_scaling_min = 1.0
cfg.sl_scaling_max = 1.0
cfg.sl_repeats = 1
cfg.sl_extra_scaling_steps = 2
cfg.placer = tm.Placer.Barycenter
cfg.merger = tm.Merger.LocalBiconnected
cfg.merger_factor = 2.0
cfg.merger_adjustment = 0
cfg.fme_iterations = 1000
cfg.sl_scaling_type = tm.ScalingType.RelativeToDesiredLength
cfg.node_size = 1 / 37
cfg.mmm_repeats = 1


# fingerprints = [tm.VectorUint(enc.encode(s)) for s in list_main_smi]
fingerprints = [enc.from_weight_array(fp.tolist(), method="I2CWS") for fp in tqdm(pfas_fp)]
lf.batch_add(fingerprints)
lf.index()
x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=cfg)

faerun = Faerun(coords=False, view='front', title='VDR_OECD_MAP')
faerun.add_scatter(
    "docking_score",
    {"x": x,
     "y": y,
     "c": [list_main_num],
     "labels": labels},
    point_scale=6.5,
    max_point_size=10,
    colormap=cmap,
    has_legend=True,
    categorical=False,
    legend_labels=VDR_labels,
    shader='sphere'

)
faerun.add_tree("docking_score_tree", {"from": s, "to": t}, point_helper="docking_score")

faerun.plot('Bert_VDR_OECD_Map', template="smiles")
