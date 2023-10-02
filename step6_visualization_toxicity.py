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


colors = ['LightSkyBlue', 'DeepSkyBlue', 'MediumSpringGreen', 'Lime', 'GreenYellow', 'yellow', 'Gold', 'Orange', 'OrangeRed', 'red', 'gray']
cmap = matplotlib.colors.ListedColormap(colors)

hit_labels = [
    (0, '0'),
    (1, 'A: 0 ~ 0.05'),
    (2, 'A: 0.05 ~ 0.1'),
    (3, 'A: 0.1 ~ 0.15'),
    (4, 'A: 0.15 ~ 0.2'),
    (5, 'A: 0.2 ~ 0.25'),
    (6, 'A: 0.25 ~ 0.3'),
    (7, 'A: 0.3 ~ 0.35'),
    (8, 'A: 0.35 ~ 0.4'),
    (9, 'A: 0.4 ~ 0.45'),
    (10, 'no activity data')
]

#df_tmap = pd.read_csv('input_data/OECD_4000/step3_PFASOECD_class_0801.csv')
df_tmap = pd.read_csv('/home/pfas2.0/0921/apply_PFAS/150_toxicity/150_OECD_0812_review.csv')
list_main_smi = df_tmap['RDKIT_SMILES'].tolist()
list_main_num = df_tmap['range of hit_ratio'].tolist()

labels = []
for index, i in df_tmap.iterrows():
    if not math.isnan(i['hit_ratio']):
        label = (
                i['SMILES'] +
                "__<h1>" + i['Second_Class'] + '<br />' +
                'hit_ratio:' + '<br />' +
                str(round(i['hit_ratio'], 2)) + "</h1>"
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

faerun = Faerun(coords=False, view='front', title='150hit_OECD_MAP')
faerun.add_scatter(
    "hit_ratio",
    {"x": x,
     "y": y,
     "c": [list_main_num],
     "labels": labels},
    point_scale=6.5,
    max_point_size=10,
    colormap=cmap,
    has_legend=True,
    categorical=False,
    legend_labels=hit_labels,
    shader='sphere'

)
faerun.add_tree("hit_ratio_tree", {"from": s, "to": t}, point_helper="hit_ratio")

faerun.plot('Bert_150hit_OECD_Map', template="smiles")
