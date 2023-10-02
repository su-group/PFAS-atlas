import numpy as np
import pandas as pd
import tmap as tm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from tqdm import tqdm
from faerun import Faerun
import matplotlib
from rxnfp.transformer_fingerprints import (
    RXNBERTFingerprintGenerator, get_default_model_and_tokenizer, generate_fingerprints
)
import os
# os.environ["CUDA_VISIBLE_DEVICES"] = "-1"



colors = ["Fuchsia", "Turquoise", "green", 'Orange', 'red']
cmap = matplotlib.colors.ListedColormap(colors)
OECD_labels = [
    (5, "Not PFAS"),
    (4, "Other PFASs"),
    (3, "Polyfluoroalkyl acids"),
    (2, "PFAA precursors"),
    (1, "PFAAs"),

]

OECD_dict = {
    "PFAAs": 1,
    "PFAAs, cyclic": 1,
    "PFAA precursors": 2,
    "PFAA precursors, cyclic": 2,
    'Polyfluoroalkyl acids': 3,
    "Polyfluoroalkyl acids, cyclic": 3,
    "Other PFASs": 4,
    "Other PFASs, cyclic": 4,
    "Not PFAS": 5,

}


df_tmap = pd.read_csv('/home/pfas2.0/0921/input_data/OECD_4000/step3_OECD_Class_0812.csv')
list_main_smi = df_tmap['RDKIT_SMILES'].tolist()
list_main_num = df_tmap['First_Class'].apply(lambda m: OECD_dict[m]).tolist()


labels = []
for index, i in df_tmap.iterrows():
    label = (
            i['SMILES'] +
            "__<h1>" + i['Second_Class'] + "</h1>"
    )
    labels.append(label)


# compute reaction fingerprint
# model, tokenizer = get_default_model_and_tokenizer('pfas_config')
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

faerun = Faerun(coords=False, view='front', title='OECD_MAP')
faerun.add_scatter(
    "Class_OECD",
    {"x": x,
     "y": y,
     "c": [list_main_num],
     "labels": labels},
    point_scale=6.5,
    max_point_size=10,
    colormap='rainbow',
    has_legend=True,
    categorical=False,
    legend_labels=OECD_labels,
    shader='sphere'

)
faerun.add_tree("Class_OECD_tree", {"from": s, "to": t}, point_helper="Class_OECD")

faerun.plot('Bert_OECD_Map', template="smiles")
