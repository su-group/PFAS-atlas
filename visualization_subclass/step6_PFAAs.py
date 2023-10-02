import rdkit
from rdkit import Chem
import pandas as pd
import numpy as np
import pandas as pd
import tmap as tm
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from tqdm import tqdm
from mhfp.encoder import MHFPEncoder
from faerun import Faerun
import matplotlib
from matplotlib.colors import ListedColormap
from rxnfp.transformer_fingerprints import (
    RXNBERTFingerprintGenerator, get_default_model_and_tokenizer, generate_fingerprints
)
import matplotlib.pyplot as plt

df_tmap = pd.read_csv('/home/ubuntu/pythonProject/PFAS-2.0/input_data/OECD_4000/step3_OECD_Class_0812.csv')
PFAAs_labels = [
    (11, "PFCA-ester derivatives"),
    (10, 'PFSA derivatives'),
    (9, "PFCA-anhydrides"),
    (8, "PFdiCAs"),
    (7, "PFSIAs"),
    (6, "PFPIAs"),
    (5, "PFECAs"),
    (4, "PFESAs"),
    (3, "PFPAs"),
    (2, "PFSAs"),
    (1, "PFCAs"),
]
PFAAs_dict = {
    "PFCAs": 1,
    "PFSAs": 2,
    'PFSAs, cyclic': 2,
    "PFSIAs": 7,
    "PFECAs": 5,
    'PFECAs, cyclic': 5,
    "PFESAs": 4,
    "PFPAs": 3,
    "PFPIAs": 6,
    "PFdiCAs": 8,
    "PFCA-anhydrides": 9,
    'PFSA derivatives': 10,
    'PFSA derivatives, cyclic': 10,
    "PFCA-ester derivatives": 11,
    'PFCA-ester derivatives, cyclic': 11,

}
df_PFAAs_subclass = df_tmap[(df_tmap['First_Class'] == 'PFAAs') | (df_tmap['First_Class'] == 'PFAAs, cyclic')]
list_main_smi = df_PFAAs_subclass['RDKIT_SMILES'].tolist()
list_main_num = df_PFAAs_subclass['Second_Class'].apply(lambda m: PFAAs_dict[m]).tolist()


labels = []
for index, i in df_PFAAs_subclass.iterrows():
    label = (
            i['SMILES'] +
            "__<h1>" + i['Second_Class'] + "</h1>"
    )
    labels.append(label)


#c_list = df['EGP'].tolist()
# compute reaction fingerprint
model, tokenizer = get_default_model_and_tokenizer('pfas_config')
rxnfp_generator = RXNBERTFingerprintGenerator(model, tokenizer)
pfas_fp = generate_fingerprints(list_main_smi, rxnfp_generator, batch_size=16)
np.savez_compressed('pfas', fps=pfas_fp)
pfas_fp = np.load('pfas.npz')['fps']



dims = 256   # 256
enc = tm.Minhash(dims)

lf = tm.LSHForest(dims, 128)   # 128
cfg = tm.LayoutConfiguration()
cfg.k = 30
cfg.kc = 30  # 50
cfg.sl_scaling_min = 1.0
cfg.sl_scaling_max = 1.0
cfg.sl_repeats = 1
cfg.sl_extra_scaling_steps = 2
cfg.placer = tm.Placer.Barycenter
cfg.merger = tm.Merger.LocalBiconnected
cfg.merger_factor = 2.0
cfg.merger_adjustment = 0
cfg.fme_iterations = 1000    # 1000
cfg.sl_scaling_type = tm.ScalingType.RelativeToDesiredLength
cfg.node_size = 1 / 50       # 37
cfg.mmm_repeats = 1



# fingerprints = [tm.VectorUint(enc.encode(s)) for s in list_main_smi]
fingerprints = [enc.from_weight_array(fp.tolist(), method="I2CWS") for fp in tqdm(pfas_fp)]
lf.batch_add(fingerprints)
lf.index()
x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=cfg)

faerun = Faerun(coords=False, view='front', title='PFAAs_MAP')
faerun.add_scatter(
    "PFAAs",
    {"x": x,
     "y": y,
     "c": [list_main_num],
     "labels": labels},
    point_scale=10,
    max_point_size=10,
    colormap='rainbow',
    has_legend=True,
    categorical=False,
    legend_labels=PFAAs_labels,
    #series_title="dGsolv",
    shader='sphere'

)
faerun.add_tree("PFAAs_tree", {"from": s, "to": t}, point_helper="PFAAs")
faerun.plot('Bert_PFAAs_Map', template="smiles")
