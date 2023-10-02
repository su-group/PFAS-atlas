import numpy as np
import pandas as pd
import tmap as tm
from tqdm import tqdm
from faerun import Faerun
from rxnfp.transformer_fingerprints import (
    RXNBERTFingerprintGenerator, get_default_model_and_tokenizer, generate_fingerprints
)


df_tmap = pd.read_csv('/home/ubuntu/pythonProject/PFAS-2.0/input_data/OECD_4000/step3_OECD_Class_0812.csv')
others_labels = [
    (12, "others"),
    (11, 'PFPEs'),
    (10, "PFAS derivatives"),
    (9, '(Hg,Sn,Ge,Sb,Se,B) PFASs'),
    (8, "Si PFASs"),
    (7, "Amide derivatives"),
    (6, "Polyfluoroalkyl ethers and derivatives"),
    (5, "Aromatic PFASs"),
    (4, "Perfluoroalkyl-tert-amines"),
    (3, "Perfluoroalkylethers"),
    (2, "Polyfluoroalkanes"),
    (1, "Perfluoroalkanes")
]
others_dict = {
    "Perfluoroalkanes": 1,
    'Perfluoroalkanes, cyclic': 1,
    "Polyfluoroalkanes": 2,
    'Polyfluoroalkanes, cyclic': 2,
    "Perfluoroalkylethers": 3,
    'Perfluoroalkylethers, cyclic': 3,
    "Perfluoroalkyl-tert-amines": 4,
    'Perfluoroalkyl-tert-amines, cyclic': 4,
    "Aromatic PFASs": 5,
    "Polyfluoroalkyl ethers and derivatives": 6,
    'Polyfluoroalkyl ethers and derivatives, cyclic': 6,
    "Amide derivatives": 7,
    'Amide derivatives, cyclic': 7,
    "Si PFASs": 8,
    'Si PFASs, cyclic': 8,
    "(Hg,Sn,Ge,Sb,Se,B) PFASs": 9,
    "PFAS derivatives": 10,
    'PFPEs': 11,
    'PFPEs, cyclic': 11,
    "others": 12,
    'others, cyclic': 12

}
df_other_subclass = df_tmap[(df_tmap['First_Class'] == 'Other PFASs') | (df_tmap['First_Class'] == 'Other PFASs, cyclic')]
list_main_smi = df_other_subclass['RDKIT_SMILES'].tolist()
list_main_num = df_other_subclass['Second_Class'].apply(lambda m: others_dict[m]).tolist()



labels = []
for index, i in df_other_subclass.iterrows():
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

faerun = Faerun(coords=False, view='front', title='other_PFASs_MAP')
faerun.add_scatter(
    "Other_PFASs",
    {"x": x,
     "y": y,
     "c": [list_main_num],
     "labels": labels},
    point_scale=9,
    max_point_size=10,
    colormap='rainbow',
    has_legend=True,
    categorical=False,
    legend_labels=others_labels,
    shader='sphere'

)
faerun.add_tree("Other_PFASs_tree", {"from": s, "to": t}, point_helper="Other_PFASs")

faerun.plot('Bert_other_PFASs_Map', template="smiles")
