import numpy as np
import pandas as pd
import tmap as tm
from tqdm import tqdm
from faerun import Faerun
from rxnfp.transformer_fingerprints import (
    RXNBERTFingerprintGenerator, get_default_model_and_tokenizer, generate_fingerprints
)


df_tmap = pd.read_csv('/home/ubuntu/pythonProject/PFAS-2.0/input_data/OECD_4000/step3_OECD_Class_0812.csv')
PFAA_precursors_labels = [

    (24, "n:2 fluorotelomer-based substances"),
    (23, "PASF-based substances"),
    (22, 'PolyFAC derivatives'),
    (21, "PFAK derivatives"),
    (20, "Sulfonyl chloride"),
    (19, "Acid chloride"),
    (18, "PolyFEAenes"),
    (17, "PolyFAenes"),
    (16, "PolyFEACs"),
    (15, "n:1 FTOHs"),
    (14, 'PolyFACs'),
    (13, "PFAenes"),
    (12, "SFAenes"),

    (11, 'PAECFs'),
    (10, "PFAKs"),
    (9, "PFALs"),
    (8, "PFAIs"),
    (7, "PFACs"),
    (6, "PASFs"),
    (5, "PACFs"),
    (4, "SFAs"),
    (3, "HFEs"),
    (2, "HFCs"),
    (1, "HFOs"),

]

PFAA_precursors_dict = {
    "HFOs": 1,
    "HFCs": 2,
    "HFCs, cyclic": 2,
    'HFEs': 3,
    'HFEs, cyclic': 3,
    'SFAs': 4,
    'PACFs': 5,
    'PACFs, cyclic': 5,
    'PASFs': 6,
    "PASFs, cyclic": 6,
    'PFACs': 7,
    'PFAIs': 8,
    'PFALs': 9,
    'PFAKs': 10,
    "PAECFs": 11,
    "PAECFs, cyclic": 11,
    'SFAenes': 12,
    'PFAenes': 13,
    'PFAenes, cyclic': 13,
    'PolyFACs': 14,
    'n:1 FTOHs, cyclic': 15,
    'n:1 FTOHs': 15,
    'PolyFEACs': 16,
    'PolyFEACs, cyclic': 16,
    "PolyFAenes": 17,
    "PolyFAenes, cyclic": 17,
    'PolyFEAenes': 18,
    'PolyFEAenes, cyclic': 18,
    'Acid chloride': 19,
    "Sulfonyl chloride": 20,
    'PFAK derivatives': 21,
    "PFAK derivatives, cyclic": 21,
    'PolyFAC derivatives': 22,
    'PolyFAC derivatives, cyclic': 22,
    'PASF-based substances': 23,
    "PASF-based substances, cyclic": 23,
    'n:2 fluorotelomer-based substances': 24,
    'n:2 fluorotelomer-based substances, cyclic': 24,



}
df_PFAA_subclass = df_tmap[(df_tmap['First_Class'] == 'PFAA precursors') | (df_tmap['First_Class'] == 'PFAA precursors, cyclic')]
list_main_smi = df_PFAA_subclass['RDKIT_SMILES'].tolist()
list_main_num = df_PFAA_subclass['Second_Class'].apply(lambda m: PFAA_precursors_dict[m]).tolist()



labels = []
for index, i in df_PFAA_subclass.iterrows():
    label = (
            i['SMILES'] +
            "__<h1>" + i['Second_Class'] + "</h1>"
    )
    labels.append(label)



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

faerun = Faerun(coords=False, view='front', title='PFAA_precursors_MAP')
faerun.add_scatter(
    "PFAA_precursors",
    {"x": x,
     "y": y,
     "c": [list_main_num],
     "labels": labels},
    point_scale=9,
    max_point_size=10,
    colormap='rainbow',
    has_legend=True,
    categorical=False,
    legend_labels=PFAA_precursors_labels,
    shader='sphere'

)
faerun.add_tree("PFAA_precursors_tree", {"from": s, "to": t}, point_helper="PFAA_precursors")

faerun.plot('Bert_PFAA_precursors_Map', template="smiles")
