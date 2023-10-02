import numpy as np
import pandas as pd
import tmap as tm
from tqdm import tqdm
from mhfp.encoder import MHFPEncoder
from faerun import Faerun
from rxnfp.transformer_fingerprints import (
    RXNBERTFingerprintGenerator, get_default_model_and_tokenizer, generate_fingerprints
)


df_tmap = pd.read_csv('/home/ubuntu/pythonProject/PFAS-2.0/input_data/OECD_4000/step3_OECD_Class_0812.csv')
PolyFAAs_labels = [
    (6, "PolyFCA derivatives"),
    (5, "PolyFSA derivatives"),
    (4, "PolyFESAs"),
    (3, "PolyFECAs"),
    (2, "PolyFSAs"),
    (1, "PolyFCAs"),


]
PolyFAAs_dict = {
    "PolyFCAs": 1,
    'PolyFCAs, cyclic': 1,
    'PolyFSAs': 2,
    'PolyFECAs': 3,
    'PolyFESAs': 4,
    'PolyFSA derivatives': 5,
    'PolyFCA derivatives': 6,
    'PolyFCA derivatives, cyclic': 6,

}



df_poly_subclass = df_tmap[(df_tmap['First_Class'] == 'Polyfluoroalkyl acids') | (df_tmap['First_Class'] == 'Polyfluoroalkyl acids, cyclic')]
list_main_smi = df_poly_subclass['RDKIT_SMILES'].tolist()
list_main_num = df_poly_subclass['Second_Class'].apply(lambda m: PolyFAAs_dict[m]).tolist()



labels = []
for index, i in df_poly_subclass.iterrows():
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

faerun = Faerun(coords=False, view='front', title='PolyFAA_acids_MAP')
faerun.add_scatter(
    "Polyfluoroalkyl_acids",
    {"x": x,
     "y": y,
     "c": [list_main_num],
     "labels": labels},
    point_scale=10,
    max_point_size=10,
    colormap='rainbow',
    has_legend=True,
    categorical=False,
    legend_labels=PolyFAAs_labels,
    #series_title="dGsolv",
    shader='sphere'

)
faerun.add_tree("Polyfluoroalkyl_acids_tree", {"from": s, "to": t}, point_helper="Polyfluoroalkyl_acids")

faerun.plot('Bert_PolyFAAs_Map', template="smiles")
