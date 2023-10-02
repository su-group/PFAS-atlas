import pandas as pd
from mhfp.encoder import MHFPEncoder
from rdkit import Chem


def calculation_mhfp(df_list):
    df_mhfp = []
    df_fault = []
    mhfp_encoder = MHFPEncoder(1024)

    for index, i in df_list.iterrows():
        x = Chem.MolFromSmiles(i['RDKIT_SMILES'])
        #print(str(index) +'  '+i['RDKIT_SMILES'] +'   ',end ='')
        if not x:
            df_fault.append(i['RDKIT_SMILES'])
        else:
            df_mhfp.append([i['RDKIT_SMILES'], i['SMILES'], mhfp_encoder.encode_mol(x)])
    df_mhfps = pd.DataFrame(df_mhfp, columns=['RDKIT_SMILES', 'SMILES', 'MHFP'])
    return df_mhfps




def calculation_mhfp_stru(df_list):
    df_mhfp = []
    df_fault = []
    mhfp_encoder = MHFPEncoder(1024)

    for index, i in df_list.iterrows():
        x = Chem.MolFromSmiles(i['RDKIT_SMILES'])
        # print(str(index) +'  '+i['RDKIT_SMILES'] +'   ',end ='')
        if not x:
            df_fault.append(i['RDKIT_SMILES'])
        else:
            df_mhfp.append([i['RDKIT_SMILES'], i['SMILES'], i['Structure_function'], mhfp_encoder.encode_mol(x)])
            df_mhfps = pd.DataFrame(df_mhfp, columns=['RDKIT_SMILES', 'SMILES', 'Structure_function', 'MHFP'])
    return df_mhfps



# if __name__ == "__main__":
#
#     df = pd.read_csv('apply_PFAS/supp_hit/supp_hit_pro_0830.csv')
#     df_getmhfp = calculation_mhfp(df)
#     df_getmhfp.to_csv('apply_PFAS/supp_hit/supp_hit_mhfp_0830.csv', index_label=False, index=False)
