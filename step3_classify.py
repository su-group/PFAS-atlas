import pandas as pd
from classification_helper import classify_pfas_molecule




def get_classify(df_classify):
    df = df_classify.copy()
    df.reset_index(inplace=True)
    df['Classification'] = df['RDKIT_SMILES'].map(lambda x: classify_pfas_molecule(x))
    df_classification = df[['RDKIT_SMILES', 'Classification']].copy()
    df_classification['First_Class'] = df_classification['Classification'].map(lambda x: x[0])
    df_classification['Second_Class'] = df_classification['Classification'].map(lambda x: x[1])
    df_join = pd.merge(df_classify, df_classification, on='RDKIT_SMILES', how='left')
    df_join.drop(columns='Classification', inplace=True)
    return df_join

# if __name__ == "__main__":
#
#     df_classify = pd.read_csv('apply_PFAS/supp_hit/supp_hit_mhfp_0830.csv')
#     df_getclass = get_classify(df_classify)
#     df_getclass.to_csv('apply_PFAS/supp_hit/supp_hit_class_0830.csv', index_label=False, index=False)
