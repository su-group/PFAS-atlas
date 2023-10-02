import pandas as pd
from rdkit_helper import rdkit_smiles_from_input_smiles, mol_from_smiles, rdkit_smiles_from_mol
import numpy as np


# data process(PFAS claass）
def prepare_rdkit_smiles(df):
    # STEP ONE Remove PFASs with NO SMILES
    df = df[['SMILES']].replace(['', '-'], np.nan)
    df.dropna(inplace=True)
    # STEP TWO Convert Input SMILES to RDKit SMILES (including removing salts and extract PFAS units from complexes)
    df['RDKIT_SMILES'] = df['SMILES'].map(lambda x: rdkit_smiles_from_input_smiles(x))
    # STEP THREE Remove PFASs with INVALID SMILES (CANNOT CONVERT TO VALID RDKIT SMILES)
    df.dropna(axis=0, inplace=True)
    # STEP FOUR Remove Duplicates
    df.drop_duplicates(subset=['RDKIT_SMILES'], inplace=True)
    # STEP FIVE Remove Complexes containing more than one type of PFAS units
    df = df[~df['RDKIT_SMILES'].str.contains("\.")]
    # STEP Six Remove Isomeric SMILES (To simplify classification)
    df = df[~df['SMILES'].str.contains('@')]  # Remove isomeric smiles
    df['RDKIT_SMILES_help'] = df['RDKIT_SMILES'].apply(lambda x: mol_from_smiles(x))
    df['RDKIT_SMILES'] = df['RDKIT_SMILES_help'].apply(lambda x: rdkit_smiles_from_mol(x))
    return df


# data process(Structure_function）
def prepare_rdkit_stru(df):
    # STEP ONE Remove PFASs with NO SMILES
    df = df[['SMILES', 'Structure_function']]
    df[['SMILES']].replace(['', '-'], np.nan)
    df.dropna(inplace=True)
    # STEP TWO Convert Input SMILES to RDKit SMILES (including removing salts and extract PFAS units from complexes)
    df['RDKIT_SMILES'] = df['SMILES'].map(lambda x: rdkit_smiles_from_input_smiles(x))
    # STEP THREE Remove PFASs with INVALID SMILES (CANNOT CONVERT TO VALID RDKIT SMILES)
    df.dropna(axis=0, inplace=True)
    # STEP FOUR Remove Duplicates
    df.drop_duplicates(subset=['RDKIT_SMILES'], inplace=True)
    # STEP FIVE Remove Complexes containing more than one type of PFAS units
    df = df[~df['RDKIT_SMILES'].str.contains("\.")]
    # STEP Six Remove Isomeric SMILES (To simplify classification)
    df = df[~df['SMILES'].str.contains('@')]  # Remove isomeric smiles
    df['RDKIT_SMILES_help'] = df['RDKIT_SMILES'].apply(lambda x: mol_from_smiles(x))
    df['RDKIT_SMILES'] = df['RDKIT_SMILES_help'].apply(lambda x: rdkit_smiles_from_mol(x))
    return df


#if __name__ == "__main__":
    # # create rdkit smiles
    #
    # df = pd.read_csv('/home/dell/PycharmProjects/PFAS-2.0/apply_PFAS/Neurotoxicity_pro.csv')
    # # df['Structure_function'].fillna('no activity data', inplace=True)
    # df_smiles = prepare_rdkit_smiles(df)
    # df_smiles.to_csv('/home/dell/PycharmProjects/PFAS-2.0/apply_PFAS/Neurotoxicity_pro.csv', index_label=False, index=False)



