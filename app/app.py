import os
import streamlit as st
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image
import io

st.set_page_config(
    page_title="Ersilia Model Hub - Merge Datasets",
    page_icon=":molecule:",
    layout="wide",
    initial_sidebar_state="expanded",
)

ROOT = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(ROOT, "..", "data")

MAX_COLUMNS = 10

PRIMARY_MODEL_IDS = [
    "eos78ao",
    "eos2db3",
]

AUXILIARY_MODEL_IDS = [
    "eosXXX",
    "eosYYY"
]

st.title("Merge datasets from the Ersilia Model Hub")
st.markdown("This is a simple app to merge datasets from the Ersilia Model Hub. Feel free to use Excel or any specialized cheminformatics tool for more advanced usage")

st.sidebar.title("Upload datasets as `CSV`")

st.sidebar.header("Input file (SMILES)")

st.sidebar.header("Out")

st.sidebar.header("Dataset 3")

st.sidebar.header("Dataset 4")

st.sidebar.header("Dataset 5")

st.cache_data()
def upload_data():
    example_dir = os.path.join(DATA_DIR, "example")
    df = pd.read_csv(os.path.join(example_dir, "input.csv"))
    columns = df.columns.tolist()
    if "input" in columns:
        df = df.rename(columns={"input": "smiles"})
    for fn in os.listdir(example_dir):
        if not fn.startswith("eos") and fn.endswith(".csv"):
            continue
        model_id = fn.split(".")[0]
        df_ = pd.read_csv(os.path.join(example_dir, fn))
        columns = df_.columns.tolist()
        if "input" in columns:
            df_ = df_.drop(columns=["input"])
        if "key" in columns:
            df_ = df_.drop(columns=["key"])
        columns = df_.columns.tolist()
        for col in columns:
            if "drugbank_approved_percentile" in col:
                df_ = df_.drop(columns=[col])
        rename = dict((k, f"{k} - {model_id}") for k in columns)
        df_ = df_.rename(columns=rename)
        df = pd.concat([df, df_], axis=1)
    return df

df = upload_data()

columns = df.columns.tolist()
columns = [col for col in columns if col != "smiles"]

selected_columns = st.multiselect("Selected columns", options=columns)

if len(selected_columns) > MAX_COLUMNS:
    st.warning(f"You can only select up to {MAX_COLUMNS} columns for your selection. The first {MAX_COLUMNS} columns are kept for downstream analysis.")


st.cache_data()
def molecule_dictionary(df):
    mol_dict = {}
    for smiles in df["smiles"].tolist():
        mol = Chem.MolFromSmiles(smiles)
        mol_dict[smiles] = mol
    return mol_dict

mol_dict = molecule_dictionary(df)

df_ = df.copy()
df_ = df_[["smiles"] + selected_columns]

cols = st.columns(int(MAX_COLUMNS/2))
for i, column in enumerate(selected_columns[:int(MAX_COLUMNS/2)]):
    min_val, max_val = np.min(df[column]), np.max(df[column])
    col = cols[i]
    min_val, max_val = col.slider(column, min_val, max_val, (min_val, max_val))
    df_ = df_[df_[column] >= min_val]
    df_ = df_[df_[column] <= max_val]

for i, column in enumerate(selected_columns[int(MAX_COLUMNS/2):]):
    min_val, max_val = np.min(df[column]), np.max(df[column])
    col = cols[i]
    min_val, max_val = col.slider(column, min_val, max_val, (min_val, max_val))
    df_ = df_[df_[column] >= min_val]
    df_ = df_[df_[column] <= max_val]


st.dataframe(df_)

st.header("Molecule grid")

smiles_list = df_["smiles"].tolist()
mols = [mol_dict[smiles] for smiles in smiles_list]

def mols_to_grid_image(mols, molsPerRow=5, subImgSize=(300,300)):
    img = Draw.MolsToGridImage(
        mols, 
        molsPerRow=molsPerRow, 
        subImgSize=subImgSize, 
        legends=[str(i) for i in range(len(mols))]
    )
    return img

if mols:
    grid_img = mols_to_grid_image(mols[:25]) 
    buf = io.BytesIO()
    grid_img.save(buf, format="PNG")
    st.image(buf.getvalue())
else:
    st.info("No molecules to display.")