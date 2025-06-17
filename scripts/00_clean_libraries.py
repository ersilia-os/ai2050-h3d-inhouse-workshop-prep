import os
from rdkit import Chem
import pandas as pd

DATAPATH = "../data"

supplier = Chem.SDMolSupplier(os.path.join(DATAPATH, "libraries", "Anti-Infective-Library-18716.sdf"))
data = []
for mol in supplier:
    if mol is None:
        continue
    smiles = Chem.MolToSmiles(mol)
    name = mol.GetProp("IDNUMBER")
    inchikey = mol.GetProp("InChI Key")
    data.append({"name": name,"inchikey": inchikey, "smiles": smiles})

df = pd.DataFrame(data)
df.to_csv(os.path.join(DATAPATH, "libraries", "chemdiv_19k_antiinfective.csv"), index=False)

#For Ersilia
df = df[["smiles"]]
df_ = pd.read_csv(os.path.join(DATAPATH, "libraries", "chemdiv_100k_generalistic.csv"))
df_ = df_[["smiles"]]

df.to_csv(os.path.join(DATAPATH, "libraries", "antiinfective_smiles.csv"), index=False)
df_.to_csv(os.path.join(DATAPATH, "libraries", "generalistic_smiles.csv"), index=False)