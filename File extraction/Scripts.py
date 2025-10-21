
"""
fragmented_search.py
--------------------
Fragment a complex molecule and perform similarity searches on PubChem and ZINC15
for each fragment. Results are combined into a single dataset.
"""

import os
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, BRICS, rdMolDescriptors as Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from tqdm import tqdm
import urllib.parse

# ==========================
# CONFIGURATION
# ==========================
QUERY_FILE = "FK506_tacro_6536850 (3).sdf"  # or SMILES string
OUTPUT_FILE = "fragmented_pubchem_zinc_results.csv"

# ==========================
# LOAD QUERY MOLECULE
# ==========================
def load_smiles(input_path_or_smiles: str) -> str:
    """Load molecule from SMILES or SDF."""
    if os.path.isfile(input_path_or_smiles) and input_path_or_smiles.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(input_path_or_smiles)
        for mol in suppl:
            if mol is not None:
                return Chem.MolToSmiles(mol)
        raise ValueError("No valid molecule found in SDF.")
    else:
        return input_path_or_smiles.strip()

# ==========================
# FRAGMENTATION
# ==========================
def fragment_molecule(smiles: str, max_frags=6):
    """Generate Murcko scaffold + BRICS fragments."""
    mol = Chem.MolFromSmiles(smiles)
    fragments = set()

    # Murcko scaffold
    try:
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        fragments.add(Chem.MolToSmiles(scaffold))
    except Exception:
        pass

    # BRICS fragments
    try:
        brics_frags = BRICS.BRICSDecompose(mol)
        for f in list(brics_frags)[:max_frags]:
            fragments.add(f)
    except Exception:
        pass

    print(f"🧩 Generated {len(fragments)} fragments for searching.")
    return list(fragments)

# ==========================
# PUBCHEM SEARCH
# ==========================
def search_pubchem(smiles, threshold=85, max_results=20):
    encoded = urllib.parse.quote(smiles)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{encoded}/JSON?Threshold={threshold}"
    r = requests.get(url)
    if r.status_code != 200:
        return pd.DataFrame()

    try:
        cids = r.json()['IdentifierList']['CID'][:max_results]
        cid_str = ",".join(map(str, cids))
        prop_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/property/CanonicalSMILES,InChIKey,MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount/CSV"
        df = pd.read_csv(pd.compat.StringIO(requests.get(prop_url).text))
        df["Source"] = "PubChem"
        df["FragmentUsed"] = smiles
        return df
    except Exception:
        return pd.DataFrame()

# ==========================
# ZINC SEARCH
# ==========================
def search_zinc(smiles, threshold=0.8):
    encoded = urllib.parse.quote(smiles)
    url = "https://zinc15.docking.org/substances.txt"
    params = {"smiles": encoded, "maxsimilarity": threshold}
    r = requests.get(url, params=params)
    if r.status_code == 200 and r.text.strip():
        with open("zinc_tmp.txt", "w", encoding="utf-8") as f:
            f.write(r.text)
        try:
            df = pd.read_csv("zinc_tmp.txt", sep="\t")
            df["Source"] = "ZINC15"
            df["FragmentUsed"] = smiles
            return df
        except Exception:
            return pd.DataFrame()
    return pd.DataFrame()

# ==========================
# MAIN EXECUTION
# ==========================
if __name__ == "__main__":
    query_smiles = load_smiles(QUERY_FILE)
    print(f"✅ Original SMILES:\n{query_smiles}\n")

    fragments = fragment_molecule(query_smiles)

    all_results = []
    for frag in tqdm(fragments, desc="Searching fragments"):
        pubchem_df = search_pubchem(frag)
        zinc_df = search_zinc(frag)
        for df in [pubchem_df, zinc_df]:
            if not df.empty:
                all_results.append(df)

    if all_results:
        combined = pd.concat(all_results, ignore_index=True)
        combined.drop_duplicates(subset=["CanonicalSMILES"], inplace=True, ignore_index=True)
        combined.to_csv(OUTPUT_FILE, index=False)
        print(f"💾 Combined dataset saved: {OUTPUT_FILE} ({len(combined)} unique compounds)")
    else:
        print("⚠️ No similar compounds found for any fragments.")
