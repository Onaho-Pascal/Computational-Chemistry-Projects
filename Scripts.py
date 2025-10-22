from chembl_webresource_client.new_client import new_client
import pandas as pd
import requests
import os
from tqdm import tqdm

# ===============================
# CONFIGURATION
# ===============================
QUERY_SMILES = "C=CCC1/C=C(\\C)CC(C)CC(OC)C2OC(O)(C(=O)C(=O)N3CCCCC3C(=O)OC(/C(C)=C/C3CCC(O)C(OC)C3)C(C)C(O)CC1=O)C(C)CC2OC"
IMAGE_DIR = "chembl_images"
OUT_FILE = "chembl_full_results.csv"
os.makedirs(IMAGE_DIR, exist_ok=True)

print(f"✅ Query SMILES: {QUERY_SMILES}")
print("\n🔍 Searching ChEMBL for similar and substructure compounds...")

# ===============================
# STEP 1: SEARCH COMPOUNDS
# ===============================
similarity_search = new_client.similarity
substructure_search = new_client.substructure

# similarity > 70% Tanimoto
similar_hits = similarity_search.filter(smiles=QUERY_SMILES, similarity=70).only(
    ["molecule_chembl_id", "pref_name", "molecule_type", "canonical_smiles"]
)

# substructure search
substructure_hits = substructure_search.filter(smiles=QUERY_SMILES).only(
    ["molecule_chembl_id", "pref_name", "molecule_type", "canonical_smiles"]
)

results = list(similar_hits) + list(substructure_hits)
if not results:
    print("❌ No results found in ChEMBL.")
    exit()

df = pd.DataFrame(results).drop_duplicates(subset=["molecule_chembl_id"])
print(f"✅ Retrieved {len(df)} unique compounds.")

# ===============================
# STEP 2: DOWNLOAD IMAGES
# ===============================
print("\n🧩 Downloading 2D structure images...")
for chembl_id in tqdm(df["molecule_chembl_id"], desc="Downloading images"):
    img_url = f"https://www.ebi.ac.uk/chembl/api/data/image/{chembl_id}.svg"
    img_path = os.path.join(IMAGE_DIR, f"{chembl_id}.svg")
    try:
        r = requests.get(img_url, timeout=10)
        if r.status_code == 200:
            with open(img_path, "wb") as f:
                f.write(r.content)
        else:
            pass  # skip failed downloads
    except Exception:
        pass

# ===============================
# STEP 3: FETCH BIOACTIVITY DATA
# ===============================
print("\n🧪 Fetching bioactivity data (IC50, Ki, etc.)...")
bioactivity = new_client.activity

bio_data = []
for chembl_id in tqdm(df["molecule_chembl_id"], desc="Fetching bioactivity"):
    acts = bioactivity.filter(molecule_chembl_id=chembl_id).only(
        ["standard_type", "standard_value", "standard_units", "target_chembl_id", "target_pref_name"]
    )
    for a in acts:
        a["molecule_chembl_id"] = chembl_id
        bio_data.append(a)

if bio_data:
    bio_df = pd.DataFrame(bio_data)
    merged = pd.merge(df, bio_df, on="molecule_chembl_id", how="left")
else:
    merged = df
    print("⚠️ No bioactivity data found for these molecules.")

# ===============================
# STEP 4: SAVE RESULTS
# ===============================
merged.to_csv(OUT_FILE, index=False)
print(f"\n💾 Full results (structures + bioactivity) saved to: {OUT_FILE}")
print(f"🖼️ Images saved in folder: {IMAGE_DIR}/")
