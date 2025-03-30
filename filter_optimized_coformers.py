import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

input_file = "data/optimized_coformers.csv"
output_file = "data/final_coformers.csv"

df = pd.read_csv(input_file)
print(f"🔹 Загружено {len(df)} оптимизированных молекул.")

def is_valid_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

filtered_molecules = []
for _, row in df.iterrows():
    original_smiles = row["SMILES"]
    optimized_smiles = row["Optimized_SMILES"]
    
    if not is_valid_smiles(optimized_smiles):
        continue
    mol_original = Chem.MolFromSmiles(original_smiles)
    mol_optimized = Chem.MolFromSmiles(optimized_smiles)
    
    mw_orig, logp_orig, sa_orig = Descriptors.MolWt(mol_original), Descriptors.MolLogP(mol_original), row["SA"]
    mw_opt, logp_opt = Descriptors.MolWt(mol_optimized), Descriptors.MolLogP(mol_optimized)
    
    if abs(mw_opt - mw_orig) <= 50 and (0 <= logp_opt <= 5) and (sa_orig <= 3):
        filtered_molecules.append({
            "SMILES": original_smiles,
            "Optimized_SMILES": optimized_smiles,
            "MW_Original": mw_orig,
            "MW_Optimized": mw_opt,
            "LogP_Original": logp_orig,
            "LogP_Optimized": logp_opt,
            "SA": sa_orig
        })

df_final = pd.DataFrame(filtered_molecules)

df_final.to_csv(output_file, index=False)
print(f"Финальные отобранные молекулы сохранены в {output_file}")
print(f"Осталось {len(df_final)} молекул после фильтрации.")
