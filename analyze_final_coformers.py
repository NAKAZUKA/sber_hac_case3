import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw, AllChem, DataStructs
from classifier import Classifier

final_file = "data/final_coformers.csv"
filtered_file = "data/filtered_coformers.csv"
output_dir = "sources"
os.makedirs(output_dir, exist_ok=True)

df_final = pd.read_csv(final_file)
df_filtered = pd.read_csv(filtered_file)

print(f"Оптимизированных молекул (final_coformers): {len(df_final)}")
print(f"Исходных молекул (filtered_coformers): {len(df_filtered)}")

clf = Classifier()

def median_probability(smiles):
    try:
        data = clf.create_clf_dataframe("CN1C=NC2=C1C(=O)N(C(=O)C2)C", smiles)
        p1 = clf.gbc_unobstructed.predict_proba(data[clf.features_unobstructed])[:,1][0]
        p2 = clf.gbc_orthogonal_planes.predict_proba(data[clf.features_orthogonal_planes])[:,1][0]
        p3 = clf.gbc_h_bond_bridging.predict_proba(data[clf.features_h_bond_bridging])[:,0][0]
        return round(np.median([p1, p2, p3]), 3)
    except:
        return np.nan
original_smiles_set = set(df_filtered["SMILES"])

def is_new(smiles):
    return smiles not in original_smiles_set

def novelty_tanimoto(original_smi, new_smi):
    mol1 = Chem.MolFromSmiles(original_smi)
    mol2 = Chem.MolFromSmiles(new_smi)
    if not mol1 or not mol2:
        return None
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, 2)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, 2)
    tanimoto = DataStructs.TanimotoSimilarity(fp1, fp2)
    return round(1 - tanimoto, 3)

df_final["Score"] = -abs(df_final["MW_Optimized"] - 250) \
                    - abs(df_final["LogP_Optimized"] - 2.5) \
                    - abs(df_final["SA"] - 2.0)

df_top10 = df_final.sort_values("Score", ascending=False).head(10).copy()

df_top10["MedianProb"] = df_top10["Optimized_SMILES"].apply(median_probability)

df_top10["IsNew"] = df_top10["Optimized_SMILES"].apply(is_new)

df_top10["Novelty"] = df_top10.apply(lambda row: novelty_tanimoto(row["SMILES"], row["Optimized_SMILES"]), axis=1)

df_top10.to_csv(os.path.join(output_dir, "top10_molecules.csv"), index=False)

plt.figure(figsize=(10, 6))
plt.hist(df_final["MW_Original"], bins=20, alpha=0.6, label="MW Original", color="steelblue")
plt.hist(df_final["MW_Optimized"], bins=20, alpha=0.6, label="MW Optimized", color="orange")
plt.xlabel("MW")
plt.ylabel("Частота")
plt.title("MW до/после оптимизации")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "hist_MW_comparison.png"))
plt.close()

plt.figure(figsize=(10, 6))
plt.hist(df_final["LogP_Original"], bins=20, alpha=0.6, label="LogP Original", color="cadetblue")
plt.hist(df_final["LogP_Optimized"], bins=20, alpha=0.6, label="LogP Optimized", color="coral")
plt.xlabel("LogP")
plt.ylabel("Частота")
plt.title("LogP до/после оптимизации")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "hist_LogP_comparison.png"))
plt.close()
print("Сохранена гистограмма LogP.")
try:
    from rdkit.Chem import Draw
    top_mols = []
    legends = []
    for idx, row in df_top10.iterrows():
        mol = Chem.MolFromSmiles(row["Optimized_SMILES"])
        if mol:
            legends.append(f"{idx+1}. Score={row['Score']:.2f}, Prob={row['MedianProb']:.2f}")
            top_mols.append(mol)

    if top_mols:
        img = Draw.MolsToGridImage(top_mols, legends=legends, molsPerRow=5, subImgSize=(200, 200))
        img.save(os.path.join(output_dir, "top10_molecules.png"))
except Exception as e:
    print("Не удалось визуализировать топ-10 молекул:", e)

print("Готово! См. top10_molecules.csv в папке sources/, а также графики и изображения.")
