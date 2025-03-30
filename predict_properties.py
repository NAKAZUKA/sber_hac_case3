import pandas as pd
from classifier import Classifier

input_file = "data/filtered_coformers.csv"
output_file = "data/predicted_coformers.csv"

df = pd.read_csv(input_file)

THEOPHYLLINE_SMILES = "CN1C=NC2=C1C(=O)N(C(=O)C2)C"

clf = Classifier()

results = []

for _, row in df.iterrows():
    smiles = row["SMILES"]
    
    try:
        predictions = clf.predict_properties(THEOPHYLLINE_SMILES, smiles, 
                                             properties=["unobstructed", "orthogonal_planes", "h_bond_bridging"])
        
        results.append({
            "SMILES": smiles,
            "MW": row["MW"],
            "LogP": row["LogP"],
            "SA": row["SA"],
            "Unobstructed": predictions.get("unobstructed", None),
            "Orthogonal_Planes": predictions.get("orthogonal_planes", None),
            "H_Bond_Bridging": predictions.get("h_bond_bridging", None)
        })
    except Exception as e:
        print(f"Ошибка при обработке {smiles}: {e}")

df_results = pd.DataFrame(results)

df_results.to_csv(output_file, index=False)

print(f"Файл с предсказанными механическими свойствами сохранён: {output_file}")
