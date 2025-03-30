import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.BRICS import BRICSDecompose

ccdc_path = "data/database_CCDC.csv"
gan_path = "data/database_GAN.csv"
output_path = "data/filtered_coformers.csv"

MAX_MW = 500
MAX_LOGP = 5
MAX_SA = 3

def is_valid_smiles(smiles: str) -> bool:
    """ Проверяет, корректен ли SMILES (можно ли его распарсить в RDKit). """
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

def synthetic_accessibility(mol):
    """ Рассчитывает SA Score через разложение BRICS. """
    if mol is None:
        return 10

    try:
        num_fragments = len(list(BRICSDecompose(mol)))
    except:
        num_fragments = 5

    num_heavy_atoms = Chem.Lipinski.HeavyAtomCount(mol)
    sa_score = num_fragments + (num_heavy_atoms / 10)

    return min(sa_score, 10)

def calculate_molecular_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"MW": None, "LogP": None, "SA": None}

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    sa_score = synthetic_accessibility(mol)

    return {"MW": mw, "LogP": logp, "SA": sa_score}

def process_dataset(file_path: str, dataset_name: str) -> pd.DataFrame:
    print(f" Загружаем {dataset_name} из {file_path}...")

    df = pd.read_csv(file_path)
    total_molecules = len(df)

    df.columns = ["SMILES"] if df.shape[1] == 1 else df.columns

    df = df.drop_duplicates(subset="SMILES")
    after_dedup = len(df)
    print(f"   🛠 Убрано дубликатов: {total_molecules - after_dedup}")

    df["valid"] = df["SMILES"].apply(is_valid_smiles)
    df = df[df["valid"]].drop(columns=["valid"])  
    after_valid = len(df)
    print(f"Убрано невалидных SMILES: {after_dedup - after_valid}")

    print(f"Вычисляем молекулярные свойства для {len(df)} молекул...")
    df[["MW", "LogP", "SA"]] = df["SMILES"].apply(lambda x: pd.Series(calculate_molecular_properties(x)))

    df = df.dropna(subset=["MW", "LogP", "SA"])
    after_props = len(df)
    print(f"Убрано молекул без рассчитанных свойств: {after_valid - after_props}")

    df = df[(df["MW"] <= MAX_MW) & (df["LogP"] <= MAX_LOGP) & (df["SA"] <= MAX_SA)]
    after_filtering = len(df)
    print(f"Убрано по фильтрам (MW, LogP, SA): {after_props - after_filtering}")

    print(f"Очищенный {dataset_name}: {len(df)} молекул после фильтрации SA ≤ {MAX_SA}.")
    return df

df_ccdc = process_dataset(ccdc_path, "CCDC")
df_gan = process_dataset(gan_path, "GAN")

df_filtered = pd.concat([df_ccdc, df_gan], ignore_index=True)
print(f" Итоговое количество молекул после фильтрации SA ≤ {MAX_SA}: {len(df_filtered)}")

df_filtered.to_csv(output_path, index=False)
print(f"Файл с обновленными коформерами сохранен в {output_path}")
