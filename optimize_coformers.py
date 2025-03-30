import random
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, BRICS, rdMolTransforms
from deap import base, creator, tools, algorithms

input_file = "data/predicted_coformers.csv"
output_file = "data/optimized_coformers.csv"

df = pd.read_csv(input_file)

df_filtered = df[(df["Unobstructed"] == 1) & 
                 (df["Orthogonal_Planes"] == 1) & 
                 (df["H_Bond_Bridging"] == 0)].copy()

print(f"🔹 Найдено {len(df_filtered)} кандидатов для оптимизации.")

def fitness_function(individual):
    smiles = individual[0]
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return -100,

    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    
    try:
        sa = len(list(BRICS.BRICSDecompose(mol)))
    except:
        sa = 5
    
    score = 100 - abs(mw - 250) - abs(logp - 2.5) - abs(sa - 2.0)
    return score,

def mutate(ind):
    smiles = ind[0]
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return ind,
    
    fragments = list(BRICS.BRICSDecompose(mol))
    if len(fragments) < 2:
        return ind,
    
    fragment_mols = [Chem.MolFromSmiles(f) for f in fragments if Chem.MolFromSmiles(f)]
    if len(fragment_mols) < 2:
        return ind,
    
    new_mol = next(BRICS.BRICSBuild(fragment_mols), None)
    
    if new_mol:
        new_smiles = Chem.MolToSmiles(new_mol)
        ind[0] = new_smiles
    return ind,

def safe_crossover(ind1, ind2):
    if len(ind1) > 1 and len(ind2) > 1:
        tools.cxTwoPoint(ind1, ind2)
    else:
        tools.cxUniform(ind1, ind2, indpb=0.5)
    return ind1, ind2

creator.create("FitnessMax", base.Fitness, weights=(1.0,))
creator.create("Individual", list, fitness=creator.FitnessMax)

toolbox = base.Toolbox()
toolbox.register("attr_smiles", random.choice, df_filtered["SMILES"].tolist())
toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_smiles, n=1)
toolbox.register("population", tools.initRepeat, list, toolbox.individual)

toolbox.register("evaluate", fitness_function)
toolbox.register("mate", safe_crossover)
toolbox.register("mutate", mutate)
toolbox.register("select", tools.selTournament, tournsize=3)

optimized_smiles_list = []

for smiles in df_filtered["SMILES"]:
    print(f"Оптимизация для {smiles}...")
    
    population = toolbox.population(n=10)
    CXPB, MUTPB, NGEN = 0.5, 0.3, 10
    
    for gen in range(NGEN):
        offspring = algorithms.varAnd(population, toolbox, cxpb=CXPB, mutpb=MUTPB)
        fits = list(map(toolbox.evaluate, offspring))
        
        for ind, fit in zip(offspring, fits):
            ind.fitness.values = fit
        
        population[:] = toolbox.select(offspring, k=len(population))
    
    best_ind = tools.selBest(population, k=1)[0]
    optimized_smiles_list.append(best_ind[0])

df_filtered["Optimized_SMILES"] = [s if Chem.MolFromSmiles(s) else "INVALID" for s in optimized_smiles_list]
df_filtered.to_csv(output_file, index=False)

print(f"Оптимизированные молекулы сохранены в {output_file}")
