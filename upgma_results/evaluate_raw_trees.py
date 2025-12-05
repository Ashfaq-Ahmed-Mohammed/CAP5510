import pandas as pd
import numpy as np
from Bio import Phylo

# --- 1. Helper to get Taxonomy ---
def get_taxonomy_mapping(csv_file):
    df = pd.read_csv(csv_file)
    mapping = {}
    for _, row in df.iterrows():
        mapping[row['species_1']] = row['superorder_1']
        mapping[row['species_2']] = row['superorder_2']
    return mapping

# --- 2. IMPROVED Intrinsic Metric: Robinson-Foulds (Manual Split Calculation) ---
def get_clade_splits(tree):
    """
    Generates a set of splits (bipartitions) for a tree.
    Each split is represented by a frozenset of taxa names on one side of the branch.
    """
    splits = set()
    taxa = sorted([t.name for t in tree.get_terminals()])
    all_taxa = frozenset(taxa)
    
    for clade in tree.find_clades():
        # Get all terminals under this clade
        clade_taxa = frozenset([t.name for t in clade.get_terminals()])
        
        # A split is defined by the set of taxa (normalized to be the smaller side)
        # Ignore trivial splits (leaf nodes or full tree)
        if len(clade_taxa) > 1 and len(clade_taxa) < len(all_taxa):
            # Normalize: always take the smaller set to represent the split
            if len(clade_taxa) > len(all_taxa) / 2:
                splits.add(all_taxa - clade_taxa)
            else:
                splits.add(clade_taxa)
    return splits

def calculate_rf_distance(tree_file_1, tree_file_2):
    t1 = Phylo.read(tree_file_1, "newick")
    t2 = Phylo.read(tree_file_2, "newick")
    
    splits1 = get_clade_splits(t1)
    splits2 = get_clade_splits(t2)
    
    # RF Distance = (Splits in A but not B) + (Splits in B but not A)
    diff1 = len(splits1 - splits2)
    diff2 = len(splits2 - splits1)
    
    return diff1 + diff2

# --- 3. ROBUST Extrinsic Metric: F1-Score for Clades ---
def evaluate_clustering_robust(tree_file, tax_mapping):
    """
    Instead of finding one MRCA (fragile), we evaluate every node in the tree.
    For each Superorder, we find the 'Best Matching Clade' (highest F1 Score).
    """
    tree = Phylo.read(tree_file, "newick")
    superorders = set(tax_mapping.values())
    
    print(f"\n--- Evaluation: {tree_file} ---")
    print(f"{'Superorder':<20} | {'Precision':<10} | {'Recall':<10} | {'F1 Score':<10}")
    print("-" * 60)
    
    avg_f1 = 0
    count = 0
    
    for so in superorders:
        target_species = {sp for sp, tax in tax_mapping.items() if tax == so}
        total_targets = len(target_species)
        
        best_f1 = 0.0
        best_prec = 0.0
        best_rec = 0.0
        
        # Scan every clade in the tree to find the "Home" of this Superorder
        for clade in tree.find_clades():
            terminals = {t.name for t in clade.get_terminals()}
            if len(terminals) < 2: continue # Skip leaves
            
            # Intersection: How many targets are in this clade?
            tp = len(terminals.intersection(target_species)) # True Positives
            fp = len(terminals) - tp                         # False Positives
            fn = total_targets - tp                          # False Negatives
            
            if tp == 0: continue
            
            precision = tp / (tp + fp)
            recall = tp / (tp + fn)
            
            if (precision + recall) > 0:
                f1 = 2 * (precision * recall) / (precision + recall)
            else:
                f1 = 0
                
            if f1 > best_f1:
                best_f1 = f1
                best_prec = precision
                best_rec = recall
        
        print(f"{so:<20} | {best_prec:.2f}       | {best_rec:.2f}       | {best_f1:.2f}")
        avg_f1 += best_f1
        count += 1
        
    final_score = avg_f1 / count
    return final_score

# --- MAIN ---
if __name__ == "__main__":
    try:
        taxonomy = get_taxonomy_mapping("../nw_full_results.csv")
    except:
        # Placeholder for testing if file missing
        print("Warning: CSV not found. Ensure 'distances_global.csv' is present.")
        exit()

    files = {
        "Global": "nw_final_tree.nwk",
        "Local": "sw_final_tree.nwk",
        "SemiGlobal": "sg_final_tree.nwk"
    }
    
    scores = {}
    for name, f in files.items():
        try:
            s = evaluate_clustering_robust(f, taxonomy)
            scores[name] = s
            print(f">> {name} Macro-F1 Score: {s:.3f}")
        except Exception as e:
            print(f"Error {name}: {e}")

    
