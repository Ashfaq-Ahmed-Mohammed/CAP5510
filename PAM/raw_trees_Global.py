import pandas as pd
import matplotlib.pyplot as plt
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import io

def create_tree_from_csv(csv_path, algorithm_name, method='nj'):
    """
    Reads pairwise distances from CSV, builds a DistanceMatrix, 
    constructs a tree, and saves the plot.
    """
    print(f"--- Processing {algorithm_name} ({csv_path}) ---")
    
    # 1. Read Data
    df = pd.read_csv(csv_path)
    
    # 2. Extract unique IDs (Species names preferred for tree labels)
    # We combine species_1 and species_2 columns to get all unique species
    all_species = pd.concat([df['species_1'], df['species_2']]).unique()
    species_list = sorted(all_species.tolist())
    n = len(species_list)
    
    # Map species name to index for matrix construction
    spec_to_idx = {spec: i for i, spec in enumerate(species_list)}
    
    # 3. Initialize Lower Triangular Matrix
    # Biopython DistanceMatrix expects a lower triangular list of lists
    # matrix[i] has i+1 elements (from 0 to i)
    matrix = [[0.0] * (i + 1) for i in range(n)]
    
    # 4. Fill Matrix
    # We iterate through the dataframe and populate the matrix
    for _, row in df.iterrows():
        s1 = row['species_1']
        s2 = row['species_2']
        dist = float(row['distance'])
        
        idx1 = spec_to_idx[s1]
        idx2 = spec_to_idx[s2]
        
        # Ensure we fill lower triangle (row index >= col index)
        if idx1 >= idx2:
            matrix[idx1][idx2] = dist
        else:
            matrix[idx2][idx1] = dist

    # 5. Create Biopython DistanceMatrix Object
    dm = DistanceMatrix(names=species_list, matrix=matrix)
    
    # 6. Construct Tree
    constructor = DistanceTreeConstructor()
    if method == 'nj':
        tree = constructor.nj(dm)
    else:
        tree = constructor.upgma(dm)
        
    # 7. Better Visualization (ASCII + Image)
    print(f"\n{algorithm_name} Tree (ASCII):")
    Phylo.draw_ascii(tree)
    
    # Save Tree Object (Newick format)
    tree_file = f"tree_{algorithm_name.lower().replace(' ', '_')}.nwk"
    Phylo.write(tree, tree_file, "newick")
    print(f"Tree saved to {tree_file}")

    # Plotting
    plt.figure(figsize=(15, 10))
    
    # Remove internal node labels (confusing numbers) for cleaner plot
    for clade in tree.find_clades():
        if not clade.is_terminal():
            clade.name = None
            
    # Ladderize makes the tree look 'standard' (shortest branches up top)
    tree.ladderize()
    
    Phylo.draw(tree, do_show=False)
    plt.title(f"Phylogenetic Tree - {algorithm_name} ({method.upper()})")
    plt.savefig(f"plot_{algorithm_name.lower().replace(' ', '_')}.png")
    plt.close()
    print(f"Plot saved to plot_{algorithm_name.lower().replace(' ', '_')}.png\n")
    
    return tree

# --- MAIN PIPELINE ---
if __name__ == "__main__":
    # File paths (CHANGE THESE to your actual filenames)
    files = {
        "Global Alignment": "nw_full_results.csv"
    }
    
    trees = {}
    
    for name, filepath in files.items():
        try:
            # Using Neighbor Joining ('nj')
            t = create_tree_from_csv(filepath, name, method='nj')
            trees[name] = t
        except Exception as e:
            print(f"Error processing {name}: {e}")

    print("Pipeline Complete!")
