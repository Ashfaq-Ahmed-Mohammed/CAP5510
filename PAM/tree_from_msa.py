import pandas as pd
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from Bio import Phylo
import matplotlib.pyplot as plt

def calculate_msa_distance_matrix(msa_file, mapping_csv, output_csv="distances_msa.csv"):
    # 1. Load ID -> Species Mapping
    print("Loading Species Mapping...")
    df_map = pd.read_csv(mapping_csv)
    # We need a dict: {str(Protein_ID) : Species_Name}
    # Check which columns hold IDs. Assuming 'protein_id_1' and 'species_1'
    id_to_species = {}
    for _, row in df_map.iterrows():
        id_to_species[str(row['protein_id_1'])] = row['species_1']
        id_to_species[str(row['protein_id_2'])] = row['species_2']

    # 2. Load MSA
    print(f"Loading MSA: {msa_file}...")
    records = list(SeqIO.parse(msa_file, "fasta"))
    
    # Convert IDs to Species Names immediately
    # This ensures the Matrix and Tree use Species Names
    final_ids = []
    for rec in records:
        pid = str(rec.id)
        if pid in id_to_species:
            final_ids.append(id_to_species[pid])
        else:
            print(f"Warning: ID {pid} not found in mapping CSV. Keeping ID.")
            final_ids.append(pid)
            
    seqs = [str(rec.seq) for rec in records]
    n = len(seqs)
    length = len(seqs[0]) 
    
    print(f"Calculating Distance Matrix for {n} sequences...")
    
    matrix_data = [[0.0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(i + 1, n):
            s1 = seqs[i]
            s2 = seqs[j]
            
            valid_positions = 0
            diffs = 0
            
            for k in range(length):
                c1 = s1[k]
                c2 = s2[k]
                if c1 == '-' and c2 == '-':
                    continue
                valid_positions += 1
                if c1 != c2:
                    diffs += 1
            
            if valid_positions > 0:
                dist = diffs / valid_positions
            else:
                dist = 1.0
                
            matrix_data[i][j] = dist
            matrix_data[j][i] = dist

    # 3. Build Tree with SPECIES NAMES
    lower_tri_matrix = [matrix_data[i][:i+1] for i in range(n)]
    
    dm = DistanceMatrix(names=final_ids, matrix=lower_tri_matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    
    # 4. Save Tree
    Phylo.write(tree, "tree_msa.nwk", "newick")
    print("Tree saved as 'tree_msa.nwk' (with Species Names)")



def calculate_msa_distance_matrix_muscle(msa_file, mapping_csv, output_csv="distances_msa_muscle.csv"):
    # 1. Load ID -> Species Mapping
    print("Loading Species Mapping...")
    df_map = pd.read_csv(mapping_csv)
    # We need a dict: {str(Protein_ID) : Species_Name}
    # Check which columns hold IDs. Assuming 'protein_id_1' and 'species_1'
    id_to_species = {}
    for _, row in df_map.iterrows():
        id_to_species[str(row['protein_id_1'])] = row['species_1']
        id_to_species[str(row['protein_id_2'])] = row['species_2']

    # 2. Load MSA
    print(f"Loading MSA: {msa_file}...")
    records = list(SeqIO.parse(msa_file, "fasta"))
    
    # Convert IDs to Species Names immediately
    # This ensures the Matrix and Tree use Species Names
    final_ids = []
    for rec in records:
        pid = str(rec.id)
        if pid in id_to_species:
            final_ids.append(id_to_species[pid])
        else:
            print(f"Warning: ID {pid} not found in mapping CSV. Keeping ID.")
            final_ids.append(pid)
            
    seqs = [str(rec.seq) for rec in records]
    n = len(seqs)
    length = len(seqs[0]) 
    
    print(f"Calculating Distance Matrix for {n} sequences...")
    
    matrix_data = [[0.0] * n for _ in range(n)]
    
    for i in range(n):
        for j in range(i + 1, n):
            s1 = seqs[i]
            s2 = seqs[j]
            
            valid_positions = 0
            diffs = 0
            
            for k in range(length):
                c1 = s1[k]
                c2 = s2[k]
                if c1 == '-' and c2 == '-':
                    continue
                valid_positions += 1
                if c1 != c2:
                    diffs += 1
            
            if valid_positions > 0:
                dist = diffs / valid_positions
            else:
                dist = 1.0
                
            matrix_data[i][j] = dist
            matrix_data[j][i] = dist

    # 3. Build Tree with SPECIES NAMES
    lower_tri_matrix = [matrix_data[i][:i+1] for i in range(n)]
    
    dm = DistanceMatrix(names=final_ids, matrix=lower_tri_matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)
    
    # 4. Save Tree
    Phylo.write(tree, "tree_muscle_msa.nwk", "newick")
    print("Tree saved as 'tree_muscle_msa.nwk' (with Species Names)")

if __name__ == "__main__":
    # Pass your GLOBAL CSV as the mapping source
    calculate_msa_distance_matrix("msa_from_csv.fasta", "nw_full_results.csv")
    calculate_msa_distance_matrix_muscle("msa_muscle.fasta", "nw_full_results.csv")
    
