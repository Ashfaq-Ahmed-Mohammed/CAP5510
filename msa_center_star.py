import pandas as pd
from Bio import SeqIO

# -------------------------------------------------------------------
# HELPER: Insert Gaps into MSA
# -------------------------------------------------------------------
def insert_gap_column(msa_dict, index):
    """Inserts a gap '-' at 'index' for EVERY sequence currently in the MSA."""
    for seq_id in msa_dict:
        msa_dict[seq_id].insert(index, '-')

# -------------------------------------------------------------------
# CENTER-STAR MSA BUILDER (Using CSV Data)
# -------------------------------------------------------------------
def build_msa_from_csv(csv_path, sequences_file):
    # 1. Load Sequences to initialize the Center
    # We need the raw unaligned sequence for the center to start
    seq_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(sequences_file, "fasta")}
    
    print(f"Loading CSV: {csv_path}...")
    df = pd.read_csv(csv_path)
    
    # 2. Identify Center Sequence
    # Calculate sum of distances for each protein
    # Use protein_id_1 / protein_id_2
    dist_sums = {}
    
    # Optimization: Iterate once to find sums
    for _, row in df.iterrows():
        id1 = str(row['protein_id_1'])
        id2 = str(row['protein_id_2'])
        d = row['distance']
        dist_sums[id1] = dist_sums.get(id1, 0) + d
        dist_sums[id2] = dist_sums.get(id2, 0) + d
        
    center_id = min(dist_sums, key=dist_sums.get)
    print(f"Center Sequence: {center_id} (Min Dist Sum: {dist_sums[center_id]:.4f})")
    
    # 3. Initialize MSA with Center
    if center_id not in seq_dict:
        print(f"Error: Center ID {center_id} not found in FASTA file.")
        return
        
    msa = {center_id: list(seq_dict[center_id])}
    
    # 4. Filter CSV for rows involving the Center
    # We need alignments where one side is the center.
    # Note: The CSV stores unique pairs. Center could be col 1 OR col 2.
    
    # Create a subset for speed
    center_rows = df[(df['protein_id_1'].astype(str) == center_id) | (df['protein_id_2'].astype(str) == center_id)]
    
    processed_count = 0
    total_others = len(seq_dict) - 1
    
    # Get list of all other sequences we need to add
    other_ids = [sid for sid in seq_dict if sid != center_id]
    
    for other_id in other_ids:
        processed_count += 1
        if processed_count % 20 == 0:
            print(f"  Merging {processed_count}/{total_others}...")
            
        # Find the row for (Center, Other)
        # Row where (id1=Center AND id2=Other) OR (id1=Other AND id2=Center)
        
        row = center_rows[
            ((center_rows['protein_id_1'].astype(str) == center_id) & (center_rows['protein_id_2'].astype(str) == other_id)) |
            ((center_rows['protein_id_1'].astype(str) == other_id) & (center_rows['protein_id_2'].astype(str) == center_id))
        ]
        
        if row.empty:
            print(f"Warning: No alignment found for {center_id} vs {other_id}. Skipping.")
            continue
            
        row = row.iloc[0] # Take the series
        
        # extract alignment strings
        # We need to know which column corresponds to Center and which to Other
        if str(row['protein_id_1']) == center_id:
            aln_center = row['aln1']
            aln_other = row['aln2']
        else:
            aln_center = row['aln2']
            aln_other = row['aln1']
            
        # 5. Merge Logic (Same as before)
        current_msa_center = msa[center_id]
        new_row_for_other = []
        
        msa_idx = 0
        pair_idx = 0
        
        # Iterate through the pairwise alignment
        while pair_idx < len(aln_center):
            char_center = aln_center[pair_idx]
            char_other = aln_other[pair_idx]
            
            if char_center != '-':
                # Match/Mismatch: Map to existing MSA column
                # Advance MSA pointer past any existing gaps
                while msa_idx < len(current_msa_center) and current_msa_center[msa_idx] == '-':
                    new_row_for_other.append('-')
                    msa_idx += 1
                
                new_row_for_other.append(char_other)
                msa_idx += 1
                pair_idx += 1
            else:
                # Gap in Center (Insertion in Other): Needs new MSA column
                insert_gap_column(msa, msa_idx)
                
                # Refresh reference
                current_msa_center = msa[center_id]
                
                new_row_for_other.append(char_other)
                msa_idx += 1
                pair_idx += 1
                
        # Fill trailing gaps
        while msa_idx < len(current_msa_center):
            new_row_for_other.append('-')
            msa_idx += 1
            
        msa[other_id] = new_row_for_other

    # 6. Export
    final_msa_str = {k: "".join(v) for k, v in msa.items()}
    
    out_name = "msa_from_csv.fasta"
    with open(out_name, "w") as f:
        for sid, seq in final_msa_str.items():
            f.write(f">{sid}\n{seq}\n")
            
    print(f"Done! MSA saved to {out_name}")
    return final_msa_str

# -------------------------------------------------------------------
# RUN
# -------------------------------------------------------------------
if __name__ == "__main__":
    # ADJUST FILENAMES AS NEEDED
    CSV_FILE = "nw_full_results.csv"
    FASTA_FILE = "hemoglobin_209_species_final_fasta.fasta"
    
    build_msa_from_csv(CSV_FILE, FASTA_FILE)
