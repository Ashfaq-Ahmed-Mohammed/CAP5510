# unified_phylogenetic_pipeline.py

import csv
import json
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple

# All 26 amino acid codes (standard 20 + ambiguous codes)
AMINO_ACIDS = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V',
               'B', 'Z', 'J', 'X', 'U', 'O']

# Map ambiguous codes to standard amino acids
AMBIGUOUS_MAPPING = {
    'B': ['D', 'N'],
    'Z': ['E', 'Q'],
    'J': ['L', 'I'],
    'X': ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
          'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'],
    'U': ['C'],
    'O': ['K'],
    '*': [],
    '-': [],
}

# Standard 20 amino acids for profile scoring
STANDARD_AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 
               'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']


class Profile:
    """Represents a profile (frequency table) for aligned sequences"""
    
    def __init__(self, aligned_sequences):
        """
        Args:
            aligned_sequences: List of aligned sequence strings
        """
        self.alignment = aligned_sequences
        self.length = len(aligned_sequences[0]) if aligned_sequences else 0
        self.profile = self.calculate_profile()
    
    def calculate_profile(self):
        """Create frequency table for each position"""
        profile = []
        
        for pos in range(self.length):
            position_counts = {aa: 0.0 for aa in STANDARD_AA}
            total = 0.0
            
            for seq in self.alignment:
                if pos < len(seq):
                    aa = seq[pos].upper()
                    
                    if aa == '-' or aa == '*':
                        continue
                    
                    if aa in AMBIGUOUS_MAPPING:
                        mapped = AMBIGUOUS_MAPPING[aa]
                        if mapped:
                            for mapped_aa in mapped:
                                if mapped_aa in position_counts:
                                    position_counts[mapped_aa] += 1.0 / len(mapped)
                                    total += 1.0 / len(mapped)
                    elif aa in position_counts:
                        position_counts[aa] += 1.0
                        total += 1.0
            
            position_freq = {}
            for aa in STANDARD_AA:
                position_freq[aa] = position_counts[aa] / total if total > 0 else 0.0
            
            profile.append(position_freq)
        
        return profile
    
    def score_amino_acid(self, position, amino_acid, blosum62):
        """Score amino acid against profile at position"""
        score = 0.0
        aa = amino_acid.upper()
        
        if aa == '-' or aa == '*':
            return -4
        
        if aa in AMBIGUOUS_MAPPING:
            mapped = AMBIGUOUS_MAPPING[aa]
            if mapped:
                partial_sum = 0.0
                for mapped_aa in mapped:
                    if mapped_aa in STANDARD_AA:
                        for other_aa in STANDARD_AA:
                            freq = self.profile[position].get(other_aa, 0)
                            blosum_score = blosum62.get((mapped_aa, other_aa), -4)
                            partial_sum += freq * blosum_score
                score = partial_sum / len(mapped)
            else:
                score = -4
        elif aa in STANDARD_AA:
            for other_aa in STANDARD_AA:
                freq = self.profile[position].get(other_aa, 0)
                blosum_score = blosum62.get((aa, other_aa), -4)
                score += freq * blosum_score
        else:
            score = -4
        
        return score


class UnifiedPhylogeneticPipeline:
    """
    Unified pipeline that takes pairwise alignment output and builds phylogenetic trees.
    Works with NW, SW, SG alignment outputs - adapts to the column structure provided.
    
    PHASE 1: Load scores → Convert to distances
    PHASE 2: Build guide tree (UPGMA)
    PHASE 3&4: Progressive MSA with profiles
    PHASE 5-6: Final tree from MSA distances
    """
    
    def __init__(self, algorithm_name: str, blosum62: Dict = None):
        """
        Initialize pipeline
        
        Args:
            algorithm_name: Name of algorithm (e.g., "NW", "SW", "SG")
            blosum62: Dictionary of BLOSUM62 matrix (optional)
        """
        self.algorithm_name = algorithm_name
        self.sequence_names = []
        self.num_sequences = 0
        self.protein_ids = {}
        self.species_info = {}
        
        # Storage for outputs at each stage
        self.raw_scores_dict = {}
        self.distance_matrix = None
        self.guide_tree = None
        self.msa = None
        self.msa_profiles = None
        self.final_tree = None
        self.alignments = {}
        
        # BLOSUM62 scoring matrix
        if blosum62 is None:
            self.blosum62 = self._get_default_blosum62()
        else:
            self.blosum62 = blosum62
        
        print(f"[PIPELINE] Initialized {algorithm_name} pipeline")
    
    def _get_default_blosum62(self) -> Dict:
        """Get default BLOSUM62 scoring matrix (full 20x20)"""
        blosum62_full = {
            ('A', 'A'): 4, ('A', 'R'): -1, ('A', 'N'): -2, ('A', 'D'): -2,
            ('A', 'C'): 0, ('A', 'Q'): -1, ('A', 'E'): -1, ('A', 'G'): 0,
            ('A', 'H'): -2, ('A', 'I'): -1, ('A', 'L'): -1, ('A', 'K'): -1,
            ('A', 'M'): -1, ('A', 'F'): -2, ('A', 'P'): -1, ('A', 'S'): 1,
            ('A', 'T'): 0, ('A', 'W'): -3, ('A', 'Y'): -2, ('A', 'V'): 0,
            
            ('R', 'A'): -1, ('R', 'R'): 5, ('R', 'N'): 0, ('R', 'D'): -2,
            ('R', 'C'): -3, ('R', 'Q'): 1, ('R', 'E'): 0, ('R', 'G'): -2,
            ('R', 'H'): 0, ('R', 'I'): -3, ('R', 'L'): -2, ('R', 'K'): 2,
            ('R', 'M'): -1, ('R', 'F'): -3, ('R', 'P'): -2, ('R', 'S'): -1,
            ('R', 'T'): -1, ('R', 'W'): -3, ('R', 'Y'): -2, ('R', 'V'): -3,
            
            ('N', 'A'): -2, ('N', 'R'): 0, ('N', 'N'): 6, ('N', 'D'): 1,
            ('N', 'C'): -3, ('N', 'Q'): 0, ('N', 'E'): 0, ('N', 'G'): 0,
            ('N', 'H'): 1, ('N', 'I'): -3, ('N', 'L'): -3, ('N', 'K'): 0,
            ('N', 'M'): -2, ('N', 'F'): -3, ('N', 'P'): -2, ('N', 'S'): 1,
            ('N', 'T'): 0, ('N', 'W'): -4, ('N', 'Y'): -2, ('N', 'V'): -3,
            
            ('D', 'A'): -2, ('D', 'R'): -2, ('D', 'N'): 1, ('D', 'D'): 6,
            ('D', 'C'): -3, ('D', 'Q'): 0, ('D', 'E'): 2, ('D', 'G'): -1,
            ('D', 'H'): -1, ('D', 'I'): -3, ('D', 'L'): -4, ('D', 'K'): -1,
            ('D', 'M'): -3, ('D', 'F'): -3, ('D', 'P'): -1, ('D', 'S'): 0,
            ('D', 'T'): -1, ('D', 'W'): -4, ('D', 'Y'): -3, ('D', 'V'): -3,
            
            ('C', 'A'): 0, ('C', 'R'): -3, ('C', 'N'): -3, ('C', 'D'): -3,
            ('C', 'C'): 9, ('C', 'Q'): -3, ('C', 'E'): -4, ('C', 'G'): -3,
            ('C', 'H'): -3, ('C', 'I'): -1, ('C', 'L'): -1, ('C', 'K'): -3,
            ('C', 'M'): -1, ('C', 'F'): -2, ('C', 'P'): -3, ('C', 'S'): -1,
            ('C', 'T'): -1, ('C', 'W'): -2, ('C', 'Y'): -2, ('C', 'V'): -1,
            
            ('Q', 'A'): -1, ('Q', 'R'): 1, ('Q', 'N'): 0, ('Q', 'D'): 0,
            ('Q', 'C'): -3, ('Q', 'Q'): 5, ('Q', 'E'): 2, ('Q', 'G'): -2,
            ('Q', 'H'): 0, ('Q', 'I'): -3, ('Q', 'L'): -2, ('Q', 'K'): 1,
            ('Q', 'M'): 0, ('Q', 'F'): -3, ('Q', 'P'): -1, ('Q', 'S'): 0,
            ('Q', 'T'): -1, ('Q', 'W'): -2, ('Q', 'Y'): -1, ('Q', 'V'): -2,
            
            ('E', 'A'): -1, ('E', 'R'): 0, ('E', 'N'): 0, ('E', 'D'): 2,
            ('E', 'C'): -4, ('E', 'Q'): 2, ('E', 'E'): 5, ('E', 'G'): -2,
            ('E', 'H'): 0, ('E', 'I'): -3, ('E', 'L'): -3, ('E', 'K'): 1,
            ('E', 'M'): -2, ('E', 'F'): -3, ('E', 'P'): -1, ('E', 'S'): 0,
            ('E', 'T'): -1, ('E', 'W'): -3, ('E', 'Y'): -2, ('E', 'V'): -2,
            
            ('G', 'A'): 0, ('G', 'R'): -2, ('G', 'N'): 0, ('G', 'D'): -1,
            ('G', 'C'): -3, ('G', 'Q'): -2, ('G', 'E'): -2, ('G', 'G'): 6,
            ('G', 'H'): -2, ('G', 'I'): -4, ('G', 'L'): -4, ('G', 'K'): -2,
            ('G', 'M'): -3, ('G', 'F'): -3, ('G', 'P'): -2, ('G', 'S'): 0,
            ('G', 'T'): -2, ('G', 'W'): -2, ('G', 'Y'): -3, ('G', 'V'): -3,
            
            ('H', 'A'): -2, ('H', 'R'): 0, ('H', 'N'): 1, ('H', 'D'): -1,
            ('H', 'C'): -3, ('H', 'Q'): 0, ('H', 'E'): 0, ('H', 'G'): -2,
            ('H', 'H'): 8, ('H', 'I'): -3, ('H', 'L'): -3, ('H', 'K'): -1,
            ('H', 'M'): -2, ('H', 'F'): -1, ('H', 'P'): -2, ('H', 'S'): -1,
            ('H', 'T'): -2, ('H', 'W'): -2, ('H', 'Y'): 2, ('H', 'V'): -3,
            
            ('I', 'A'): -1, ('I', 'R'): -3, ('I', 'N'): -3, ('I', 'D'): -3,
            ('I', 'C'): -1, ('I', 'Q'): -3, ('I', 'E'): -3, ('I', 'G'): -4,
            ('I', 'H'): -3, ('I', 'I'): 4, ('I', 'L'): 2, ('I', 'K'): -3,
            ('I', 'M'): 1, ('I', 'F'): 0, ('I', 'P'): -3, ('I', 'S'): -2,
            ('I', 'T'): -1, ('I', 'W'): -3, ('I', 'Y'): -1, ('I', 'V'): 3,
            
            ('L', 'A'): -1, ('L', 'R'): -2, ('L', 'N'): -3, ('L', 'D'): -4,
            ('L', 'C'): -1, ('L', 'Q'): -2, ('L', 'E'): -3, ('L', 'G'): -4,
            ('L', 'H'): -3, ('L', 'I'): 2, ('L', 'L'): 4, ('L', 'K'): -2,
            ('L', 'M'): 2, ('L', 'F'): 0, ('L', 'P'): -3, ('L', 'S'): -2,
            ('L', 'T'): -1, ('L', 'W'): -2, ('L', 'Y'): -1, ('L', 'V'): 1,
            
            ('K', 'A'): -1, ('K', 'R'): 2, ('K', 'N'): 0, ('K', 'D'): -1,
            ('K', 'C'): -3, ('K', 'Q'): 1, ('K', 'E'): 1, ('K', 'G'): -2,
            ('K', 'H'): -1, ('K', 'I'): -3, ('K', 'L'): -2, ('K', 'K'): 5,
            ('K', 'M'): -1, ('K', 'F'): -3, ('K', 'P'): -1, ('K', 'S'): 0,
            ('K', 'T'): -1, ('K', 'W'): -3, ('K', 'Y'): -2, ('K', 'V'): -2,
            
            ('M', 'A'): -1, ('M', 'R'): -1, ('M', 'N'): -2, ('M', 'D'): -3,
            ('M', 'C'): -1, ('M', 'Q'): 0, ('M', 'E'): -2, ('M', 'G'): -3,
            ('M', 'H'): -2, ('M', 'I'): 1, ('M', 'L'): 2, ('M', 'K'): -1,
            ('M', 'M'): 5, ('M', 'F'): 0, ('M', 'P'): -2, ('M', 'S'): -1,
            ('M', 'T'): -1, ('M', 'W'): -1, ('M', 'Y'): -1, ('M', 'V'): 1,
            
            ('F', 'A'): -2, ('F', 'R'): -3, ('F', 'N'): -3, ('F', 'D'): -3,
            ('F', 'C'): -2, ('F', 'Q'): -3, ('F', 'E'): -3, ('F', 'G'): -3,
            ('F', 'H'): -1, ('F', 'I'): 0, ('F', 'L'): 0, ('F', 'K'): -3,
            ('F', 'M'): 0, ('F', 'F'): 6, ('F', 'P'): -4, ('F', 'S'): -2,
            ('F', 'T'): -2, ('F', 'W'): 1, ('F', 'Y'): 3, ('F', 'V'): -1,
            
            ('P', 'A'): -1, ('P', 'R'): -2, ('P', 'N'): -2, ('P', 'D'): -1,
            ('P', 'C'): -3, ('P', 'Q'): -1, ('P', 'E'): -1, ('P', 'G'): -2,
            ('P', 'H'): -2, ('P', 'I'): -3, ('P', 'L'): -3, ('P', 'K'): -1,
            ('P', 'M'): -2, ('P', 'F'): -4, ('P', 'P'): 7, ('P', 'S'): -1,
            ('P', 'T'): -1, ('P', 'W'): -4, ('P', 'Y'): -3, ('P', 'V'): -2,
            
            ('S', 'A'): 1, ('S', 'R'): -1, ('S', 'N'): 1, ('S', 'D'): 0,
            ('S', 'C'): -1, ('S', 'Q'): 0, ('S', 'E'): 0, ('S', 'G'): 0,
            ('S', 'H'): -1, ('S', 'I'): -2, ('S', 'L'): -2, ('S', 'K'): 0,
            ('S', 'M'): -1, ('S', 'F'): -2, ('S', 'P'): -1, ('S', 'S'): 4,
            ('S', 'T'): 1, ('S', 'W'): -3, ('S', 'Y'): -2, ('S', 'V'): -2,
            
            ('T', 'A'): 0, ('T', 'R'): -1, ('T', 'N'): 0, ('T', 'D'): -1,
            ('T', 'C'): -1, ('T', 'Q'): -1, ('T', 'E'): -1, ('T', 'G'): -2,
            ('T', 'H'): -2, ('T', 'I'): -1, ('T', 'L'): -1, ('T', 'K'): -1,
            ('T', 'M'): -1, ('T', 'F'): -2, ('T', 'P'): -1, ('T', 'S'): 1,
            ('T', 'T'): 5, ('T', 'W'): -2, ('T', 'Y'): -2, ('T', 'V'): 0,
            
            ('W', 'A'): -3, ('W', 'R'): -3, ('W', 'N'): -4, ('W', 'D'): -4,
            ('W', 'C'): -2, ('W', 'Q'): -2, ('W', 'E'): -3, ('W', 'G'): -2,
            ('W', 'H'): -2, ('W', 'I'): -3, ('W', 'L'): -2, ('W', 'K'): -3,
            ('W', 'M'): -1, ('W', 'F'): 1, ('W', 'P'): -4, ('W', 'S'): -3,
            ('W', 'T'): -2, ('W', 'W'): 11, ('W', 'Y'): 2, ('W', 'V'): -3,
            
            ('Y', 'A'): -2, ('Y', 'R'): -2, ('Y', 'N'): -2, ('Y', 'D'): -3,
            ('Y', 'C'): -2, ('Y', 'Q'): -1, ('Y', 'E'): -2, ('Y', 'G'): -3,
            ('Y', 'H'): 2, ('Y', 'I'): -1, ('Y', 'L'): -1, ('Y', 'K'): -2,
            ('Y', 'M'): -1, ('Y', 'F'): 3, ('Y', 'P'): -3, ('Y', 'S'): -2,
            ('Y', 'T'): -2, ('Y', 'W'): 2, ('Y', 'Y'): 7, ('Y', 'V'): -1,
            
            ('V', 'A'): 0, ('V', 'R'): -3, ('V', 'N'): -3, ('V', 'D'): -3,
            ('V', 'C'): -1, ('V', 'Q'): -2, ('V', 'E'): -2, ('V', 'G'): -3,
            ('V', 'H'): -3, ('V', 'I'): 3, ('V', 'L'): 1, ('V', 'K'): -2,
            ('V', 'M'): 1, ('V', 'F'): -1, ('V', 'P'): -2, ('V', 'S'): -2,
            ('V', 'T'): 0, ('V', 'W'): -3, ('V', 'Y'): -1, ('V', 'V'): 4,
        }
        
        result = blosum62_full.copy()
        
        # Add default scores for ambiguous codes
        for amb_code, mapped_aas in AMBIGUOUS_MAPPING.items():
            if amb_code in ['-', '*']:
                continue
            if mapped_aas:
                for other_aa in STANDARD_AA:
                    scores = []
                    for mapped_aa in mapped_aas:
                        if (mapped_aa, other_aa) in blosum62_full:
                            scores.append(blosum62_full[(mapped_aa, other_aa)])
                    if scores:
                        avg_score = sum(scores) / len(scores)
                        result[(amb_code, other_aa)] = avg_score
                        result[(other_aa, amb_code)] = avg_score
        
        return result
    
    # ============================================
    # STAGE 1: Load Pairwise Alignment Results
    # ============================================
    
    def stage1_load_alignment_results(self, csv_filepath: str) -> Dict:
        """
        PHASE 1: Load pairwise alignment results from CSV file.
        """
        print(f"\n[STAGE 1 - PHASE 1] Loading alignment results from {csv_filepath}...")
        
        df = pd.read_csv(csv_filepath)
        
        print(f"✓ Loaded {len(df)} pairwise alignments")
        print(f"  Columns: {list(df.columns)}")
        
        # Extract unique protein IDs
        all_protein_ids = set(df['protein_id_1'].unique()) | set(df['protein_id_2'].unique())
        protein_id_list = sorted(list(all_protein_ids))
        
        print(f"  Rows in CSV: {len(df)}")
        print(f"  Unique protein_id_1: {df['protein_id_1'].nunique()}")
        print(f"  Unique protein_id_2: {df['protein_id_2'].nunique()}")
        print(f"  Total unique protein IDs (union): {len(all_protein_ids)}")

        # Create mapping: protein_id -> index
        self.protein_ids = {pid: idx for idx, pid in enumerate(protein_id_list)}
        self.num_sequences = len(protein_id_list)
        
        print(f"✓ Found {self.num_sequences} unique protein sequences")
        
        # Extract species information
        for idx, row in df.iterrows():
            pid1 = row['protein_id_1']
            pid2 = row['protein_id_2']
            
            if pid1 not in self.species_info:
                self.species_info[pid1] = {
                    'species': row['species_1'],
                    'order': row['order_1'],
                    'superorder': row['superorder_1']
                }
            
            if pid2 not in self.species_info:
                self.species_info[pid2] = {
                    'species': row['species_2'],
                    'order': row['order_2'],
                    'superorder': row['superorder_2']
                }
        
        # Create sequence names
        self.sequence_names = []
        for pid in sorted(self.protein_ids.keys(), key=lambda x: self.protein_ids[x]):
            species = self.species_info[pid]['species']
            self.sequence_names.append(species)
        
        print(f"✓ Species: {self.sequence_names[:5]}...")
        
        # Extract alignment scores and alignments
        scores_dict = {}
        alignments_dict = {}
        
        for idx, row in df.iterrows():
            pid1 = row['protein_id_1']
            pid2 = row['protein_id_2']
            
            idx1 = self.protein_ids[pid1]
            idx2 = self.protein_ids[pid2]
            
            score = float(row['score'])
            
            # Store both directions
            scores_dict[(idx1, idx2)] = score
            scores_dict[(idx2, idx1)] = score
            
            # Store alignments
            alignments_dict[(idx1, idx2)] = (row['aln1'], row['aln2'])
            alignments_dict[(idx2, idx1)] = (row['aln2'], row['aln1'])
        
        self.raw_scores_dict = scores_dict
        self.alignments = alignments_dict
        
        # Add self-scores (diagonal) - use maximum score found
        max_score = max(scores_dict.values()) if scores_dict else 100
        self_score = max_score * 1.5
        
        for i in range(self.num_sequences):
            self.raw_scores_dict[(i, i)] = self_score
        
        print(f"✓ Alignment data loaded")
        print(f"  Min score: {min(scores_dict.values()):.2f}")
        print(f"  Max score: {max_score:.2f}")
        print(f"  Self-score: {self_score:.2f}")
        
        return self.raw_scores_dict
    
    # ============================================
    # STAGE 2: Convert Scores to Distance Matrix
    # ============================================
    
    def stage2_scores_to_distances(self) -> np.ndarray:
        """
        PHASE 1 (cont): Convert alignment scores to normalized distance matrix.
        """
        print(f"\n[STAGE 2 - PHASE 1] Converting scores to distance matrix...")
        
        if not self.raw_scores_dict:
            raise ValueError("Scores not loaded. Call stage1_load_alignment_results first.")
        
        n = self.num_sequences
        distances = np.zeros((n, n))
        
        all_scores = list(self.raw_scores_dict.values())
        min_score = min(all_scores)
        max_score = max(all_scores)
        score_range = max_score - min_score
        
        print(f"  Score range: [{min_score:.2f}, {max_score:.2f}]")
        print(f"  Score span: {score_range:.2f}")
        
        for i in range(n):
            for j in range(n):
                if i == j:
                    distances[i][j] = 0
                else:
                    raw_score = self.raw_scores_dict.get((i, j), min_score)
                    
                    if score_range > 0:
                        normalized_score = (raw_score - min_score) / score_range
                    else:
                        normalized_score = 0.5
                    
                    normalized_score = max(0, min(1, normalized_score))
                    distances[i][j] = 1 - normalized_score
        
        self.distance_matrix = distances
        
        upper_triangle = distances[np.triu_indices_from(distances, k=1)]
        print(f"✓ Distance matrix created ({n}×{n})")
        print(f"  Min distance: {upper_triangle.min():.4f}")
        print(f"  Max distance: {upper_triangle.max():.4f}")
        print(f"  Mean distance: {upper_triangle.mean():.4f}")
        print(f"  Std distance: {upper_triangle.std():.4f}")
        
        return self.distance_matrix
    
    # ============================================
    # STAGE 3: Build Guide Tree using UPGMA
    # ============================================
    
    def stage3_build_guide_tree(self) -> Dict:
        """
        PHASE 2: Build UPGMA tree from distance matrix.
        """
        print(f"\n[STAGE 3 - PHASE 2] Building guide tree (UPGMA)...")
        
        if self.distance_matrix is None:
            raise ValueError("Distance matrix not computed. Call stage2_scores_to_distances first.")
        
        clusters = {}
        for i in range(self.num_sequences):
            clusters[i] = {
                'id': i,
                'label': self.sequence_names[i],
                'size': 1
            }
        
        distances = self.distance_matrix.copy()
        node_counter = self.num_sequences
        
        iteration = 0
        while len(clusters) > 1:
            iteration += 1
            
            min_dist = float('inf')
            min_pair = None
            
            cluster_ids = list(clusters.keys())
            
            for i_idx in range(len(cluster_ids)):
                for j_idx in range(i_idx + 1, len(cluster_ids)):
                    i = cluster_ids[i_idx]
                    j = cluster_ids[j_idx]
                    
                    dist = distances[i][j]
                    if dist < min_dist:
                        min_dist = dist
                        min_pair = (i, j)
            
            if min_pair is None:
                break
            
            ci, cj = min_pair
            ni = clusters[ci]['size']
            nj = clusters[cj]['size']
            
            new_id = node_counter
            node_counter += 1
            
            # Store actual node dictionaries
            clusters[new_id] = {
                'id': new_id,
                'label': f"Node_{new_id}",
                'size': ni + nj,
                'left': clusters[ci],
                'right': clusters[cj],
                'left_distance': min_dist / 2,
                'right_distance': min_dist / 2
            }
            
            # Update distance matrix
            max_id = max(clusters.keys()) + 1
            new_distances = np.zeros((max_id, max_id))
            
            for k1 in clusters:
                for k2 in clusters:
                    if k1 == k2:
                        new_distances[k1][k2] = 0
                    elif k1 < k2:
                        if k1 == new_id or k2 == new_id:
                            if k1 == new_id:
                                other = k2
                            else:
                                other = k1
                            
                            if other not in [ci, cj, new_id]:
                                new_dist = (ni * distances[ci][other] + 
                                           nj * distances[cj][other]) / (ni + nj)
                                new_distances[k1][k2] = new_dist
                                new_distances[k2][k1] = new_dist
                        else:
                            new_distances[k1][k2] = distances[k1][k2]
                            new_distances[k2][k1] = distances[k1][k2]
            
            distances = new_distances
            del clusters[ci]
            del clusters[cj]
            
            if iteration % 10 == 0 or len(clusters) <= 2:
                print(f"  Iteration {iteration}: {len(clusters)} clusters remaining, min_dist={min_dist:.6f}")
        
        root_id = list(clusters.keys())[0]
        self.guide_tree = clusters[root_id]
        
        print(f"✓ Guide tree built")
        print(f"  Total iterations: {iteration}")
        
        return self.guide_tree
    
    # ============================================
    # STAGE 4: Progressive MSA with Profiles
    # ============================================
    
    def stage4_progressive_msa(self) -> List[str]:
        """
        PHASE 3&4: Traverse guide tree and perform progressive alignment using profiles.
        """
        print(f"\n[STAGE 4 - PHASE 3&4] Building progressive MSA using guide tree...")
        
        def align_node(node):
            """Recursively align subtrees and return profile"""
            
            if 'left' not in node:
                # Leaf node: single sequence
                seq_id = node['id']
                seq = None
                
                # Try to find alignment for this sequence
                for j in range(self.num_sequences):
                    if j != seq_id:
                        if (seq_id, j) in self.alignments:
                            seq = self.alignments[(seq_id, j)][0]
                            break
                        elif (j, seq_id) in self.alignments:
                            seq = self.alignments[(j, seq_id)][1]
                            break
                
                if not seq:
                    seq = 'M' + '-' * 149
                
                return Profile([seq])
            
            # Internal node: align left and right subtrees
            left_profile = align_node(node['left'])
            right_profile = align_node(node['right'])
            
            # Perform profile-to-profile alignment
            alignment = self._profile_alignment(left_profile, right_profile)
            
            # Create new profile from alignment
            return Profile(alignment)
        
        # Process entire tree
        final_profile = align_node(self.guide_tree)
        self.msa = final_profile.alignment
        self.msa_profiles = final_profile.profile
        
        print(f"✓ Progressive MSA constructed")
        print(f"  MSA sequences: {len(self.msa)}")
        print(f"  Alignment length: {len(self.msa[0]) if self.msa else 0}")
        
        return self.msa
    
    def _profile_alignment(self, profile1: Profile, profile2: Profile) -> List[str]:
        """Profile-to-profile alignment using dynamic programming"""
        m = profile1.length
        n = profile2.length
        gap_penalty = -2
        
        # Initialize DP table
        dp = [[0.0] * (n + 1) for _ in range(m + 1)]
        for i in range(m + 1):
            dp[i][0] = i * gap_penalty
        for j in range(n + 1):
            dp[0][j] = j * gap_penalty
        
        # Fill DP table
        for i in range(1, m + 1):
            for j in range(1, n + 1):
                profile_score = 0.0
                for aa1 in STANDARD_AA:
                    freq1 = profile1.profile[i-1].get(aa1, 0)
                    if freq1 == 0:
                        continue
                    for aa2 in STANDARD_AA:
                        freq2 = profile2.profile[j-1].get(aa2, 0)
                        if freq2 == 0:
                            continue
                        blosum_score = self.blosum62.get((aa1, aa2), -1)
                        profile_score += freq1 * freq2 * blosum_score
                
                match = dp[i-1][j-1] + profile_score
                delete = dp[i-1][j] + gap_penalty
                insert = dp[i][j-1] + gap_penalty
                
                dp[i][j] = max(match, delete, insert)
        
        # Traceback
        alignment1 = []
        alignment2 = []
        
        i, j = m, n
        while i > 0 or j > 0:
            if i == 0:
                alignment1.append('-' * len(profile1.alignment[0]))
                alignment2.append(profile2.alignment[0][j-1] if j <= len(profile2.alignment[0]) else '-')
                j -= 1
            elif j == 0:
                alignment1.append(profile1.alignment[0][i-1] if i <= len(profile1.alignment[0]) else '-')
                alignment2.append('-' * len(profile2.alignment[0]))
                i -= 1
            else:
                profile_score = 0.0
                for aa1 in STANDARD_AA:
                    freq1 = profile1.profile[i-1].get(aa1, 0)
                    if freq1 == 0:
                        continue
                    for aa2 in STANDARD_AA:
                        freq2 = profile2.profile[j-1].get(aa2, 0)
                        if freq2 == 0:
                            continue
                        blosum_score = self.blosum62.get((aa1, aa2), -1)
                        profile_score += freq1 * freq2 * blosum_score
                
                if dp[i][j] == dp[i-1][j-1] + profile_score:
                    alignment1.append(profile1.alignment[0][i-1] if i <= len(profile1.alignment[0]) else '-')
                    alignment2.append(profile2.alignment[0][j-1] if j <= len(profile2.alignment[0]) else '-')
                    i -= 1
                    j -= 1
                elif dp[i][j] == dp[i-1][j] + gap_penalty:
                    alignment1.append(profile1.alignment[0][i-1] if i <= len(profile1.alignment[0]) else '-')
                    alignment2.append('-')
                    i -= 1
                else:
                    alignment1.append('-')
                    alignment2.append(profile2.alignment[0][j-1] if j <= len(profile2.alignment[0]) else '-')
                    j -= 1
        
        alignment1.reverse()
        alignment2.reverse()
        
        result = profile1.alignment + profile2.alignment
        return result
    
    # ============================================
    # STAGE 5: Calculate Distances from MSA
    # ============================================
    
    def stage5_msa_distances(self) -> np.ndarray:
        """Calculate pairwise distances from MSA"""
        print(f"\n[STAGE 5] Calculating distances from MSA...")
        
        if self.msa is None:
            raise ValueError("MSA not built. Call stage4_progressive_msa first.")
        
        n = len(self.msa)
        msa_distances = np.zeros((n, n))
        
        for i in range(n):
            for j in range(n):
                if i == j:
                    msa_distances[i][j] = 0
                else:
                    seq_i = self.msa[i]
                    seq_j = self.msa[j]
                    
                    matches = 0
                    total = 0
                    
                    for k in range(min(len(seq_i), len(seq_j))):
                        if seq_i[k] != '-' and seq_j[k] != '-':
                            total += 1
                            if seq_i[k] == seq_j[k]:
                                matches += 1
                    
                    if total > 0:
                        identity = matches / total
                    else:
                        identity = 0
                    
                    msa_distances[i][j] = 1 - identity
        
        print(f"✓ MSA distances calculated ({n}×{n})")
        
        return msa_distances
    
    # ============================================
    # STAGE 6: Build Final Tree from MSA
    # ============================================
    
    def stage6_final_tree(self, msa_distances: np.ndarray) -> Dict:
        """Build final phylogenetic tree from MSA-based distances"""
        print(f"\n[STAGE 6] Building final phylogenetic tree...")
        
        clusters = {}
        for i in range(self.num_sequences):
            clusters[i] = {
                'id': i,
                'label': self.sequence_names[i],
                'size': 1
            }
        
        distances = {}
        for i in range(self.num_sequences):
            for j in range(self.num_sequences):
                distances[(i, j)] = msa_distances[i][j]
        
        node_counter = self.num_sequences
        
        iteration = 0
        while len(clusters) > 1:
            iteration += 1
            min_dist = float('inf')
            min_pair = None
            
            cluster_ids = list(clusters.keys())
            
            for i_idx in range(len(cluster_ids)):
                for j_idx in range(i_idx + 1, len(cluster_ids)):
                    i = cluster_ids[i_idx]
                    j = cluster_ids[j_idx]
                    
                    dist = distances.get((i, j), distances.get((j, i), float('inf')))
                    if dist < min_dist:
                        min_dist = dist
                        min_pair = (i, j)
            
            if min_pair is None:
                break
            
            ci, cj = min_pair
            ni = clusters[ci]['size']
            nj = clusters[cj]['size']
            
            new_id = node_counter
            node_counter += 1
            
            clusters[new_id] = {
                'id': new_id,
                'label': f"Node_{new_id}",
                'size': ni + nj,
                'left': clusters[ci],
                'right': clusters[cj],
                'left_distance': min_dist / 2,
                'right_distance': min_dist / 2
            }
            
            for k in clusters:
                if k not in [ci, cj, new_id]:
                    dist_ci_k = distances.get((ci, k), distances.get((k, ci), float('inf')))
                    dist_cj_k = distances.get((cj, k), distances.get((k, cj), float('inf')))
                    
                    if dist_ci_k != float('inf') and dist_cj_k != float('inf'):
                        new_dist = (ni * dist_ci_k + nj * dist_cj_k) / (ni + nj)
                        distances[(new_id, k)] = new_dist
                        distances[(k, new_id)] = new_dist
            
            keys_to_remove = [(x, y) for x, y in distances.keys() 
                              if x in [ci, cj] or y in [ci, cj]]
            for key in keys_to_remove:
                del distances[key]
            
            del clusters[ci]
            del clusters[cj]
            
            if iteration % 10 == 0 or len(clusters) <= 2:
                print(f"  Iteration {iteration}: {len(clusters)} clusters remaining")
        
        root_id = list(clusters.keys())[0]
        self.final_tree = clusters[root_id]
        
        print(f"✓ Final phylogenetic tree built")
        print(f"  Total iterations: {iteration}")
        
        return self.final_tree
    
    # ============================================
    # OUTPUT FUNCTIONS
    # ============================================
    
    def tree_to_newick(self, node: Dict = None, include_distances: bool = True) -> str:
        """Convert tree to Newick format by directly traversing tree structure"""
        if node is None:
            node = self.final_tree
        
        def format_node(n):
            if 'left' not in n:
                return n['label']
            
            left_str = format_node(n['left'])
            right_str = format_node(n['right'])
            
            if include_distances:
                left_dist = n.get('left_distance', 0)
                right_dist = n.get('right_distance', 0)
                return f"({left_str}:{left_dist:.6f},{right_str}:{right_dist:.6f})"
            else:
                return f"({left_str},{right_str})"
        
        return format_node(node) + ";"
    
    def save_tree_newick(self, filepath: str):
        """Save final tree in Newick format"""
        newick = self.tree_to_newick()
        with open(filepath, 'w') as f:
            f.write(newick)
        print(f"✓ Tree saved to {filepath}")
        return newick
    
    def save_msa_fasta(self, filepath: str):
        """Save MSA in FASTA format"""
        if self.msa is None:
            print("⚠ MSA not available to save")
            return
        
        with open(filepath, 'w') as f:
            for i, seq in enumerate(self.msa):
                f.write(f">{self.sequence_names[i]}\n{seq}\n")
        print(f"✓ MSA saved to {filepath}")
    
    def save_distances_csv(self, distances: np.ndarray, filepath: str):
        """Save distance matrix to CSV"""
        with open(filepath, 'w', newline='') as f:
            writer = csv.writer(f)
            header = [''] + self.sequence_names
            writer.writerow(header)
            
            for i in range(len(self.sequence_names)):
                row = [self.sequence_names[i]] + list(distances[i])
                writer.writerow(row)
        print(f"✓ Distances saved to {filepath}")
    
    # ============================================
    # MAIN PIPELINE RUNNER
    # ============================================
    
    def run_complete_pipeline(self,
                             alignment_csv: str,
                             output_prefix: str = None) -> Dict:
        """Run complete pipeline from alignment results to final tree"""
        if output_prefix is None:
            output_prefix = self.algorithm_name
        
        print(f"\n{'='*70}")
        print(f"UNIFIED PHYLOGENETIC PIPELINE - {self.algorithm_name}")
        print(f"PHASE 1: Scores → Distances")
        print(f"PHASE 2: Guide Tree (UPGMA)")
        print(f"PHASE 3&4: Progressive MSA with Profiles (26 amino acids)")
        print(f"PHASE 5-6: MSA-based Final Tree")
        print(f"{'='*70}")
        
        # Execute all stages
        self.stage1_load_alignment_results(alignment_csv)
        self.stage2_scores_to_distances()
        self.stage3_build_guide_tree()
        self.stage4_progressive_msa()
        msa_distances = self.stage5_msa_distances()
        self.stage6_final_tree(msa_distances)
        
        # Save outputs
        newick_tree = self.save_tree_newick(f"{output_prefix}_final_tree.nwk")
        self.save_msa_fasta(f"{output_prefix}_msa.fasta")
        self.save_distances_csv(self.distance_matrix, f"{output_prefix}_distance_matrix.csv")
        self.save_distances_csv(msa_distances, f"{output_prefix}_msa_distance_matrix.csv")
        
        print(f"\n{'='*70}")
        print(f"PIPELINE COMPLETE")
        print(f"✓ All {self.num_sequences} sequences processed successfully")
        print(f"{'='*70}\n")
        
        return {
            'algorithm': self.algorithm_name,
            'final_tree': self.final_tree,
            'tree_newick': newick_tree,
            'msa': self.msa,
            'distance_matrix': self.distance_matrix,
            'msa_distances': msa_distances,
            'guide_tree': self.guide_tree,
            'sequence_names': self.sequence_names,
            'num_sequences': self.num_sequences
        }


# ============================================
# USAGE EXAMPLE
# ============================================


if __name__ == "__main__":
    print("Starting Phylogenetic Analysis Pipeline\n")
    
    results = {}
    
    # Run for each algorithm
    for algo in ['NW', 'SW', 'SG']:
        try:
            pipeline = UnifiedPhylogeneticPipeline(algo)
            if('algo'=='NW'):
                po= 'nw_full_results.csv'
            elif('algo'=='SW'):
                po= 'smith_waterman_alignments.csv'
            else:
                po= 'alignment_statistics_detailed.csv'
            result = pipeline.run_complete_pipeline(
                alignment_csv=po,
                output_prefix=f'upgma_results/{algo.lower()}'
            )
            results[algo] = result
        except Exception as e:
            print(f"✗ {algo} pipeline failed: {e}\n")
            import traceback
            traceback.print_exc()
    
    # Print summary
    print("\n" + "="*70)
    print("FINAL SUMMARY")
    print("="*70 + "\n")
    
    for algo, result in results.items():
        print(f"{algo} Algorithm:")
        print(f"  Sequences: {result['num_sequences']}")
        print(f"  Newick length: {len(result['tree_newick'])} characters")
        print()
