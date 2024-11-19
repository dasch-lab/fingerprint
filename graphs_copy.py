import ast
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Definitions for regions
chothia_H_definition = [
    {"name": "FRH1", "start": 1, "end": 25},
    {"name": "CDRH1", "start": 26, "end": 32},
    {"name": "FRH2", "start": 33, "end": 52},
    {"name": "CDRH2", "start": 53, "end": 55},
    {"name": "FRH3", "start": 56, "end": 95},
    {"name": "CDRH3", "start": 96, "end": 101},
    {"name": "FRH4", "start": 102, "end": 113},
]

chothia_L_definition = [
    {"name": "FRL1", "start": 1, "end": 25},
    {"name": "CDRL1", "start": 26, "end": 32},
    {"name": "FRL2", "start": 33, "end": 49},
    {"name": "CDRL2", "start": 50, "end": 52},
    {"name": "FRL3", "start": 53, "end": 90},
    {"name": "CDRL3", "start": 91, "end": 96},
    {"name": "FRL4", "start": 97, "end": 109},
]

colors_new = ['#341c3c', '#e43444', '#ad1758', '#701f57', '#f37650', '#f6b48f', '#801848', '#6c1874']

# Function to convert plddtlist into a list of scores
def parse_plddtlist(plddtlist_str):
    plddtlist = ast.literal_eval(plddtlist_str)
    return [score for _, score in plddtlist]

def position_list(plddtlist_str):
    plddtlist = ast.literal_eval(plddtlist_str)
    return [k[1] for k, _ in plddtlist]

def anarci_list(plddtlist_str):
    plddtlist = ast.literal_eval(plddtlist_str)
    return [k[2] for k, _ in plddtlist]

def calculate_max_region_lengths(df, regions_definition):
    region_lengths = {region['name']: [] for region in regions_definition}
    for _, row in df.iterrows():
        anarci_numbers = row['anarci']
        for region in regions_definition:
            region_name = region['name']
            start, end = region['start'], region['end']
            region_count = sum(1 for num in anarci_numbers if start <= num <= end)
            region_lengths[region_name].append(region_count)
    return {region_name: max(lengths) for region_name, lengths in region_lengths.items()}

def fill_matrix(df, matrix, regions_definition, max_amino_acids):
    for i, row in df.iterrows():
        anarci_numbers = row['anarci']
        scores = row['scores']
        index = row['id']
        shift = 0
        for region in regions_definition:
            region_name = region['name']
            max_length = max_amino_acids[region_name]
            region_scores = [(idx, score) for idx, anarci, score in zip(index, anarci_numbers, scores)
                             if region['start'] <= anarci <= region['end']]
            for idx, score in region_scores:
                matrix_col = idx - 1
                if matrix_col < matrix.shape[1] and matrix_col + shift < matrix.shape[1]:
                    matrix[i, matrix_col + shift] = score
            shift += (max_length - len(region_scores))
    return matrix

def add_cdr_boxes(ax, cdrs, heatmap_height):
    for cdr in cdrs:
        start = cdr['start'] - 1
        end = cdr['end']
        rect = Rectangle((start, 0), end - start, heatmap_height, linewidth=2, edgecolor='red', facecolor='none')
        ax.add_patch(rect)

# Remove specified PDB IDs
#pdbids_to_remove = ["8h07", "8g3q", "7lxz", "7sk8", "5w1k"]

# Load VH and VL data from CSV files
df_vh = pd.read_csv('/disk1/fingerprint/vh_plddt_results_more.csv')
df_vl = pd.read_csv('/disk1/fingerprint/vl_plddt_results_more.csv')

# Filter out unwanted PDB IDs
#df_vh = df_vh[~df_vh['pdbid'].isin(pdbids_to_remove)].reset_index(drop=True)
#df_vl = df_vl[~df_vl['pdbid'].isin(pdbids_to_remove)].reset_index(drop=True)

# Parse the plddtlist for VH and VL
df_vh['scores'] = df_vh['plddtlist'].apply(parse_plddtlist)
df_vl['scores'] = df_vl['plddtlist'].apply(parse_plddtlist)
df_vh['id'] = df_vh['plddtlist'].apply(position_list)
df_vl['id'] = df_vl['plddtlist'].apply(position_list)
df_vh['anarci'] = df_vh['plddtlist'].apply(anarci_list)
df_vl['anarci'] = df_vl['plddtlist'].apply(anarci_list)

max_amino_acids_vh = calculate_max_region_lengths(df_vh, chothia_H_definition)
max_amino_acids_vl = calculate_max_region_lengths(df_vl, chothia_L_definition)
max_id_vh = sum(max_amino_acids_vh.values())
max_id_vl = sum(max_amino_acids_vl.values())

vh_matrix = np.zeros((len(df_vh), max_id_vh))
vl_matrix = np.zeros((len(df_vl), max_id_vl))

fill_matrix(df_vh, vh_matrix, chothia_H_definition, max_amino_acids_vh)
fill_matrix(df_vl, vl_matrix, chothia_L_definition, max_amino_acids_vl)

# Adjust CDRs for plotting
chothia_H_cdrs_new = []
current_position_vh = 1
for region_name, max_length in max_amino_acids_vh.items():
    if "CDRH" in region_name:
        chothia_H_cdrs_new.append({
            "name": region_name,
            "start": current_position_vh,
            "end": current_position_vh + max_length - 1
        })
    current_position_vh += max_length

chothia_L_cdrs_new = []
current_position_vl = 1
for region_name, max_length in max_amino_acids_vl.items():
    if "CDRL" in region_name:
        chothia_L_cdrs_new.append({
            "name": region_name,
            "start": current_position_vl,
            "end": current_position_vl + max_length - 1
        })
    current_position_vl += max_length

# Plotting the heatmaps
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
#custom_cmap = sns.color_palette("rocket", as_cmap=True)
custom_cmap = sns.color_palette("mako", as_cmap=True)

sns.heatmap(vh_matrix, ax=ax1, cmap=custom_cmap, cbar_kws={'label': 'pLDDT score'})
ax1.set_title('VH pLDDT Scores')
ax1.set_ylabel('PDB ID')
ax1.set_xlabel('Residue Position')
ax1.set_yticks(np.arange(vh_matrix.shape[0]) + 0.5)
ax1.set_yticklabels(df_vh['pdbid'].values, rotation=0)
add_cdr_boxes(ax1, chothia_H_cdrs_new, vh_matrix.shape[0])

sns.heatmap(vl_matrix, ax=ax2, cmap=custom_cmap, cbar_kws={'label': 'pLDDT score'})
ax2.set_title('VL pLDDT Scores')
ax2.set_ylabel('PDB ID')
ax2.set_xlabel('Residue Position')
ax2.set_yticks(np.arange(vl_matrix.shape[0]) + 0.5)
ax2.set_yticklabels(df_vl['pdbid'].values, rotation=0)
add_cdr_boxes(ax2, chothia_L_cdrs_new, vl_matrix.shape[0])

plt.tight_layout()
plt.savefig('/disk1/fingerprint/img/heatmap_cdrs_tot.png', transparent=True)
plt.show()

import numpy as np
from scipy.stats import ttest_ind


def calculate_statistics(df, regions_definition, max_amino_acids):
    sequence_statistics = {
        'whole_sequence': [],
        'non_cdr': [],
        'cdrs': {region['name']: [] for region in regions_definition if "CDR" in region['name']}
    }

    # Process each sequence row in the dataframe
    for i, row in df.iterrows():
        anarci_numbers = row['anarci']  # ANARCI numbering for the sequence
        scores = row['scores']  # Corresponding pLDDT scores
        index = row['id']  # Indices in the sequence (aligned with ANARCI)

        # Lists to accumulate scores
        whole_seq_scores = []
        non_cdr_scores = []
        cdr_scores = {region['name']: [] for region in regions_definition if "CDR" in region['name']}
        
        for region in regions_definition:
            region_name = region['name']
            max_length = max_amino_acids[region_name]

            # Collect scores for the region
            region_scores = [score for idx, anarci, score in zip(index, anarci_numbers, scores)
                             if region['start'] <= anarci <= region['end']]
            
            # Append scores to the appropriate list
            whole_seq_scores.extend(region_scores)  # For the whole sequence
            if "CDR" in region_name:
                cdr_scores[region_name].extend(region_scores)
            else:
                non_cdr_scores.extend(region_scores)

        # Calculate means for the current sequence
        sequence_statistics['whole_sequence'].append(np.mean(whole_seq_scores))
        sequence_statistics['non_cdr'].append(np.mean(non_cdr_scores))
        
        for cdr_name, scores in cdr_scores.items():
            if scores:  # Only add if there are scores for this region
                sequence_statistics['cdrs'][cdr_name].append(np.mean(scores))

    # Compute aggregate statistics
    overall_mean = np.mean(sequence_statistics['whole_sequence'])
    overall_std = np.std(sequence_statistics['whole_sequence'])

    non_cdr_mean = np.mean(sequence_statistics['non_cdr'])
    non_cdr_std = np.std(sequence_statistics['non_cdr'])

    cdr_stats = {
        cdr_name: {
            'mean': np.mean(scores),
            'std': np.std(scores),
            'min': np.min(scores),
            'max': np.max(scores)
        } for cdr_name, scores in sequence_statistics['cdrs'].items()
    }

    # Printing results
    print(f"Whole Sequence - Mean: {overall_mean:.2f}, Std: {overall_std:.2f}")
    print(f"Non-CDR Regions - Mean: {non_cdr_mean:.2f}, Std: {non_cdr_std:.2f}")
    print("CDR Region Statistics:")
    for cdr, stats in cdr_stats.items():
        print(f"{cdr}: Mean={stats['mean']:.2f}, Std={stats['std']:.2f}, Min={stats['min']:.2f}, Max={stats['max']:.2f}")

    # Return aggregated statistics as a dictionary
    return {
        'overall': {'mean': overall_mean, 'std': overall_std},
        'non_cdr': {'mean': non_cdr_mean, 'std': non_cdr_std},
        'cdr_stats': cdr_stats
    }


def statistics_values(dictionary_plddt, chain):
    # Initialize dictionaries to store mean and std across all sequences for each region
    sequence_mean_whole_seq = {}
    sequence_mean_wtout_cdrs = {}
    
    # Define which CDRs to use based on chain type
    if chain == 'H':
        cdr_regions = ["CDRH1", "CDRH2", "CDRH3"]
    else:
        cdr_regions = ["CDRL1", "CDRL2", "CDRL3"]

    cdr_statistics = {cdr: [] for cdr in cdr_regions}
    
    # Iterate over each sequence (assuming values for all regions are lists of the same length)
    for i in range(len(next(iter(dictionary_plddt.values())))):  
        sequence_whole_seq = []
        sequence_wtout_cdrs = []
        
        # Collect values for each region in the current sequence
        for region, values in dictionary_plddt.items():
            region_values = values[i]  # Values for this sequence in the region
            
            # Whole sequence accumulation
            sequence_whole_seq.extend(region_values)
            
            # Separate CDR and non-CDR (framework) regions
            if region in cdr_regions:
                cdr_statistics[region].extend(region_values)
            else:
                sequence_wtout_cdrs.extend(region_values)
        
        # Store mean for the entire sequence and non-CDR regions
        sequence_mean_whole_seq[i] = np.mean(sequence_whole_seq)
        sequence_mean_wtout_cdrs[i] = np.mean(sequence_wtout_cdrs)
    
    # Calculate overall statistics for whole sequence and non-CDR regions
    overall_mean = np.mean(list(sequence_mean_whole_seq.values()))
    overall_std = np.std(list(sequence_mean_whole_seq.values()))
    overall_min = min(sequence_mean_whole_seq.values())
    overall_max = max(sequence_mean_whole_seq.values())

    non_cdr_mean = np.mean(list(sequence_mean_wtout_cdrs.values()))
    non_cdr_std = np.std(list(sequence_mean_wtout_cdrs.values()))
    non_cdr_min = min(sequence_mean_wtout_cdrs.values())
    non_cdr_max = max(sequence_mean_wtout_cdrs.values())
    
    # Calculate statistics for each CDR region across all sequences
    cdr_stats = {}
    for cdr, values in cdr_statistics.items():
        cdr_stats[cdr] = {
            "mean": np.mean(values),
            "std": np.std(values),
            "min": np.min(values),
            "max": np.max(values)
        }
    
    # Results dictionary
    results = {
        'overall_mean': overall_mean,
        'overall_std': overall_std,
        'overall_min': overall_min,
        'overall_max': overall_max,
        'non_cdr_mean': non_cdr_mean,
        'non_cdr_std': non_cdr_std,
        'non_cdr_min': non_cdr_min,
        'non_cdr_max': non_cdr_max,
        'cdr_stats': cdr_stats
    }
    
    # Print results in a readable format
    print("Overall Sequence Statistics:")
    print(f"Mean: {overall_mean:.2f}, Std: {overall_std:.2f}, Min: {overall_min:.2f}, Max: {overall_max:.2f}")
    
    print("\nFramework (Non-CDR) Region Statistics:")
    print(f"Mean: {non_cdr_mean:.2f}, Std: {non_cdr_std:.2f}, Min: {non_cdr_min:.2f}, Max: {non_cdr_max:.2f}")
    
    print("\nCDR Region Statistics:")
    for cdr, stats in cdr_stats.items():
        print(f"{cdr}: Mean={stats['mean']:.2f}, Std={stats['std']:.2f}, Min={stats['min']:.2f}, Max={stats['max']:.2f}")
    
    return results

import numpy as np
from scipy.stats import ttest_ind
from scipy.stats import wilcoxon

def calculate_statistics_and_compare(df, regions_definition, max_amino_acids):
    sequence_statistics = {
        'whole_sequence': [],
        'non_cdr': [],
        'cdrs': {region['name']: [] for region in regions_definition if "CDR" in region['name']}
    }

    # Process each sequence row in the dataframe
    for _, row in df.iterrows():
        anarci_numbers = row['anarci']  # ANARCI numbering for the sequence
        scores = row['scores']  # Corresponding pLDDT scores
        index = row['id']  # Indices in the sequence (aligned with ANARCI)

        # Lists to accumulate scores
        whole_seq_scores = []
        non_cdr_scores = []
        cdr_scores = {region['name']: [] for region in regions_definition if "CDR" in region['name']}
        
        for region in regions_definition:
            region_name = region['name']
            max_length = max_amino_acids[region_name]

            # Collect scores for the region
            region_scores = [score for idx, anarci, score in zip(index, anarci_numbers, scores)
                             if region['start'] <= anarci <= region['end']]
            
            # Append scores to the appropriate list
            whole_seq_scores.extend(region_scores)  # For the whole sequence
            if "CDR" in region_name:
                cdr_scores[region_name].extend(region_scores)
            else:
                non_cdr_scores.extend(region_scores)

        # Calculate means for the current sequence
        sequence_statistics['whole_sequence'].append(np.mean(whole_seq_scores))
        sequence_statistics['non_cdr'].append(np.mean(non_cdr_scores))
        
        for cdr_name, scores in cdr_scores.items():
            if scores:  # Only add if there are scores for this region
                sequence_statistics['cdrs'][cdr_name].append(np.mean(scores))

    # Compute aggregate statistics
    overall_mean = np.mean(sequence_statistics['whole_sequence'])
    overall_std = np.std(sequence_statistics['whole_sequence'])

    non_cdr_mean = np.mean(sequence_statistics['non_cdr'])
    non_cdr_std = np.std(sequence_statistics['non_cdr'])

    cdr_stats = {
        cdr_name: {
            'mean': np.mean(scores),
            'std': np.std(scores),
            'min': np.min(scores),
            'max': np.max(scores)
        } for cdr_name, scores in sequence_statistics['cdrs'].items()
    }

    # Statistical Tests: Whole Sequence vs Non-CDR
    stat, p_value = wilcoxon( sequence_statistics['whole_sequence'], sequence_statistics['non_cdr'])
    print(p_value)

    # Print statistics
    print(f"\nWhole Sequence - Mean: {overall_mean:.2f}, Std: {overall_std:.2f}")
    print(f"Non-CDR Regions - Mean: {non_cdr_mean:.2f}, Std: {non_cdr_std:.2f}")
    print("\nCDR Region Statistics:")
    for cdr, stats in cdr_stats.items():
        print(f"{cdr}: Mean={stats['mean']:.2f}, Std={stats['std']:.2f}, Min={stats['min']:.2f}, Max={stats['max']:.2f}")
    
    # Print statistical comparison
    print(f"\nStatistical Comparison (Whole Sequence vs Non-CDR):")
    print(f"Wilcoxon Signed-Rank Test: statistic={stat}, p-value={p_value:.3f}")

    # Return statistics and test results
    return {
        'overall': {'mean': overall_mean, 'std': overall_std},
        'non_cdr': {'mean': non_cdr_mean, 'std': non_cdr_std},
        'cdr_stats': cdr_stats,
        'statistical_tests': {
            'whole_vs_non_cdr': {
                'Wilcoxon Signed-Rank Test': stat,
                'p_value': p_value
            }
        }
    }




res_h = calculate_statistics_and_compare(df_vh, chothia_H_definition, max_amino_acids_vh)
res_l = calculate_statistics_and_compare(df_vl, chothia_L_definition, max_amino_acids_vl) 

#res_h = calculate_statistics(df_vh, chothia_H_definition, max_amino_acids_vh)
#res_l = calculate_statistics(df_vl, chothia_L_definition, max_amino_acids_vl)
"""

################# Project statistics

import matplotlib.pyplot as plt
import numpy as np

# Example data: replace these with your actual mean and std values
vh_stats = {
    "Whole Sequence": {"mean": 0.85, "std": 0.04},
    "Non-CDR": {"mean": 0.88, "std": 0.03},
    "CDRH1": {"mean": 0.82, "std": 0.07},
    "CDRH2": {"mean": 0.79, "std": 0.11},
    "CDRH3": {"mean": 0.68, "std": 0.15}
}

vl_stats = {
    "Whole Sequence": {"mean": 0.86, "std": 0.04},
    "Non-CDR": {"mean": 0.87, "std": 0.04},
    "CDRL1": {"mean": 0.80, "std": 0.10},
    "CDRL2": {"mean": 0.85, "std": 0.07},
    "CDRL3": {"mean": 0.79, "std": 0.07}
}

# Function to plot VH or VL statistics with custom colors
def plot_cdr_stats(ax, stats, title, colors):
    labels = list(stats.keys())
    means = [stats[label]["mean"] for label in labels]
    stds = [stats[label]["std"] for label in labels]
    
    # Plotting with individual color and label for each bar to create a legend
    bars = []
    for i, (label, mean, std, color) in enumerate(zip(labels, means, stds, colors)):
        bar = ax.bar(label, mean, yerr=std, capsize=5, color=color, alpha=0.8, edgecolor='black', label=label)
        bars.append(bar)
    
    # Legend and axis settings
    ax.set_title(title)
    ax.set_ylabel("pLDDT Score")
    ax.set_ylim(0, 1)  # Adjust based on pLDDT score range
    ax.set_xlabel("Regions")
    #ax.legend(title="Regions")  # Adding the legend with a title
    ax.grid(axis='y', linestyle='--', alpha=0.7)

# Plotting
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 14), sharey=True)

# VH Plot with custom colors
plot_cdr_stats(ax1, vh_stats, title="VH Sequence pLDDT Scores", colors=colors_new)

# VL Plot with custom colors
plot_cdr_stats(ax2, vl_stats, title="VL Sequence pLDDT Scores", colors=colors_new)

# Overall adjustments
plt.suptitle("Mean pLDDT Scores with Standard Deviation for VH and VL Regions")
plt.tight_layout(rect=[0, 0.03, 1, 0.95])  # Adjust layout for title
plt.show()

plt.savefig('/disk1/fingerprint/img/heatmap_stats.png', transparent=True)#, format='eps', transparent=True)

"""