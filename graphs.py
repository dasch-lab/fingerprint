import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Data from the table
features = [
    'geom, chem, flex', 'geom, chem', 'geom, flex', 
    'chem, flex', 'geom', 'chem', 'flex'
]

# Mean values for each feature combination for each column
means = {
    'PPI flexibility': [0.783, 0.597, 0.581, 0.742, 0.538, 0.582, 0.529],
    'PPI iterative flexibility': [0.778, 0.613, 0.543, 0.751, 0.555, 0.592, 0.523],
    'Ab-Ag inference flexibility': [0.765, 0.595, 0.626, 0.738, 0.502, 0.579, 0.619],
    'Ab-Ag inference iterative flexibility': [0.765, 0.598, 0.604, 0.747, 0.524, 0.585, 0.592],
    'Ab-Ag fine flexibility': [0.863, 0.669, 0.655, 0.810, 0.567, 0.587, 0.619],
    'Ab-Ag fine iterative flexibility': [0.866, 0.562, 0.599, 0.832, 0.602, 0.501, 0.551],
    'Ab-Ag flexibility': [0.895, 0.565, 0.654, 0.886, 0.541, 0.523, 0.651],
    'Ab-Ag iterative flexibility': [0.910, 0.629, 0.794, 0.905, 0.473, 0.640, 0.800]
}

# Standard deviations for each feature combination for each column
std_devs = {
    'PPI flexibility': [0.003, 0.006, 0.026, 0.001, 0.025, 0.034, 0.018],
    'PPI iterative flexibility': [0.008, 0.005, 0.024, 0.008, 0.028, 0.009, 0.025],
    'Ab-Ag inference flexibility': [0.017, 0.009, 0.011, 0.008, 0.023, 0.031, 0.011],
    'Ab-Ag inference iterative flexibility': [0.011, 0.014, 0.052, 0.002, 0.016, 0.036, 0.078],
    'Ab-Ag fine flexibility': [0.006, 0.106, 0.007, 0.024, 0.018, 0.184, 0.005],
    'Ab-Ag fine iterative flexibility': [0.002, 0.104, 0.022, 0.002, 0.018, 0.086, 0.030],
    'Ab-Ag flexibility': [0.002, 0.110, 0.040, 0.001, 0.028, 0.068, 0.052],
    'Ab-Ag iterative flexibility': [0.002, 0.082, 0.037, 0.006, 0.047, 0.072, 0.051]
}

#colors = ['#9999ce', '#867375', '#eb7537', '#0f9dd5', '#c40505', '#feda15', '#666666', '#0053d7']
colors = ['#341c3c', '#e43444', '#ad1758', '#701f57', '#f37650', '#f6b48f', '#801848', '#6c1874']
# Create a figure with 8 subplots
fig, axs = plt.subplots(2, 4, figsize=(30, 15))  # 4 rows and 2 columns
fig.suptitle('Ablation Study for Different Feature Combinations', fontsize=20)

# Plotting the bar plots with error bars
for i, (key, ax) in enumerate(zip(means.keys(), axs.ravel())):
    ax.bar(features, means[key], yerr=std_devs[key], capsize=5, color = colors[i])
    ax.set_ylim(0, 1)  # Set y-axis limit to a maximum of 1
    ax.set_title(key, fontsize=16)
    ax.set_ylabel('Performance Magnitude')
    ax.set_xlabel('Feature Combination')
    ax.set_xticklabels(features, rotation=45, ha='right')

# Adjust layout for better spacing
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
plt.savefig('/disk1/fingerprint/img/ablation_study_plots.png', transparent=True)
#plt.savefig('/disk1/fingerprint/img/ablation_study_plots.eps', format='eps', transparent=True)

#############
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import ast

# Load the CSV file
df = pd.read_csv('/disk1/fingerprint/protein_analysis_results_residual_5layer.csv')

# Convert the 'plddt' column (string of lists) to actual lists
df['plddt'] = df['plddt'].apply(ast.literal_eval)

# Extract the first and second elements of the plddt list into separate columns (Antigen and Antibody)
df['plddt_x'] = df['plddt'].apply(lambda x: x[0])  # Antigen
df['plddt_y'] = df['plddt'].apply(lambda x: x[1])  # Antibody

# Define PLDDT thresholds for both Antigen and Antibody
def classify_plddt(plddt_value):
    if plddt_value > 90:
        return 'Very High'
    elif 70 < plddt_value <= 90:
        return 'Confident'
    elif 50 < plddt_value <= 70:
        return 'Low'
    else:
        return 'Very Low'

# Classify both Antigen (plddt_x) and Antibody (plddt_y)
df['plddt_x_category'] = df['plddt_x'].apply(classify_plddt)
df['plddt_y_category'] = df['plddt_y'].apply(classify_plddt)

# Calculate confusion matrix elements (TP, TN, FP, FN)
def calculate_confusion_elements(row):
    if row['gt'] == 1 and row['probability'] > 0.5:  # True Positive
        return 'TP'
    elif row['gt'] == 0 and row['probability'] <= 0.5:  # True Negative
        return 'TN'
    elif row['gt'] == 0 and row['probability'] > 0.5:  # False Positive
        return 'FP'
    elif row['gt'] == 1 and row['probability'] <= 0.5:  # False Negative
        return 'FN'

# Apply confusion calculation to the dataframe
df['confusion'] = df.apply(calculate_confusion_elements, axis=1)

# Function to compute normalized confusion matrix elements for a given PLDDT category
def compute_confusion_matrix(df, plddt_category_column):
    categories = ['Very High', 'Confident', 'Low', 'Very Low']
    confusion_counts = {category: {'TP': 0, 'TN': 0, 'FP': 0, 'FN': 0} for category in categories}
    total_counts = {category: 0 for category in categories}

    # Count TP, TN, FP, FN and total samples for each PLDDT category
    for category in categories:
        category_data = df[df[plddt_category_column] == category]
        total_counts[category] = len(category_data)
        confusion_summary = category_data['confusion'].value_counts()
        confusion_counts[category]['TP'] = confusion_summary.get('TP', 0)
        confusion_counts[category]['TN'] = confusion_summary.get('TN', 0)
        confusion_counts[category]['FP'] = confusion_summary.get('FP', 0)
        confusion_counts[category]['FN'] = confusion_summary.get('FN', 0)

    # Normalize the confusion counts
    for category in confusion_counts.keys():
        for key in confusion_counts[category]:
            if total_counts[category] > 0:
                confusion_counts[category][key] /= total_counts[category]  # Normalize by total count

    # Create a DataFrame for easy plotting
    confusion_df = pd.DataFrame(confusion_counts).T
    return confusion_df, total_counts

# Compute confusion matrices for Antigen and Antibody
confusion_df_x, total_counts_x = compute_confusion_matrix(df, 'plddt_x_category')
confusion_df_y, total_counts_y = compute_confusion_matrix(df, 'plddt_y_category')
colors_new = ['#341c3c', '#e43444', '#ad1758', '#701f57', '#f37650', '#f6b48f', '#801848', '#6c1874']

# Plotting
fig, axes = plt.subplots(1, 2, figsize=(15, 6), sharey=True)

# Plot for Antigen (plddt_x)
confusion_df_x.plot(kind='bar', ax=axes[0], color=colors_new[:len(confusion_df_x)])
axes[0].set_title('Antigen')
axes[0].set_xlabel('PLDDT Category')
axes[0].set_ylabel('Normalized Count')
axes[0].set_ylim(0, 1)
axes[0].grid(axis='y')
axes[0].set_xticklabels(axes[0].get_xticklabels(), rotation=0)

# Plot for Antibody (plddt_y)
confusion_df_y.plot(kind='bar', ax=axes[1], color=colors_new[:len(confusion_df_y)])
axes[1].set_title('Antibody')
axes[1].set_xlabel('PLDDT Category')
axes[1].set_ylim(0, 1)
axes[1].grid(axis='y')
axes[1].set_xticklabels(axes[1].get_xticklabels(), rotation=0)

# Adjust layout and show the plots
plt.tight_layout()
plt.show()

# Print the number of points for each category for both Antigen and Antibody
print("Number of points for each pLDDT category (Antigen):")
for category in total_counts_x:
    print(f"{category}: {total_counts_x[category]}")

print("\nNumber of points for each pLDDT category (Antibody):")
for category in total_counts_y:
    print(f"{category}: {total_counts_y[category]}")
    
# Save the plots as a file
plt.savefig('/disk1/fingerprint/img/antigen_antibody_correlation_normalized.png', transparent=True)


plt.savefig('/disk1/fingerprint/img/correlation_normalized_pypx.eps', format='eps', transparent=True)

############### CDRS
import ast
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap

# Function to convert plddtlist into a list of scores
def parse_plddtlist(plddtlist_str):
    plddtlist = ast.literal_eval(plddtlist_str)  # Safely evaluate the string
    scores = [score for _, score in plddtlist]  # Extract pLDDT scores
    return scores

def position_list(plddtlist_str):
    plddtlist = ast.literal_eval(plddtlist_str)  # Safely evaluate the string
    position = [k[1] for k, _ in plddtlist]   # Extract pLDDT scores
    return position

# Load VH and VL data from CSV files
df_vh = pd.read_csv('/disk1/fingerprint/vh_plddt_results_more.csv')
df_vh = df_vh.head(80)
df_vl = pd.read_csv('/disk1/fingerprint/vl_plddt_results_more.csv')
df_vl = df_vl.head(80)

# Parse the plddtlist for VH and VL
df_vh['scores'] = df_vh['plddtlist'].apply(parse_plddtlist)
df_vl['scores'] = df_vl['plddtlist'].apply(parse_plddtlist)
df_vh['id'] = df_vh['plddtlist'].apply(position_list)
df_vl['id'] = df_vl['plddtlist'].apply(position_list)

# Determine the maximum sequence length for heatmap consistency
max_id_vh = df_vh['id'].apply(max).max()  # Get the maximum ID for VH
max_id_vl = df_vl['id'].apply(max).max()  # Get the maximum ID for VL

vh_matrix = np.zeros((len(df_vh), max_id_vh))
vl_matrix = np.zeros((len(df_vl), max_id_vl))

for i, row in df_vh.iterrows():
    for (position, score) in zip(row['id'],row['scores']):
        vh_matrix[i, position-1] = score  # Assign score to the corresponding position

for i, row in df_vl.iterrows():
    for (position, score) in zip(row['id'],row['scores']):
        vl_matrix[i, position-1] = score  # Assign score to the corresponding position


# Pad scores with NaN to ensure each row has the same length
#df_vh['scores_padded'] = df_vh.apply(lambda row: pad_scores_with_positions(row['scores'], row['id'], max_len_vh), axis=1)
#df_vl['scores_padded'] = df_vl.apply(lambda row: pad_scores_with_positions(row['scores'], row['id'], max_len_vl), axis=1)

# Create a heatmap matrix for VH and VL
#vh_matrix = np.array(df_vh['scores_padded'].tolist())
#vl_matrix = np.array(df_vl['scores_padded'].tolist())

# Define CDR regions (Chothia numbering)
chothia_H_cdrs = [
    {"name": "CDRH1", "start": 26, "end": 32},
    {"name": "CDRH2", "start": 53, "end": 55},
    {"name": "CDRH3", "start": 96, "end": 101},
]

chothia_L_cdrs = [
    {"name": "CDRL1", "start": 26, "end": 32},
    {"name": "CDRL2", "start": 50, "end": 52},
    {"name": "CDRL3", "start": 91, "end": 96},
]

# Function to add rectangles to highlight CDRs
def add_cdr_boxes(ax, cdrs, heatmap_height):
    for cdr in cdrs:
        start = cdr['start'] - 1  # Convert 1-based to 0-based index
        end = cdr['end']  # End is inclusive, so no need to subtract 1
        # Add rectangle for the entire CDR region
        rect = Rectangle((start, 0), end-start, heatmap_height, linewidth=2, edgecolor='red', facecolor='none')
        ax.add_patch(rect)

# Plotting the heatmaps
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
custom_cmap = sns.color_palette("mako", as_cmap=True)
tick_positions = np.arange(0, vh_matrix.shape[1], 5) + 0.5  # +0.5 for centering
tick_labels = np.arange(1, vh_matrix.shape[1] + 1, 5)  #
# Plot heatmap for VH
sns.heatmap(vh_matrix, ax=ax1, cmap=custom_cmap, cbar_kws={'label': 'pLDDT score'})
ax1.set_title('VH pLDDT Scores')
ax1.set_ylabel('PDB ID')
ax1.set_xlabel('Residue Position')
ax1.set_yticks(np.arange(vh_matrix.shape[0]) + 0.5)  # Set ticks in the middle of each row
ax1.set_yticklabels(df_vh['pdbid'].values, rotation=0)
ax1.set_xticks(tick_positions)  # Set ticks at the defined positions
ax1.set_xticklabels(tick_labels, rotation=90)  # Labels starting from 1

# Highlight CDRs in VH
add_cdr_boxes(ax1, chothia_H_cdrs, vh_matrix.shape[0])
tick_positions = np.arange(0, vl_matrix.shape[1], 5) + 0.5  # +0.5 for centering
tick_labels = np.arange(1, vl_matrix.shape[1] + 1, 5)  #
# Plot heatmap for VL
sns.heatmap(vl_matrix, ax=ax2, cmap=custom_cmap, cbar_kws={'label': 'pLDDT score'})
ax2.set_title('VL pLDDT Scores')
ax2.set_ylabel('PDB ID')
ax2.set_xlabel('Residue Position')
#ax2.set_yticklabels(df_vl['pdbid'])
ax2.set_yticks(np.arange(vl_matrix.shape[0]) + 0.5)  # Set ticks in the middle of each row
ax2.set_yticklabels(df_vl['pdbid'].values, rotation=0)
ax2.set_xticks(tick_positions)  # Set ticks at the defined positions
ax2.set_xticklabels(tick_labels, rotation=90)  # Labels starting from 1

add_cdr_boxes(ax2, chothia_L_cdrs, vl_matrix.shape[0])

plt.tight_layout()
plt.show()
plt.savefig('/disk1/fingerprint/img/heatmap_cdrs.png')#, format='eps', transparent=True)

# Function to calculate mean and std for the entire sequence, without CDRs, and for CDRs
# Function to calculate statistics for each sequence
def calculate_sequence_statistics(scores, cdrs):
    # Calculate overall mean and std for the entire sequence
    overall_mean = np.mean(scores)
    overall_std = np.std(scores)

    # Calculate the indices of CDRs
    cdr_indices = []
    for cdr in cdrs:
        start = cdr['start'] - 1  # Convert to 0-based index
        end = cdr['end']  # End is inclusive
        cdr_indices.extend(range(start, end))

    # Get the non-CDR scores
    non_cdr_scores = np.delete(scores, cdr_indices)
    non_cdr_mean = np.mean(non_cdr_scores)
    non_cdr_std = np.std(non_cdr_scores)

    # Initialize a dictionary to hold CDR statistics
    cdr_stats = {}
    
    # Calculate mean and std for each CDR
    for cdr in cdrs:
        start = cdr['start'] - 1  # Convert to 0-based index
        end = cdr['end']  # End is inclusive
        cdr_scores = scores[start:end]
        cdr_mean = np.mean(cdr_scores)
        cdr_std = np.std(cdr_scores)
        cdr_stats[cdr['name']] = {'mean': cdr_mean, 'std': cdr_std}

    return {
        'overall_mean': overall_mean,
        'overall_std': overall_std,
        'non_cdr_mean': non_cdr_mean,
        'non_cdr_std': non_cdr_std,
        'cdr_stats': cdr_stats
    }

# Calculate statistics for each VH sequence
vh_statistics_list = []
for index, row in df_vh.iterrows():
    sequence_stats = calculate_sequence_statistics(row['scores'], chothia_H_cdrs)
    vh_statistics_list.append({
        'pdbid': row['pdbid'],
        **sequence_stats
    })

# Calculate overall statistics for VH
overall_vh_statistics = {
    'overall_mean': np.mean([stats['overall_mean'] for stats in vh_statistics_list]),
    'overall_std': np.std([stats['overall_std'] for stats in vh_statistics_list]),
    'non_cdr_mean': np.mean([stats['non_cdr_mean'] for stats in vh_statistics_list]),
    'non_cdr_std': np.std([stats['non_cdr_std'] for stats in vh_statistics_list]),
    'cdr_stats': {cdr['name']: {'mean': np.mean([stats['cdr_stats'][cdr['name']]['mean'] for stats in vh_statistics_list]),
                                  'std': np.std([stats['cdr_stats'][cdr['name']]['std'] for stats in vh_statistics_list])}
                                 for cdr in chothia_H_cdrs}
}

"""# Print results for each VH sequence
print("VH Statistics per sequence:")
for stats in vh_statistics_list:
    print(f"PDB ID: {stats['pdbid']} - Overall Mean: {stats['overall_mean']:.2f}, Overall Std: {stats['overall_std']:.2f}, "
          f"Non-CDR Mean: {stats['non_cdr_mean']:.2f}, Non-CDR Std: {stats['non_cdr_std']:.2f}")
    for cdr_name, cdr_stat in stats['cdr_stats'].items():
        print(f"  {cdr_name} - Mean: {cdr_stat['mean']:.2f}, Std: {cdr_stat['std']:.2f}")
"""
# Print overall statistics for VH
print("\nOverall VH Statistics:")
print(f"Overall Mean: {overall_vh_statistics['overall_mean']:.2f}, Overall Std: {overall_vh_statistics['overall_std']:.2f}")
print(f"Non-CDR Mean: {overall_vh_statistics['non_cdr_mean']:.2f}, Non-CDR Std: {overall_vh_statistics['non_cdr_std']:.2f}")
for cdr_name, stats in overall_vh_statistics['cdr_stats'].items():
    print(f"{cdr_name} - Mean: {stats['mean']:.2f}, Std: {stats['std']:.2f}")

# Repeat similar calculations for VL
vl_statistics_list = []
for index, row in df_vl.iterrows():
    sequence_stats = calculate_sequence_statistics(row['scores'], chothia_L_cdrs)
    vl_statistics_list.append({
        'pdbid': row['pdbid'],
        **sequence_stats
    })

# Calculate overall statistics for VL
overall_vl_statistics = {
    'overall_mean': np.mean([stats['overall_mean'] for stats in vl_statistics_list]),
    'overall_std': np.std([stats['overall_std'] for stats in vl_statistics_list]),
    'non_cdr_mean': np.mean([stats['non_cdr_mean'] for stats in vl_statistics_list]),
    'non_cdr_std': np.std([stats['non_cdr_std'] for stats in vl_statistics_list]),
    'cdr_stats': {cdr['name']: {'mean': np.mean([stats['cdr_stats'][cdr['name']]['mean'] for stats in vl_statistics_list]),
                                  'std': np.std([stats['cdr_stats'][cdr['name']]['std'] for stats in vl_statistics_list])}
                                 for cdr in chothia_L_cdrs}
}

"""# Print results for each VL sequence
print("\nVL Statistics per sequence:")
for stats in vl_statistics_list:
    print(f"PDB ID: {stats['pdbid']} - Overall Mean: {stats['overall_mean']:.2f}, Overall Std: {stats['overall_std']:.2f}, "
          f"Non-CDR Mean: {stats['non_cdr_mean']:.2f}, Non-CDR Std: {stats['non_cdr_std']:.2f}")
    for cdr_name, cdr_stat in stats['cdr_stats'].items():
        print(f"  {cdr_name} - Mean: {cdr_stat['mean']:.2f}, Std: {cdr_stat['std']:.2f}")"""

# Print overall statistics for VL
print("\nOverall VL Statistics:")
print(f"Overall Mean: {overall_vl_statistics['overall_mean']:.2f}, Overall Std: {overall_vl_statistics['overall_std']:.2f}")
print(f"Non-CDR Mean: {overall_vl_statistics['non_cdr_mean']:.2f}, Non-CDR Std: {overall_vl_statistics['non_cdr_std']:.2f}")
for cdr_name, stats in overall_vl_statistics['cdr_stats'].items():
    print(f"{cdr_name} - Mean: {stats['mean']:.2f}, Std: {stats['std']:.2f}")

# Continue with your plotting code...
########## Threshold for Binary
import matplotlib.pyplot as plt
import numpy as np

# Data
labels = ['PPI flex', 'PPI rec', 'PPI flex (Ab)', 'PPI rec (Ab)', 'Ab flex fine', 'Ab rec fine', 'Ab flex', 'Ab rec']
x_values = [60, 62.5, 65, 67.5, 70, 72.5, 75, 77.5, 80]

# Residue values
residue_data = {
    'PPI flex': [0.746, 0.740, 0.734, 0.728, 0.721, 0.715, 0.706, 0.697, 0.684],
    'PPI rec': [0.725, 0.718, 0.711, 0.704, 0.696, 0.688, 0.679, 0.669, 0.661],
    'PPI flex (Ab)': [0.64, 0.641, 0.635, 0.629, 0.640, 0.627, 0.617, 0.606, 0.595],
    'PPI rec (Ab)': [0.672, 0.667, 0.660, 0.651, 0.624, 0.618, 0.613, 0.612, 0.618],
    'Ab flex fine': [0.789, 0.789, 0.788, 0.788, 0.788, 0.787, 0.785, 0.783, 0.774],
    'Ab rec fine': [0.787, 0.784, 0.781, 0.780, 0.776, 0.774, 0.766, 0.759, 0.744],
    'Ab flex': [0.849, 0.849, 0.850, 0.851, 0.851, 0.850, 0.845, 0.841, 0.831],
    'Ab rec': [0.844, 0.843, 0.842, 0.838, 0.835, 0.832, 0.824, 0.814, 0.798]
}

# Atomic values
atomic_data = {
    'PPI flex': [0.761, 0.758, 0.755, 0.752, 0.750, 0.747, 0.744, 0.739, 0.735],
    'PPI rec': [0.735, 0.732, 0.725, 0.719, 0.714, 0.706, 0.696, 0.685, 0.673],
    'PPI flex (Ab)': [0.718, 0.715, 0.716, 0.716, 0.715, 0.713, 0.712, 0.714, 0.717],
    'PPI rec (Ab)': [0.733, 0.734, 0.734, 0.735, 0.734, 0.734, 0.732, 0.727, 0.719],
    'Ab flex fine': [0.851, 0.850, 0.848, 0.847, 0.846, 0.841, 0.838, 0.834, 0.827],
    'Ab rec fine': [0.788, 0.786, 0.786, 0.786, 0.784, 0.782, 0.776, 0.769, 0.762],
    'Ab flex': [0.828, 0.827, 0.826, 0.823, 0.818, 0.811, 0.798, 0.777, 0.742],
    'Ab rec': [0.854, 0.854, 0.853, 0.850, 0.847, 0.842, 0.831, 0.815, 0.795]
}

# Create subplots
fig, axs = plt.subplots(2, 1, figsize=(10, 8))

# Plot residue data
for label in labels:
    axs[0].plot(x_values, residue_data[label], marker='o', label=label)
axs[0].set_title('Residue pLDDT')
axs[0].set_xlabel('Threshold')
axs[0].set_ylabel('Residue')
axs[0].legend()

# Plot atomic data
for label in labels:
    axs[1].plot(x_values, atomic_data[label], marker='o', label=label)
axs[1].set_title('Atomic pLDDT')
axs[1].set_xlabel('Threshold')
axs[1].set_ylabel('Atomic')
axs[1].legend()

# Adjust layout
plt.tight_layout()
plt.show()
plt.savefig('/disk1/fingerprint/img/threshold.png')

