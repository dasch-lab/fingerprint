import matplotlib.pyplot as plt
import numpy as np

# Data
models = ['original', 'flexibility (1.0)', 'modified flexibility', 'recurrent flexibility', 'weighted flexibility']
roc_auc_ppi = [0.860, 0.863, 0.861, 0.871, 0.857]
roc_auc_ab_ag_inference = [0.834, 0.830, 0.594, 0.607, 0.406]
roc_auc_ab_ag_fine_tuning = [0.893, 0.896, 0.903, 0.896, 0.899]

x = np.arange(len(models))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(figsize=(12, 8))

# Bar charts
bars1 = ax.bar(x - width, roc_auc_ppi, width, label='PPI')
bars2 = ax.bar(x, roc_auc_ab_ag_inference, width, label='Ab-Ag Inference')
bars3 = ax.bar(x + width, roc_auc_ab_ag_fine_tuning, width, label='Ab-Ag Fine-Tuning')

# Line charts
#ax.plot(x, roc_auc_ppi, marker='o', color='blue', linewidth=2)
#ax.plot(x, roc_auc_ab_ag_inference, marker='o', color='orange', linewidth=2)
#ax.plot(x, roc_auc_ab_ag_fine_tuning, marker='o', color='green', linewidth=2)

ax.set_xlabel('Models', fontsize=14)
ax.set_ylabel('ROC-AUC Scores', fontsize=14)
ax.set_title('ROC-AUC Performance Comparison Across Models', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(models, fontsize=12)

# Add some text for labels, title and custom x-axis tick labels, etc.
#ax.set_xlabel('Models')
#ax.set_ylabel('ROC-AUC Scores')
#ax.set_title('ROC-AUC Performance Comparison Across Models')
#ax.set_xticks(x)
#ax.set_xticklabels(models)
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))

# Annotate the bars with the ROC-AUC scores
def annotate_bars(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=10)

annotate_bars(bars1)
annotate_bars(bars2)
annotate_bars(bars3)

plt.tight_layout()
plt.show()
plt.savefig('prova.png')


import matplotlib.pyplot as plt
import numpy as np

# Data
features = [
    'geom + chem + flex',
    'geom + chem',
    'geom + flex',
    'chem + flex',
    'geom',
    'chem',
    'flex'
]
roc_auc_flex_1_0 = [0.863, 0.861, 0.747, 0.858, 0.748, 0.846, 0.708]
roc_auc_recurrent_flex = [0.871, 0.869, 0.751, 0.862, 0.752, 0.868, 0.752]
roc_auc_modified_flex = [0.903, 0.901, 0.848, 0.891, 0.821, 0.887, 0.806]

x = np.arange(len(features))  # the label locations
width = 0.25  # the width of the bars

fig, ax = plt.subplots(figsize=(20, 8))

# Bar charts
bars1 = ax.bar(x - width, roc_auc_flex_1_0, width, label='Flexibility 1.0 (PPI)', color='steelblue')
bars2 = ax.bar(x, roc_auc_recurrent_flex, width, label='Recurrent Flexibility (PPI)', color='lightcoral')
bars3 = ax.bar(x + width, roc_auc_modified_flex, width, label='Modified Flexibility (Ab-Ag)', color='seagreen')
plt.ylim(0, 1)

# Labels, title, and legend
ax.set_xlabel('Features', fontsize=14)
ax.set_ylabel('ROC-AUC Scores', fontsize=14)
ax.set_title('ROC-AUC Scores for Different Feature Combinations and Flexibility Types', fontsize=16)
ax.set_xticks(x)
ax.set_xticklabels(features, rotation=45, ha='right', fontsize=12)
ax.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize=12)

# Annotate the bars with the ROC-AUC scores
def annotate_bars(bars):
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.3f}', xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3), textcoords="offset points",
                    ha='center', va='bottom', fontsize=10)

annotate_bars(bars1)
annotate_bars(bars2)
annotate_bars(bars3)

# Adjust layout to make room for the legend
plt.tight_layout(rect=[0, 0, 0.8, 1])  # Adjust the right side of the figure

plt.show()

plt.savefig('ablation.png')
