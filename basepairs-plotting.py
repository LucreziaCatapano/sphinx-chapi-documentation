
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial'],
    'font.size': 12
})

df = pd.read_csv("basepairs_distance_differences.csv")
fig, axes = plt.subplots(2, 1, figsize=(12, 12))

# Boxplot of distance differences
df.boxplot(
    column='Difference_coot', 
    by='Residue_no1', 
    grid=False, 
    patch_artist=True, 
    boxprops=dict(facecolor='#1f77b4', color='black', alpha=0.8),
    medianprops=dict(color='black', linewidth=1.5),
    whiskerprops=dict(color='black', linewidth=1.5),
    capprops=dict(color='black', linewidth=1.5),
    flierprops=dict(marker='o', color='red', alpha=0.6),
    ax=axes[0]
    
)


axes[0].set_xlabel('Residue Number', fontsize=14, fontweight='bold')
axes[0].set_ylabel('Distance Difference', fontsize=14, fontweight='bold')
axes[0].set_xticklabels(df['Residue_no1'].unique(), rotation=45)
axes[0].set_title('Basepairs Distance Differences After Refinement', fontsize=16, fontweight='bold')
axes[0].grid(axis='y', linestyle='--', alpha=0.5)

# Calculate the average difference for each residue
average_differences = df.groupby('Residue_no1')['Difference_coot'].mean()

# Bar plot
average_differences.plot(kind='bar', color='#ff7f0e', alpha=0.8, ax=axes[1])
axes[1].set_xlabel('Residue Number', fontsize=14, fontweight='bold')
axes[1].set_ylabel('Average Distance Difference', fontsize=14, fontweight='bold')
axes[1].set_title('Average Distance Differences per Residue', fontsize=16, fontweight='bold')
axes[1].set_xticklabels(average_differences.index, rotation=45)
axes[1].grid(axis='y', linestyle='--', alpha=0.5)

# Adjust layout
plt.suptitle("")
plt.tight_layout()  # Adjust the top to make room for the title

# Save the figure
fig.savefig("basepairs_distance_differences.png", dpi=300, bbox_inches='tight')
plt.show()
