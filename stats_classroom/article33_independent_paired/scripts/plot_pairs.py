from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('data/ecg_rates.csv')
sns.set_theme(style='whitegrid', font_scale=1.0)
fig, ax = plt.subplots(figsize=(6.2, 6.3), dpi=180)
for _, row in df.iterrows():
    ax.plot(['Algorithm A', 'Algorithm B'], [row.A_rate, row.B_rate], color='#9bb8d3', lw=1.6, zorder=1)
ax.scatter(['Algorithm A'] * len(df), df.A_rate, s=70, color='#155b87', label='A', zorder=3)
ax.scatter(['Algorithm B'] * len(df), df.B_rate, s=70, color='#e07a5f', label='B', zorder=3)
for _, row in df.iterrows():
    ax.text(-0.055, row.A_rate, f"ID {int(row.id)}", fontsize=10, color='#155b87',ha='right',va='center')
ax.set_title('Paired heart-rate estimates', weight='bold',fontsize=14,pad=14)
ax.set_ylabel('Heart rate (beats/min)',fontsize=12)
ax.set_xlabel('')
ax.set_ylim(54, 99)
ax.set_xlim(-0.26, 1.12)
fig.tight_layout()
fig.savefig(Path('results/paired_plot.png'), dpi=300, bbox_inches='tight')
