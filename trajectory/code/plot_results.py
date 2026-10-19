"""Standalone figure of the controlled coordinate experiment, using saved draws."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
p=ROOT/'results'
old=np.load(p/'mc_published_n32768_step0.005.npz')
new=np.load(p/'mc_corrected_coordinates_n32768_step0.005.npz')
catalog=json.loads((ROOT/'data/key_clusters.json').read_text())['clusters']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(10.5,4.6),sharey=True)
for j,(ax,c) in enumerate(zip(axes,catalog)):
    threshold=c['encounter_threshold_3r50_pc']
    for data,color,label in [(old,'#aa4a30','Printed coordinates'),(new,'#14688a','Corrected coordinates only')]:
        distances=np.sort(data['minimum_pc'][:,j])
        ax.plot(distances,np.arange(1,len(distances)+1)/len(distances)*100,color=color,label=label,lw=2)
    ax.axvline(threshold,c='#333333',ls='--',lw=1.2,label='Encounter threshold: 3 r₅₀' if j==0 else None)
    ax.set_xlim(0,150);ax.set_ylim(0,60)
    ax.set_title(c['name'].replace('_',' '));ax.set_xlabel('Minimum past separation (pc)')
    ax.grid(alpha=.18)
    count=int(old['events'][:,j].sum())
    ax.text(.96,.62,f'Printed: {count}/32,768 encounters\nCorrected: 0/32,768\n95% upper: 0.00914%',transform=ax.transAxes,ha='right',va='bottom',fontsize=10,
            bbox=dict(facecolor='white',alpha=.9,edgecolor='none',pad=3))
axes[0].set_ylabel('Histories with minimum separation below x (%)')
axes[0].legend(loc='upper left',fontsize=8.8)
fig.suptitle('J0737: correcting the sky position removes sampled encounters with both clusters',fontsize=12)
fig.text(.5,.015,'Conditional two-cluster calculation; zero-centered heliocentric Gaussian-mixture radial prior.\nSame distance, proper motions, uncertainties, random draws, Galactic potential and event criterion.',ha='center',fontsize=9)
fig.tight_layout(rect=[0,.10,1,.93])
fig.savefig(p/'coordinate_effect.svg',bbox_inches='tight')
fig.savefig(p/'coordinate_effect.png',dpi=160,bbox_inches='tight')
print(p/'coordinate_effect.svg')
