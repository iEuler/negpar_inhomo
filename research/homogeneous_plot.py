"""Plot homogeneous mixing gains and estimator trajectories."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path)
    args=parser.parse_args()
    data=json.loads(args.summary.read_text())
    if data['status']!='complete':raise ValueError("Only complete experiments can be plotted")
    observables=('anisotropy','fourth_moment','cosine_difference')
    figure,axes=plt.subplots(3,len(data['cases']),figsize=(6*len(data['cases']),9),squeeze=False,constrained_layout=True)
    means,mean_axes=plt.subplots(3,len(data['cases']),figsize=(6*len(data['cases']),9),squeeze=False,constrained_layout=True)
    for column,case in enumerate(data['cases']):
        raw_reference=np.genfromtxt(args.summary.parent/f"epsilon_{case['epsilon']:g}"/'reference/trajectories.csv',delimiter=',',names=True)
        raw_initial=raw_reference[raw_reference['step']==0]
        for j,name in enumerate(observables):
            rows=[r for r in case['records'] if r['observable']==name]
            t=[r['time'] for r in rows]
            axis=axes[j,column]
            axis.plot(t,[r['variance_gain'] for r in rows],':',color='#2474b5',label='Variance gain')
            axis.plot(t,[r['reference_mse_gain'] for r in rows],color='#19814a',label='Reference MSE gain')
            axis.fill_between(t,[r['reference_mse_gain_ci95'][0] for r in rows],
                              [r['reference_mse_gain_ci95'][1] for r in rows],color='#19814a',alpha=.15,label='95% bootstrap interval')
            axis.axhline(1,color='black',linestyle='--',linewidth=1)
            axis.set_title(f"epsilon={case['epsilon']:g}: {name.replace('_',' ')}")
            axis.set_xlabel('Time');axis.set_ylabel('Gain versus better component')
            axis.grid(alpha=.2);axis.legend(fontsize=8)
            axis=mean_axes[j,column]
            for method,color in (('full','#444444'),('signed','#2474b5'),('covariance_mix','#19814a')):
                axis.plot(t,[r['methods'][method]['mean'] for r in rows],color=color,label=method.replace('_',' '))
            reference=np.array([r['reference_mean'] for r in rows]);se=np.array([r['reference_mean_se'] for r in rows])
            # MSE uses exact initial truth; plot raw reference means consistently
            # throughout time so resetting t=0 does not introduce a visual jump.
            reference[0]=raw_initial[f'full_{j}'].mean()
            se[0]=raw_initial[f'full_{j}'].std(ddof=1)/np.sqrt(len(raw_initial))
            axis.plot(t,reference,'--',color='#dd8523',label='Independent fine PIC')
            axis.fill_between(t,reference-2*se,reference+2*se,color='#dd8523',alpha=.15,label='Reference mean ± 2 SE')
            axis.set_title(f"epsilon={case['epsilon']:g}: {name.replace('_',' ')}")
            axis.set_xlabel('Time');axis.set_ylabel('Ensemble mean');axis.grid(alpha=.2);axis.legend(fontsize=8)
    figure.suptitle('Homogeneous Coulomb collisions: mixing before any resampling\nIndependent pilot weights; gains > 1 favor mixing')
    means.suptitle('Homogeneous ensemble means: reference uncertainty and estimator consistency')
    for fig,name in ((figure,'mixing_gain.png'),(means,'observable_means.png')):
        fig.savefig(args.summary.parent/name,dpi=180)
        print(args.summary.parent/name)


if __name__=='__main__':main()
