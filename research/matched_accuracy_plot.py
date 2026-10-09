"""Plot measured homogeneous cost/error frontiers and fresh validation ratios."""
import argparse
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary',type=Path)
    args=parser.parse_args()
    data=json.loads(args.summary.read_text())
    if data['status']!='complete':raise ValueError('Only complete benchmark results may be plotted')
    cases=data['cases']
    figure,axes=plt.subplots(1,len(cases),figsize=(5*len(cases),4.7),squeeze=False,sharey=True,constrained_layout=True)
    styles={'full':('PIC','#444444','o'),'pic_cv':('PIC + initial-moment CV','#9467bd','s'),
            'mixed':('HDP + mixing','#2474b5','o'),'mixed_cv':('HDP + mixing + CV','#19814a','s')}
    for axis,case in zip(axes[0],cases):
        for method,(label,color,marker) in styles.items():
            rows=[r['metrics'][method] for r in case['candidates'] if method in r['metrics']]
            axis.scatter([r['joint_relative_rmse'] for r in rows],[r['mean_compute_seconds'] for r in rows],
                         color=color,marker=marker,alpha=.3,s=18)
            frontier=[r for r in rows if not any(other['joint_relative_rmse']<=r['joint_relative_rmse'] and
                other['mean_compute_seconds']<r['mean_compute_seconds'] for other in rows)]
            frontier.sort(key=lambda r:r['joint_relative_rmse'])
            x=np.array([r['joint_relative_rmse'] for r in frontier]);y=np.array([r['mean_compute_seconds'] for r in frontier])
            low=np.array([r['joint_relative_rmse_ci95'][0] for r in frontier]);high=np.array([r['joint_relative_rmse_ci95'][1] for r in frontier])
            axis.errorbar(x,y,xerr=np.maximum(0,np.array([x-low,high-x])),fmt=marker+'-',color=color,label=label,markersize=4,capsize=2)
        for target in data['arguments']['targets']:axis.axvline(target,linestyle=':',color='black',alpha=.25)
        axis.set_xscale('log');axis.set_yscale('log');axis.grid(alpha=.2)
        axis.set_title(f"epsilon={case['epsilon']:g}")
        axis.set_xlabel('Relative trajectory RMS error of deviation')
        axis.legend(fontsize=7)
    axes[0,0].set_ylabel('Complete compute time per trajectory (s)')
    figure.suptitle('Homogeneous collisions before resampling: measured cost versus accuracy')
    output=args.summary.parent/'cost_error.png'
    figure.savefig(output,dpi=180)
    print(output)
    ratios=[];labels=[]
    for case in cases:
        for pair in case['comparisons']:
            if pair.get('complete') and pair.get('target_verified'):
                ratios.append(pair['time_ratio']);labels.append(f"eps {case['epsilon']:g}\nerror <= {pair['target']:g}\n{pair['comparison']}")
    if ratios:
        fig,axis=plt.subplots(figsize=(max(8,len(ratios)*.9),4.6),constrained_layout=True)
        bars=axis.bar(np.arange(len(ratios)),ratios,color=['#2474b5' if 'ordinary' in label else '#19814a' for label in labels])
        axis.axhline(1,color='black',linestyle='--',linewidth=1)
        axis.set_yscale('log');axis.set_xticks(np.arange(len(ratios)),labels,fontsize=8)
        axis.set_ylabel('Validated PIC time / HDP-mixing time')
        axis.set_title('Common accuracy thresholds; values > 1 favor HDP')
        for bar,value in zip(bars,ratios):axis.text(bar.get_x()+bar.get_width()/2,value,f'{value:.2g}x',ha='center',va='bottom',fontsize=8)
        output=args.summary.parent/'validated_time_ratio.png'
        fig.savefig(output,dpi=180);print(output)


if __name__=='__main__':main()
