"""Stage isolation and deterministic source-kernel integration, without resampling."""
import json
from pathlib import Path
import shutil
import subprocess
import numpy as np


def main():
    root=Path(__file__).resolve().parents[1]
    out=root/'research/runs/source_conservation_v1'
    out.mkdir(exist_ok=True)
    exe=out/'negpar_source_probe.exe'
    if not exe.exists():
        shutil.copy2(root/'build/release/Release/negpar_source_probe.exe',exe)
        for dll in (root/'build/release/Release').glob('*.dll'):shutil.copy2(dll,out/dll.name)
        shutil.copy2(Path(__file__),out/'analysis.py')
        shutil.copy2(root/'research/source_conservation.cpp',out/'probe.cpp')
        shutil.copy2(root/'src/resampling/NegativeParticleSampling.cpp',out/'sampling.cpp')
    results={'stages':[], 'quadrature':[]}
    for dt in (.01,.005):
        path=out/f'stages_dt{dt:g}.csv'
        if not path.exists():subprocess.run([str(exe),'stages',str(path),str(dt),'2048'],check=True)
        data=np.genfromtxt(path,delimiter=',',names=True)
        for fresh in (0,1):
            subset=data[data['fresh_cache']==fresh]
            row={'dt':dt,'fresh_cache':bool(fresh),'replicas':len(subset)}
            for name in ('source_mass','source_v2','transport_mass','transport_v2'):
                values=subset[name]
                row[name]={'mean':float(values.mean()),'se':float(values.std(ddof=1)/np.sqrt(len(values)))}
            combined=subset['source_v2']+subset['transport_v2']
            row['combined_v2']={'mean':float(combined.mean()),'se':float(combined.std(ddof=1)/np.sqrt(len(combined)))}
            results['stages'].append(row)
            print(row,flush=True)
    # Integrate h/M under a unit Gaussian. For a conservative Maxwellian
    # source, E[h/M-1] must vanish separately for every source velocity.
    for n in (24,40):
        nodes,weights=np.polynomial.hermite.hermgauss(n)
        xyz=np.stack(np.meshgrid(*(nodes*np.sqrt(2),)*3,indexing='ij'),axis=-1).reshape(-1,3)
        w=np.prod(np.stack(np.meshgrid(*(weights/np.sqrt(np.pi),)*3,indexing='ij'),axis=-1),axis=-1).ravel()
        for dt in (.01,.005):
            for speed in (1.,2.):
                label=f'kernel_n{n}_dt{dt:g}_speed{speed:g}'
                output=out/(label+'.csv')
                if not output.exists():
                    input_path=out/(label+'.txt')
                    np.savetxt(input_path,np.column_stack((xyz,np.tile([speed,0,0],(len(xyz),1)))))
                    subprocess.run([str(exe),'kernel',str(input_path),str(output),str(dt)],check=True)
                data=np.genfromtxt(output,delimiter=',',names=True)
                ratio=data['source']/data['maxwellian']-1
                row={'order':n,'dt':dt,'speed':speed,'source_mass_integral':float(w@ratio)}
                results['quadrature'].append(row)
                print(row,flush=True)
    (out/'summary.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
