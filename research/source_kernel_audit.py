"""Audit legacy kernel quadrature and truncated angular (Bessel) factor.

This independent vectorized translation is cross-checked against archived
C++ kernel values; variants are diagnostic, not production algorithm changes.
"""
import json
from pathlib import Path
import numpy as np


def kernel(v, source, dt, ndelta=16, cutoff_factor=1., exact_angular=False):
    u=v-source
    speed=np.linalg.norm(u,axis=1)
    r=.5*speed**3/(5*dt)
    logr=np.log(r)
    cutoff=np.exp(np.where(r>1,-.0026*logr**2-.415*logr+.9154,
                          -.0006*logr**2-.2242*logr+.4711))*cutoff_factor
    delta=cutoff[:,None]*np.arange(ndelta)[None,:]/ndelta
    zeta=(1+delta**2)**1.5
    density=np.sqrt(r[:,None]/np.pi)*zeta**1.5*np.exp(-r[:,None]*delta**2*zeta)
    density[:,[0,-1]]*=.5
    epsilon=speed[:,None]**2*delta**2/2
    perpendicular=v-(np.sum(v*u,axis=1)/speed**2)[:,None]*u
    vp=np.sum(perpendicular**2,axis=1)/2
    product=epsilon*vp[:,None]
    angular=np.i0(2*np.sqrt(product)) if exact_angular else 1+product+.25*product**2
    return 2*cutoff/ndelta*np.sum(np.exp(-epsilon)*angular*density,axis=1)-1


def main():
    out=Path(__file__).parent/'runs/source_conservation_v1'
    previous=out/'kernel_audit.json'
    results=json.loads(previous.read_text()) if previous.exists() else []
    results=[r for r in results if 'coordinates' not in r]
    for order in (() if results else (40,64)):
        nodes,weights=np.polynomial.hermite.hermgauss(order)
        v=np.stack(np.meshgrid(*(nodes*np.sqrt(2),)*3,indexing='ij'),axis=-1).reshape(-1,3)
        w=np.prod(np.stack(np.meshgrid(*(weights/np.sqrt(np.pi),)*3,indexing='ij'),axis=-1),axis=-1).ravel()
        for speed in (1.,2.):
            for name,n,cutoff,angular in [('legacy',16,1.,False),('refined_delta',128,1.,False),
                                         ('extended_delta',128,2.,False),('full_angular',128,2.,True)]:
                parts=[kernel(v[i:i+2048],np.array([speed,0,0]),.01,n,cutoff,angular)
                       for i in range(0,len(v),2048)]
                q=np.concatenate(parts)
                row=dict(order=order,speed=speed,variant=name,mass_integral=float(w@q))
                if order==40 and name=='legacy':
                    original=np.genfromtxt(out/f'kernel_n40_dt0.01_speed{speed:g}.csv',delimiter=',',names=True)
                    row['max_cpp_difference']=float(np.max(np.abs(q-(original['source']/original['maxwellian']-1))))
                    if row['max_cpp_difference']>1e-8:raise RuntimeError('Translation disagrees with C++')
                results.append(row);print(row,flush=True)
    # Gaussian product quadrature converges slowly near v=source. Resolve
    # that region explicitly in source-centered spherical coordinates.
    eps=.1/2.**np.arange(21)
    axial=np.concatenate((-3+np.arange(1,41)*(1-eps[0]+3)/40,
                          1-eps[1:],1+eps[1:],
                          3-np.arange(1,41)*(3-1-eps[0])/40))
    alpha=-float(np.min(kernel(np.column_stack((axial,np.zeros((len(axial),2)))),
                               np.array([1.,0,0]),.01)))
    for order in (64,128,256):
        radial,wr=np.polynomial.legendre.leggauss(order)
        radial=(radial+1)*8;wr=wr*8
        angle,wa=np.polynomial.legendre.leggauss(96)
        radius,cosine=np.meshgrid(radial,angle,indexing='ij')
        for speed in (1.,2.):
            v=np.column_stack(((speed+radius*cosine).ravel(),
                               (radius*np.sqrt(1-cosine**2)).ravel(),np.zeros(radius.size)))
            density=np.exp(-np.sum(v*v,axis=1)/2)/(2*np.pi)**1.5
            w=(2*np.pi*radius**2*wr[:,None]*wa[None,:]).ravel()*density
            for name,n,cutoff,angular in [('legacy',16,1.,False),('full_angular',128,2.,True)]:
                q=np.concatenate([kernel(v[i:i+2048],np.array([speed,0,0]),.01,n,cutoff,angular)
                                  for i in range(0,len(v),2048)])
                row=dict(coordinates='source_spherical',order=order,speed=speed,variant=name,
                         mass_integral=float(w@q),gaussian_integral=float(w.sum()))
                if name=='legacy':
                    row['alpha_negative']=alpha
                    row['split_omitted_mass']=float(alpha*np.sum(w[q>=alpha]))
                    row['split_mass_integral']=row['mass_integral']-row['split_omitted_mass']
                results.append(row);print(row,flush=True)
    (out/'kernel_audit.json').write_text(json.dumps(results,indent=2)+'\n')


if __name__=='__main__':main()
