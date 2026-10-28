"""Q3 tripod: N, cm, kg. No effective-length or safety factors beyond assignment."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize, minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
W=60000.; B=120.; E=7.5e6; RHO=.0028; SIGMA=15000.; FS=2.
BOUNDS=((50.,200.),(.5,10.))
NAMES=['compression','Euler buckling','minimum height','maximum height','minimum diameter','maximum diameter']
def strut_length(H):return np.sqrt(np.asarray(H)**2+B**2/3)
def area(D):return np.pi*np.asarray(D)**2/4
def inertia(D):return np.pi*np.asarray(D)**4/64
def axial_force(H):return W*strut_length(H)/(3*np.asarray(H))
def compressive_stress(x):H,D=x;return axial_force(H)/area(D)
def critical_load(x):H,D=x;return np.pi**2*E*inertia(D)/strut_length(H)**2
def mass(x):H,D=x;return 3*RHO*area(D)*strut_length(H)
def responses(x):
    H,D=np.asarray(x,dtype=float)
    if not(np.isfinite(H) and np.isfinite(D) and H>0 and D>0):raise ValueError('H,D must be finite and positive')
    return dict(length_cm=float(strut_length(H)),area_cm2=float(area(D)),inertia_cm4=float(inertia(D)),
        force_N=float(axial_force(H)),stress_N_cm2=float(compressive_stress(x)),pcr_N=float(critical_load(x)),mass_kg=float(mass(x)))
def utilization(x):
    responses(x) # validity check
    H,D=x
    return np.array([compressive_stress(x)/SIGMA,FS*axial_force(H)/critical_load(x),50/H,H/200,.5/D,D/10])
def constraints(x):return utilization(x)-1
def feasible(x,tol=1e-7):
    try:return bool(np.max(constraints(x))<=tol)
    except ValueError:return False
def diameter_compression(H):return np.sqrt(4*W*strut_length(H)/(3*np.pi*SIGMA*np.asarray(H)))
def diameter_buckling(H):return (64*FS*W*strut_length(H)**3/(3*np.pi**3*E*np.asarray(H)))**.25
def composite_diameter(H):return np.maximum.reduce([np.full_like(np.asarray(H,dtype=float),.5),diameter_compression(H),diameter_buckling(H)])
def analytical_optimum():return np.array([50.,float(composite_diameter(50.))])
def independent_check(x):
    """Recompute entirely in N-m-kg and from explicit 3D base coordinates."""
    H,D=np.asarray(x)/100
    b=B/100;radius=b/np.sqrt(3)
    feet=np.array([[radius*np.cos(a),radius*np.sin(a),0] for a in [0,2*np.pi/3,4*np.pi/3]])
    apex=np.array([0,0,H]);vectors=apex-feet;lengths=np.linalg.norm(vectors,axis=1)
    # Solve 3D equilibrium independently of the supplied scalar force expression.
    forces=np.linalg.solve((vectors/lengths[:,None]).T,np.array([0.,0.,W]))
    A=np.pi*D**2/4;I=np.pi*D**4/64;l=lengths[0]
    pcr=np.pi**2*75e9*I/l**2;stress=forces[0]/A;m=2800*A*sum(lengths)
    return dict(length_cm=float(l*100),area_cm2=float(A*1e4),inertia_cm4=float(I*1e8),force_N=float(forces[0]),
        stress_N_cm2=float(stress/1e4),pcr_N=float(pcr),mass_kg=float(m),
        base_sides_m=[float(np.linalg.norm(feet[i]-feet[(i+1)%3])) for i in range(3)],
        equilibrium_residual_N=((vectors/lengths[:,None]).T@forces-np.array([0,0,W])).tolist(),
        compression_residual_N_cm2=float(stress/1e4-SIGMA),buckling_residual_N=float(FS*forces[0]-pcr))
def optimize():
    rng=np.random.default_rng(5873)
    starts=[[50,.5],[200,10],[100,5],[50,4],[200,1]]+[[rng.uniform(50,200),rng.uniform(.5,10)] for _ in range(20)]
    runs=[]
    scale=np.array([100.,5.])
    for start in starts:
        r=minimize(lambda z:mass(z*scale),np.array(start)/scale,method='SLSQP',
            bounds=[(.5,2),(.1,2)],constraints=[dict(type='ineq',fun=lambda z:-constraints(z*scale))],
            options=dict(ftol=1e-12,maxiter=2000))
        x=r.x*scale
        runs.append(dict(start_cm=start,x_cm=x.tolist(),mass_kg=float(mass(x)),g=constraints(x).tolist(),
            success=bool(r.success),message=str(r.message),iterations=int(r.nit)))
    accepted=[r for r in runs if max(r['g'])<=1e-7 and abs(r['mass_kg']-mass(analytical_optimum()))<=1e-6]
    if not accepted:raise RuntimeError('No verified global candidate')
    # Select the smallest residual among candidates rather than favor tiny violations.
    best=min(accepted,key=lambda r:max(r['g']));x=np.array(best['x_cm']);ind=independent_check(x)
    for k,v in responses(x).items():assert np.isclose(v,ind[k],rtol=1e-12,atol=1e-10),(k,v,ind[k])
    reduced=minimize_scalar(lambda H:mass([H,composite_diameter(H)]),bounds=(50,200),method='bounded',options={'xatol':1e-11})
    g=constraints(x)
    out=dict(x_cm=x.tolist(),responses=responses(x),independent=ind,g=g.tolist(),utilization=utilization(x).tolist(),
        active=[n for n,v in zip(NAMES,g) if abs(v)<=1e-6],inactive=[n for n,v in zip(NAMES,g) if abs(v)>1e-6],
        analytic_x_cm=analytical_optimum().tolist(),reduced_scalar_x_cm=[float(reduced.x),float(composite_diameter(reduced.x))],
        reduced_mass_kg=float(reduced.fun),accepted_runs=len(accepted),solver_successes=sum(r['success'] for r in runs),
        feasibility_tolerance=1e-7,activity_tolerance=1e-6,runs=runs)
    return out
def plot():
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(11,5),layout='constrained')
    opt=analytical_optimum()
    for ax,limits in zip(axes,[(50,200,.5,10),(49,90,1,4.6)]):
        xmin,xmax,ymin,ymax=limits;H=np.linspace(max(xmin,50),xmax,1500)
        ds=diameter_compression(H);db=diameter_buckling(H);cc=composite_diameter(H)
        ax.fill_between(H,cc,10,color='#d5eadd',label='Feasible region')
        ax.fill_between(H,ds,db,color='#fff0d0',label='Strength satisfied; buckling fails')
        ax.plot(H,ds,color='#2364aa',label='Compression boundary',lw=2)
        ax.plot(H,db,color='#e17c23',label='Euler / FS boundary',lw=2.8)
        ax.plot(H,cc,'k--',lw=1.2,label='Composite boundary (buckling)')
        for m in [6.6054480656,10,20,40,80]:
            d=np.sqrt(4*m/(3*RHO*np.pi*strut_length(H)))
            ax.plot(H,d,color='.6',ls=':',lw=.8)
            idx=np.argmin(abs(H-(80 if xmin==50 else 85)))
            if ymin<d[idx]<ymax:ax.text(H[idx],d[idx]+.04,f'{m:.3f} kg' if m<7 else f'{m:g} kg',fontsize=8,color='.35')
        for h in [50,200]:ax.axvline(h,color='.3',lw=.8)
        for d in [.5,10]:ax.axhline(d,color='.3',lw=.8)
        ax.plot(*opt,'*',color='#bd2341',ms=12,clip_on=False,zorder=9)
        ax.set(xlim=(xmin,xmax),ylim=(ymin,ymax),xlabel='Height H (cm)',ylabel='Strut diameter D (cm)');ax.grid(alpha=.15)
    axes[0].set_title('Q3: minimum-mass tripod')
    axes[1].set_title('Governing boundary near the optimum')
    axes[1].annotate(f'H = 50 cm, D = {opt[1]:.6f} cm\nEuler buckling + minimum height',xy=opt,xytext=(57,4.1),arrowprops=dict(arrowstyle='->'),fontsize=9)
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=3,fontsize=8)
    fig.savefig(ROOT/'Q3_graphical.pdf');fig.savefig(ROOT/'Q3_graphical.png',dpi=200);plt.close(fig)
if __name__=='__main__':
    r=optimize();(ROOT/'verification.json').write_text(json.dumps(r,indent=2));plot()
    print(json.dumps({k:v for k,v in r.items() if k!='runs'},indent=2))
