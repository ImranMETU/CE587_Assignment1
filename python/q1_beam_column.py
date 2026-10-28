"""Q1 thin-wall model. Public structural functions take SI dimensions (m)."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent.parent
P=50000.;L=3.;E=210e9;RHO=7800.;GRAVITY=9.8;SIGMA=250e6
RMIN=.02;RMAX=.20;TMIN=.001;TMAX=.01;RATIO_MAX=20.;DELTA_MAX=.001*L
NAMES=['axial stress','self-weight deflection','R/t','minimum radius','maximum radius','minimum thickness','maximum thickness']
def cm_to_m(x):return np.asarray(x,dtype=float)/100
def m_to_cm(x):return np.asarray(x,dtype=float)*100
def _valid(x):
    R,t=np.asarray(x,dtype=float)
    if not(np.isfinite(R) and np.isfinite(t) and R>0 and t>0):raise ValueError('R,t must be finite and positive')
    return R,t
def area(x):R,t=_valid(x);return 2*np.pi*R*t
def inertia(x):R,t=_valid(x);return np.pi*R**3*t
def mass(x):return RHO*L*area(x)
def axial_stress(x):return P/area(x)
def self_weight(x):return RHO*GRAVITY*area(x)
def deflection(x):return 5*self_weight(x)*L**4/(384*E*inertia(x))
def radius_thickness_ratio(x):R,t=_valid(x);return R/t
def responses(x):
    return dict(area_m2=float(area(x)),inertia_m4=float(inertia(x)),mass_kg=float(mass(x)),
        stress_Pa=float(axial_stress(x)),self_weight_N_m=float(self_weight(x)),
        deflection_m=float(deflection(x)),radius_thickness_ratio=float(radius_thickness_ratio(x)))
def utilization(x):
    R,t=_valid(x)
    return np.array([axial_stress(x)/SIGMA,deflection(x)/DELTA_MAX,R/(RATIO_MAX*t),RMIN/R,R/RMAX,TMIN/t,t/TMAX])
def constraints(x):return utilization(x)-1
def feasible(x,tol=1e-7):
    try:return bool(np.max(constraints(x))<=tol)
    except ValueError:return False
def activity(x,tol=1e-6):return {n:('active' if abs(g)<=tol else 'inactive') for n,g in zip(NAMES,constraints(x))}
def analytic():
    k=P/(2*np.pi*SIGMA);Rstar=np.sqrt(RATIO_MAX*k);tstar=Rstar/RATIO_MAX
    rd=np.sqrt(5*RHO*GRAVITY*L**4/(192*E*DELTA_MAX))
    return dict(product_m2=k,product_cm2=k*1e4,mass_kg=RHO*L*P/SIGMA,
        optimal_R_interval_cm=[RMIN*100,Rstar*100],representative_m=[Rstar,tstar],
        representative_cm=[Rstar*100,tstar*100],deflection_radius_limit_cm=rd*100,
        mass_coefficient_kg_cm2=2*np.pi*RHO*L*1e-4)
def independent_check(x):
    # Independent N-cm geometry and unsimplified beam equation.
    R,t=m_to_cm(x);l=300.;a=2*np.pi*R*t;i=np.pi*R**3*t
    weight_cm=(7800/1e6)*a*9.8
    delta_cm=5*weight_cm*l**4/(384*(210e9/1e4)*i)
    return dict(area_m2=float(a/1e4),inertia_m4=float(i/1e8),mass_kg=float(3*(7800/1e6)*a*100),
        stress_Pa=float((50000/a)*1e4),self_weight_N_m=float(weight_cm*100),
        deflection_m=float(delta_cm/100),radius_thickness_ratio=float(R/t),
        simplified_deflection_m=float(5*RHO*GRAVITY*L**4/(192*E*x[0]**2)),
        axial_residual_Pa=float((50000/a)*1e4-SIGMA),deflection_residual_m=float(delta_cm/100-DELTA_MAX))
def optimize():
    a=analytic();rng=np.random.default_rng(5871)
    starts=[a['representative_cm'],[2,.2],[2.1,.16],[2.3,.15],[5,.5],[15,.9],[2,.1],[20,.1],[20,1],
            [2,1],[3,.11]]+[[rng.uniform(2,20),rng.uniform(.1,1)] for _ in range(14)]
    runs=[]
    for start in starts:
        r=minimize(lambda z:mass(cm_to_m(z)),start,method='SLSQP',bounds=[(2,20),(.1,1)],
            constraints=[dict(type='ineq',fun=lambda z:-constraints(cm_to_m(z)))],
            options=dict(ftol=1e-12,maxiter=1500))
        x=cm_to_m(r.x);g=constraints(x);m=float(mass(x))
        runs.append(dict(start_cm=start,initial_feasible=feasible(cm_to_m(start)),x_cm=r.x.tolist(),mass_kg=m,
            g=g.tolist(),mass_error_kg=m-a['mass_kg'],Rt_error_cm2=float(np.prod(r.x)-a['product_cm2']),
            on_analytic_segment=bool(a['optimal_R_interval_cm'][0]-1e-6<=r.x[0]<=a['optimal_R_interval_cm'][1]+1e-6),
            feasible=feasible(x),success=bool(r.success),message=str(r.message),iterations=int(r.nit)))
    accepted=[r for r in runs if r['feasible'] and abs(r['mass_error_kg'])<1e-6 and r['on_analytic_segment']]
    if not accepted:raise RuntimeError('No analytically verified optimizer candidate')
    # Choose the run nearest the requested representative; no uniqueness implied.
    rep=np.array(a['representative_cm']);best=min(accepted,key=lambda r:np.linalg.norm(np.array(r['x_cm'])-rep))
    x=cm_to_m(best['x_cm']);ind=independent_check(x)
    for k,v in responses(x).items():assert np.isclose(v,ind[k],rtol=1e-12,atol=1e-14),(k,v,ind[k])
    sample_R=np.linspace(*a['optimal_R_interval_cm'],1001);segment=[constraints(cm_to_m([R,a['product_cm2']/R])) for R in sample_R]
    out=dict(analytic=a,numerical_x_cm=best['x_cm'],responses=responses(x),independent=ind,
        utilization=utilization(x).tolist(),g=constraints(x).tolist(),activity=activity(x),
        accepted_runs=len(accepted),total_runs=len(runs),solver_successes=sum(r['success'] for r in runs),
        numerical_optimal_R_range_cm=[min(r['x_cm'][0] for r in accepted),max(r['x_cm'][0] for r in accepted)],
        segment_sample_count=1001,segment_max_g=float(np.max(segment)),feasibility_tolerance=1e-7,activity_tolerance=1e-6,runs=runs)
    return out
def plot():
    a=analytic();k=a['product_cm2'];c=a['mass_coefficient_kg_cm2'];rd=a['deflection_radius_limit_cm'];rep=a['representative_cm']
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(11,5.1),layout='constrained')
    for idx,(ax,limits) in enumerate(zip(axes,[(2,20,.1,1.),(1.5,3.3,.085,.235)])):
        xmin,xmax,ymin,ymax=limits;R=np.linspace(xmin,xmax,2400)
        lower=np.maximum.reduce([np.full_like(R,.1),k/R,R/20])
        inside=(R>=2)&(R<=20)&(lower<=1)
        ax.fill_between(R,lower,1,where=inside,color='#d5eadd',label='Feasible region')
        ax.plot(R,k/R,color='#2364aa',lw=2,label='Axial-stress boundary')
        ax.plot(R,R/20,color='#e17c23',lw=2,label='R/t = 20')
        if idx == 1:
            ax.axvline(rd,color='#a34285',ls=':',lw=1.6,label='Deflection limit (redundant)')
        else:
            ax.text(.035,.955,f'Deflection limit: $R={rd:.2f}$ cm $<R_{{\\min}}$',
                    transform=ax.transAxes,fontsize=8,color='#87386e',va='top')
        for v in [2,20]:ax.axvline(v,color='.35',lw=.8)
        for v in [.1,1]:ax.axhline(v,color='.35',lw=.8)
        ax.plot(R[inside],lower[inside],'k--',lw=1.2,label='Composite constraint boundary')
        for j,m in enumerate([4.68,5,6,8]):
            ax.plot(R,m/(c*R),color='.6',ls=':',lw=.9)
            if idx==1:
                pos=[1.74,2.08,2.55,2.92][j];y=m/(c*pos)
                ax.text(pos,y+.002,f'{m:g} kg',fontsize=8,color='.4',bbox=dict(facecolor='white',edgecolor='none',alpha=.75,pad=.3))
        optR=np.linspace(2,rep[0],500)
        ax.plot(optR,k/optR,color='#bd2341',lw=4,label='Entire optimal segment: 4.68 kg')
        ax.plot(*rep,'*',ms=13,color='#bd2341',zorder=9)
        ax.set(xlim=(xmin,xmax),ylim=(ymin,ymax),xlabel='Radius R (cm)',ylabel='Wall thickness t (cm)');ax.grid(alpha=.15)
    axes[0].set_title('Q1: thin-walled beam-column')
    axes[1].set_title('Nonunique minimum and representative design')
    axes[1].annotate('Every red point is optimal\nRt = 1/pi cm², mass = 4.68 kg',xy=(2.23,k/2.23),xytext=(2.25,.215),
        fontsize=8,arrowprops=dict(arrowstyle='->',color='.3'))
    axes[1].annotate(f'({rep[0]:.6f}, {rep[1]:.6f}) cm\nAxial stress + R/t active',xy=rep,xytext=(2.55,.101),
        fontsize=8,arrowprops=dict(arrowstyle='->',color='.3'))
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=3,fontsize=8)
    target=ROOT/'figures';target.mkdir(exist_ok=True)
    fig.savefig(target/'q1_graphical.pdf');fig.savefig(target/'q1_graphical.png',dpi=200);plt.close(fig)
if __name__=='__main__':
    out=optimize();d=ROOT/'verification';d.mkdir(exist_ok=True)
    (d/'q1_python.json').write_text(json.dumps(out,indent=2));plot()
    print(json.dumps({k:v for k,v in out.items() if k!='runs'},indent=2))
