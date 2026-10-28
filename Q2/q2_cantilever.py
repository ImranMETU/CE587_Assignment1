"""Exact annular cantilever. All model inputs/outputs use SI units."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize, brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
P=14000.; L=10.; RHO=7850.; SIGMA=165e6; TAU=50e6; RMAX=.4
GUARD=1e-8  # m: numerical guard, not a manufacturing constraint
LABELS=['bending','shear','outer maximum','inner maximum','inner nonnegative','positive thickness']
def geometry(x):
    ro,ri=np.asarray(x,dtype=float)
    if not (np.isfinite(ro) and np.isfinite(ri) and 0<=ri<ro):
        raise ValueError('Require finite 0 <= Ri < Ro')
    # Factored expressions reduce cancellation near equal radii.
    area=np.pi*(ro-ri)*(ro+ri)
    inertia=area*(ro*ro+ri*ri)/4
    return area,inertia
def mass(x): return RHO*L*geometry(x)[0]
def bending_stress(x):return P*L*x[0]/geometry(x)[1]
def shear_stress(x):
    ro,ri=x
    return P*(ro*ro+ro*ri+ri*ri)/(3*geometry(x)[1])
def constraints(x,guard=0.):
    ro,ri=x
    return np.array([bending_stress(x)/SIGMA-1,shear_stress(x)/TAU-1,
        ro/RMAX-1,ri/RMAX-1,-ri/RMAX,(ri+guard-ro)/RMAX])
def feasible(x,tol=1e-7,guard=0.):
    try:return bool(max(constraints(x,guard))<=tol and x[0]-x[1]>=guard)
    except ValueError:return False
def analytical_optimum():
    ro=RMAX
    return np.array([ro,(ro**4-4*P*L*ro/(np.pi*SIGMA))**.25])
def independent_check(x):
    # Expanded textbook expressions, evaluated in extended precision.
    ro,ri=np.array(x,dtype=np.longdouble)
    pi=np.longdouble(str(np.pi));A=pi*(ro**2-ri**2);I=pi*(ro**4-ri**4)/4
    return dict(area=float(A),inertia=float(I),mass=float(RHO*L*A),
        bending=float(P*L*ro/I),shear=float(P*(ro**2+ro*ri+ri**2)/(3*I)),
        bending_residual_Pa=float(P*L*ro/I-SIGMA),shear_residual_Pa=float(P*(ro**2+ro*ri+ri**2)/(3*I)-TAU),
        thickness_m=float(ro-ri))
def optimize(guard=GUARD):
    # Parameterization makes every trial section valid: x=[Ro, q*Ro],
    # q in [0, 1-guard/RMAX]. Physical Ro bounds imply no arbitrary minimum radius.
    solid_min=(4*P*L/(np.pi*SIGMA))**(1/3)
    def decode(z):return np.array([z[0]*RMAX,z[0]*RMAX*z[1]])
    def optimizer_constraints(z):
        x=decode(z);g=constraints(x,guard)
        # Equivalent capacity/demand constraints avoid singular stress gradients.
        g[:2]=1-1/(1+g[:2])
        return -g-1e-11
    rng=np.random.default_rng(5872)
    starts=[[.1,0],[.2,.1],[.35,.3],[.4,.399],[.15,.05]]
    starts += [[rng.uniform(solid_min,.4),0] for _ in range(15)]
    for j in range(5,len(starts)):starts[j][1]=starts[j][0]*rng.uniform(0,.95)
    runs=[]
    for x0 in starts:
        z0=[x0[0]/RMAX,x0[1]/x0[0]]
        res=minimize(lambda z:mass(decode(z))/1000,z0,method='SLSQP',
            bounds=[(solid_min/RMAX,1),(0,1-guard/RMAX)],
            constraints=[dict(type='ineq',fun=optimizer_constraints)],
            options=dict(ftol=1e-12,maxiter=2000))
        x=decode(res.x)
        runs.append(dict(start_m=x0,x_m=x.tolist(),mass_kg=mass(x),max_g=float(max(constraints(x,guard))),
                         success=bool(res.success),message=str(res.message),iterations=int(res.nit)))
    accepted=[r for r in runs if r['max_g']<=1e-7]
    if not accepted:raise RuntimeError('No independently feasible numerical result')
    strict=[r for r in accepted if r['max_g']<=0]
    best=min(strict,key=lambda r:r['mass_kg']);x=np.array(best['x_m']);ref=analytical_optimum()
    check=independent_check(x); g=constraints(x)
    assert feasible(x) and abs(check['mass']-mass(ref))<1e-4
    assert np.allclose(x,ref,atol=1e-8,rtol=0)
    out=dict(x_m=x.tolist(),x_cm=(x*100).tolist(),mass_kg=mass(x),analytical_x_m=ref.tolist(),
        responses=check,g=g.tolist(),active=[n for n,v in zip(LABELS,g) if abs(v)<=1e-6],
        inactive=[n for n,v in zip(LABELS,g) if abs(v)>1e-6],guard_m=guard,feasibility_tolerance=1e-7,
        activity_tolerance=1e-6,accepted_runs=len(accepted),strictly_feasible_runs=len(strict),
        solver_successes=sum(r['success'] for r in runs),
        near_global_runs=int(sum(abs(r['mass_kg']-mass(ref))<1e-4 for r in accepted)),runs=runs)
    return out
def plot():
    plt.rcParams.update({'font.size':10,'pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(11,5.2),layout='constrained')
    C=4*P*L/(np.pi*SIGMA);rmin=C**(1/3);opt=analytical_optimum()*100
    for ax,limits in zip(axes,[(0,40,0,40),(38,40,37.6,40)]):
        lo,hi,bot,top=limits;ro=np.linspace(max(lo/100,1e-7),hi/100,2400)
        bend=np.where(ro>=rmin,np.maximum(ro**4-C*ro,0)**.25,np.nan)
        shear=[]
        rs=np.sqrt(4*P/(3*np.pi*TAU))
        for r in ro:
            if r<rs:shear.append(np.nan)
            else:shear.append(brentq(lambda ri:shear_stress([r,ri])/TAU-1,0,r*(1-1e-12),xtol=1e-15))
        shear=np.array(shear);composite=np.minimum(bend,shear)
        ax.fill_between(ro*100,0,composite*100,color='#d5eadd',label='Feasible region')
        ax.plot(ro*100,bend*100,color='#2364aa',label='Bending limit',lw=2.8)
        ax.plot(ro*100,shear*100,color='#e17c23',label='Shear limit',lw=1.2)
        ax.plot(ro*100,composite*100,'k--',lw=1.1,label='Composite constraint boundary')
        ax.plot(ro*100,ro*100,':',color='#a34285',label='Ri = Ro (excluded)')
        for j,m in enumerate([334.4476849,500,1000,2500,5000]):
            y=np.sqrt(np.maximum(ro**2-m/(RHO*np.pi*L),0));valid=ro**2>=m/(RHO*np.pi*L)
            ax.plot(ro[valid]*100,y[valid]*100,color='.6',lw=.7)
            if lo==0:
                r=[.22,.28,.34,.37,.395][j];yend=np.sqrt(r*r-m/(RHO*np.pi*L))*100
                ax.annotate(f'{m:.2f} kg' if m<400 else f'{m:g} kg',(r*100,yend),xytext=(-5,-12),textcoords='offset points',ha='right',fontsize=7,color='.35',bbox=dict(facecolor='white',edgecolor='none',alpha=.8,pad=.3))
            elif j==0:
                # Keep the minimum-mass contour label visibly separate from the optimum star.
                r=.39;yend=np.sqrt(r*r-m/(RHO*np.pi*L))*100
                ax.annotate('334.448 kg',(r*100,yend),xytext=(-44,-22),textcoords='offset points',
                            ha='right',fontsize=7,color='.35',arrowprops=dict(arrowstyle='-',lw=.5,color='.45'),
                            bbox=dict(facecolor='white',edgecolor='none',alpha=.85,pad=.3))
        ax.axvline(40,color='.3');ax.axhline(40,color='.3');ax.axhline(0,color='.3')
        ax.plot(*opt,'*',color='#bd2341',ms=12,clip_on=False,zorder=10)
        ax.set(xlim=(lo,hi),ylim=(bot,top),xlabel='Outer radius Ro (cm)',ylabel='Inner radius Ri (cm)')
        ax.grid(alpha=.15)
    axes[0].set_title('Q2: exact annular cantilever')
    axes[1].set_title('Optimum and separated boundaries')
    axes[1].annotate(f'Optimum ({opt[0]:.6f}, {opt[1]:.6f}) cm\n334.447685 kg',opt,xytext=(38.12,38.1),arrowprops=dict(arrowstyle='->'),fontsize=9)
    h,l=axes[0].get_legend_handles_labels();fig.legend(h,l,loc='outside lower center',ncol=3,fontsize=9)
    fig.savefig(ROOT/'Q2_graphical.pdf');fig.savefig(ROOT/'Q2_graphical.png',dpi=200);plt.close(fig)
if __name__=='__main__':
    out=optimize();(ROOT/'verification.json').write_text(json.dumps(out,indent=2));plot()
    print(json.dumps({k:v for k,v in out.items() if k!='runs'},indent=2))
