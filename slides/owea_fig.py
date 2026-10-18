import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
NAVY='#1f3b5c'; RED='#c0392b'; BLUE='#1f77b4'; GREY='#777'
a,b=-9.19,0.046
X=np.arange(100,301,1.0)
def f(x):
    p=1/(1+np.exp(-(a+b*x))); return np.sqrt(p*(1-p))*np.array([1.0,x])
F=np.array([f(x) for x in X])            # N x 2
def M(S,w): return sum(wi*np.outer(f(s),f(s)) for s,wi in zip(S,w))
from scipy.optimize import minimize
def opt_weights(S):
    # optimal weights for a FIXED support: maximise log det M(w) over the simplex (Newton-type solver)
    k=len(S)
    obj=lambda w: -np.linalg.slogdet(M(S,w)+1e-12*np.eye(2))[1]
    res=minimize(obj, np.ones(k)/k, method='SLSQP', bounds=[(0,1)]*k,
                 constraints=[{'type':'eq','fun':lambda w: w.sum()-1}], options={'ftol':1e-12,'maxiter':500})
    w=np.maximum(res.x,0); w[w<1e-3]=0; w/=w.sum(); return w
def sens(S,w):
    Mi=np.linalg.inv(M(S,w)); return np.einsum('ij,jk,ik->i',F,Mi,F)-2
S=[100.0,300.0]; hist=[]
for it in range(12):
    w=opt_weights(S); d=sens(S,w); hist.append((list(S),w.copy(),d.copy()))
    if d.max()<1e-4: break
    S=sorted((set(s for s,wi in zip(S,w) if wi>0))|{X[d.argmax()]})
print(len(hist),'iterations')
for S,w,d in hist: print([f'{s:.0f}:{wi:.3f}' for s,wi in zip(S,w) if wi>1e-3], 'max d=%.3f'%d.max())
# figure: 2x2 panels, iterations 0,1,2,last
sel=[0,1,2,len(hist)-1]
fig,axes=plt.subplots(2,2,figsize=(10,6.2),sharex=True)
for ax,k in zip(axes.flat,sel):
    S,w,d=hist[k]
    ymax=max(1.0,d.max()); lo=-0.42*ymax; hi=1.28*ymax
    if k==len(hist)-1: lo=d.min()*1.35; hi=0.6; ymax=abs(lo)/0.42
    ax.set_ylim(lo,hi); ax.set_xlim(100,300)
    ax.axhline(0,color=GREY,lw=1)
    ax.plot(X,d,color=BLUE,lw=2.2)
    imax=d.argmax()
    for s_,wi in zip(S,w):
        if wi>1e-3:
            ax.vlines(s_,lo,lo+wi*0.36*ymax,color=NAVY,lw=7,alpha=.95)
            ax.text(s_,lo+wi*0.36*ymax+0.03*ymax,f'{wi:.2f}',ha='center',fontsize=9,color=NAVY,fontweight='bold')
    if k<len(hist)-1:
        ax.plot(X[imax],d[imax],'o',color=RED,ms=9)
        ax.annotate('largest $d(x)$: add %d°C'%X[imax],(X[imax],d[imax]),xytext=(0,12),textcoords='offset points',ha='center',color=RED,fontsize=10)
        title='Iteration %d — support {%s}'%(k,", ".join(f"{s_:.0f}" for s_,wi in zip(S,w) if wi>1e-3))
    else:
        title='Iteration %d — {%s}: $d(x)\\leq 0$ everywhere ⇒ optimal'%(k,", ".join(f"{s_:.0f}" for s_,wi in zip(S,w) if wi>1e-3))
    ax.set_title(title,fontsize=11,color=NAVY,loc='left')
    ax.spines[['top','right']].set_visible(False)
for ax in axes[1]: ax.set_xlabel('temperature (°C)')
fig.text(0.012,0.5,'sensitivity function $d(x)$      (bars = optimal weights of the current support)',rotation=90,va='center',fontsize=10,color=GREY)
plt.tight_layout(rect=(0.03,0,1,1)); plt.savefig('img/owea_iterations.pdf'); plt.savefig('img/owea_iterations.png',dpi=160)
plt.close()

# figure: logistic curve with optimal design vs napkin
fig,ax=plt.subplots(figsize=(8,4.2))
xx=np.linspace(100,300,400); p=1/(1+np.exp(-(a+b*xx)))
ax.plot(xx,p,color=BLUE,lw=3,label='assumed pass-rate curve')
for x,wt in [(166,.5),(233,.5)]:
    ax.vlines(x,0,1/(1+np.exp(-(a+b*x))),color=NAVY,lw=2,ls='--')
    ax.plot(x,1/(1+np.exp(-(a+b*x))),'o',color=NAVY,ms=12)
    ax.text(x,1/(1+np.exp(-(a+b*x)))+0.08,f'{x}°C\n50%',ha='center',color=NAVY,fontsize=10,fontweight='bold')
for x in (100,200,300):
    ax.plot(x,1/(1+np.exp(-(a+b*x))),'s',color=RED,ms=10,alpha=.85)
ax.plot([],[],'s',color=RED,label='"napkin" plan: 100/200/300, 40 each')
ax.plot([],[],'o',color=NAVY,label='D-optimal: 166 & 233, 60 each')
ax.set_xlabel('oven temperature (°C)'); ax.set_ylabel('P(pass)'); ax.set_ylim(-0.02,1.15)
ax.spines[['top','right']].set_visible(False); ax.legend(frameon=False,loc='upper left',fontsize=9)
plt.tight_layout(); plt.savefig('img/coating_design.pdf'); plt.savefig('img/coating_design.png',dpi=160)

# efficiency bar
fig,ax=plt.subplots(figsize=(7,3.2))
names=['D-optimal\n166 / 233','9 equally\nspaced','5 equally\nspaced','100/200/300\n(napkin)','100 & 300\nonly']
eff=[1.0,0.73,0.67,0.50,0.20]
cols=[NAVY,GREY,GREY,RED,GREY]
bars=ax.barh(names[::-1],eff[::-1],color=cols[::-1])
for bname,e in zip(bars,eff[::-1]): ax.text(e+0.02,bname.get_y()+bname.get_height()/2,f'{e:.0%}  → needs {120/e:.0f} pieces',va='center',fontsize=10)
ax.set_xlim(0,1.75); ax.set_xlabel('D-efficiency relative to the optimal design'); ax.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig('img/efficiency.pdf'); plt.savefig('img/efficiency.png',dpi=160)
