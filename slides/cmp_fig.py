import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
NAVY='#1f3b5c'; RED='#c0392b'; BLUE='#1f77b4'; GREY='#777'; GREEN='#2a9d8f'
x=np.linspace(100,300,400)
p1=1/(1+np.exp(-(-9.19+0.046*x)))
p2=1/(1+np.exp(-(-20.5+0.203*x-0.000439*x**2)))
fig,(ax,ax2)=plt.subplots(2,1,figsize=(9,5.2),gridspec_kw={'height_ratios':[3,2.2]},sharex=True)
ax.plot(x,p1,color=BLUE,lw=3,label='Model 1: straight logistic (keeps rising)')
ax.plot(x,p2,color=RED,lw=3,ls='--',label='Model 2: logistic + quadratic term (over-curing hurts)')
ax.set_ylabel('P(pass)'); ax.set_ylim(0,1.08); ax.legend(frameon=False,fontsize=9,loc='upper left'); ax.spines[['top','right']].set_visible(False)
designs=[('Optimal for Model 1 only',[(166,.5),(233,.5)],BLUE),
         ('Optimal for Model 2 only',[(129,.333),(171,.333),(300,.333)],RED),
         ('Compound design (serves both)',[(135,.272),(171,.242),(226,.264),(300,.222)],GREEN)]
for k,(name,pts,col) in enumerate(designs):
    y=2-k
    for t,w in pts:
        ax2.plot([t,t],[y-0.02,y+w*1.1],color=col,lw=7,solid_capstyle='butt')
        ax2.text(t,y+w*1.1+0.06,f'{t}',ha='center',fontsize=8,color=col)
    ax2.axhline(y,color='#ddd',lw=1,zorder=0)
ax2.set_yticks([2.2,1.2,0.2]); ax2.set_yticklabels([d[0] for d in designs],fontsize=9,fontweight='bold')
for lab,d in zip(ax2.get_yticklabels(),designs): lab.set_color(d[2])
ax2.tick_params(axis='y',length=0); ax2.set_ylim(-0.15,2.9); ax2.set_xlim(95,305); ax2.set_xlabel('oven temperature (°C)')
ax2.spines[['top','right','left']].set_visible(False)
ax2.text(302,2.75,'bar height = weight',fontsize=8,color=GREY,ha='right')
plt.tight_layout(); plt.savefig('img/compound_designs.pdf'); plt.savefig('img/compound_designs.png',dpi=160)
