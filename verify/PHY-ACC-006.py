import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
# lossless symbolic: half-cycle map
E,C,L=sp.symbols('E C L',positive=True)
V=0; s=1; peaks=[]
for n in range(1,7):
    center=s*E; V=center-(V-center); peaks.append(sp.simplify(V)); s=-s
print('V_C at zero-current instants:',peaks)
# energy given by battery in interval n (after n-th flip): increase of C energy
n=sp.symbols('n',integer=True,positive=True)
print('work:',sp.factor(C*((2*(n+1)*E)**2-(2*n*E)**2)/2), ' direct E*dq:',sp.factor(E*C*(2*n*E+2*(n+1)*E)))
print('max current after nth flip:', (2*n+1)*E*sp.sqrt(C/L))
def run(Lv,Cv,Ev,r,nflip,tend=None):
    def rhs(t,y,sg): q,i=y; return [i,(sg*Ev-q/Cv-r*i)/Lv]
    y=[0.0,0.0]; t=0.0; sg=1; out=[]; imax=[]
    for k in range(nflip):
        def ev(tt,yy,sgg): return yy[1]
        ev.terminal=True; ev.direction=-sg
        s=solve_ivp(rhs,[t,t+10*np.sqrt(Lv*Cv)],y,args=(sg,),events=ev,rtol=1e-11,atol=1e-14,max_step=np.sqrt(Lv*Cv)/50)
        # skip the trivial event at start
        te=[tt for tt in s.t_events[0] if tt>t+1e-6*np.sqrt(Lv*Cv)]
        if not te:
            raise SystemExit('no event')
        s=solve_ivp(rhs,[t,te[0]],y,args=(sg,),rtol=1e-11,atol=1e-14,max_step=np.sqrt(Lv*Cv)/200)
        imax.append(np.max(np.abs(s.y[1])))
        t=s.t[-1]; y=[s.y[0,-1],0.0]; out.append((t,y[0]/Cv)); sg=-sg
    return out,imax
Lv,Cv,Ev=2.0,0.5,1.0
out,imax=run(Lv,Cv,Ev,0.0,5)
print('lossless num (t_n/pi sqrt(LC), V_C):',[(round(a/np.pi/np.sqrt(Lv*Cv),6),round(b,6)) for a,b in out])
print('imax num',np.round(imax,6),' formula',[(2*k+1)*Ev*np.sqrt(Cv/Lv) for k in range(5)])
# lossy steady state
r=0.05
out,_=run(Lv,Cv,Ev,r,400)
Vp=abs(out[-1][1])
g=r/(2*Lv); wd=np.sqrt(1/(Lv*Cv)-g**2); kappa=np.exp(-g*np.pi/wd)
print('lossy steady Vp num',Vp,' E(1+k)/(1-k)',Ev*(1+kappa)/(1-kappa),' approx 4E/(pi r)sqrt(L/C)',4*Ev/(np.pi*r)*np.sqrt(Lv/Cv))
# check kappa: free decay (no flip) ratio of successive deviations at zero-current instants
def rhs0(t,y): q,i=y; return [i,(Ev-q/Cv-r*i)/Lv]
def e0(t,y): return y[1]
s=solve_ivp(rhs0,[0,40],[0,0],events=e0,rtol=1e-12,atol=1e-14,max_step=0.01)
d=[abs(y[0]/Cv-Ev) for y in s.y_events[0][1:5]]
print('free-decay ratios',[d[k+1]/d[k] for k in range(3)],' kappa',kappa)
