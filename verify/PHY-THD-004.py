# M52-1: orientation-dependent isobaric heating; determine M from timing; switch to isochoric at rim
import sympy as sp
P0,S,L,M,g,H,R,Ta,t1,t2,t=sp.symbols('P0 S L M g H R T_a tau1 tau2 t',positive=True)
x=sp.symbols('x',positive=True)  # x = Mg/S
Pu=P0+x; Pd=P0-x
tau1=sp.Rational(5,2)*Pu*S*L/3/H; tau2=sp.Rational(5,2)*Pd*S*L/3/H
xs=sp.solve(sp.Eq(tau2/tau1,sp.Rational(1,2)),x)[0]; print('Mg/S =',xs)
Pu_,Pd_=Pu.subs(x,xs),Pd.subs(x,xs)
Ta_=Pu_*S*L/3/R; Tb_=Pd_*S*L/3/R; print('Tb/Ta',sp.simplify(Tb_/Ta_))
tau1v=tau1.subs(x,xs)
t_rim=2*tau2.subs(x,xs); print('t_rim/tau1',sp.simplify(t_rim/tau1v))
T_rim=Pd_*S*L/R; print('T_rim/Ta',sp.simplify(T_rim/Ta_))
t_star=t_rim+sp.Rational(3,2)*S*L*(Pu_-Pd_)/H; print('t*/tau1',sp.simplify(t_star/tau1v))
T_star=Pu_*S*L/R; print('T*/Ta',sp.simplify(T_star/Ta_))
Tt=T_rim+H*(t-t_rim)/(sp.Rational(3,2)*R); print('T(t)/Ta =',sp.simplify((Tt/Ta_).subs(t,sp.Symbol('s')*tau1v)))
# energy check over whole exp2 to t*
Qtot=H*t_star; dU=sp.Rational(3,2)*R*(T_star-Tb_); Wg=Pd_*S*(2*L/3)
print('energy residual',sp.simplify(Qtot-dU-Wg))
# numeric independent stepping: march heat input in small dt with state switching
import numpy as np
P0n,Sn,Ln,Rn=1.0e5,1.0e-2,0.6,8.31; xn=P0n/3; Hn=50.0
V=Sn*Ln/3; T=(P0n-xn)*V/Rn; dt=1e-3; tt=0.0; Tb=T
while True:
    dQ=Hn*dt
    if V<Sn*Ln-1e-15:
        P=P0n-xn; dT=dQ/(2.5*Rn); V=min(Rn*(T+dT)/P,Sn*Ln); T=T+dT
    else:
        T+=dQ/(1.5*Rn)
    tt+=dt
    if V>=Sn*Ln and Rn*T/V>=P0n+xn: break
tau1n=2.5*(P0n+xn)*Sn*Ln/3/Hn; Tan=(P0n+xn)*Sn*Ln/3/Rn
print('numeric t*/tau1',tt/tau1n,' T*/Ta',T/Tan,' Tb/Ta',Tb/Tan)
