# K20202-1: boost converter, CCM/DCM
import numpy as np, sympy as sp
E,D,T,L,R,V=sp.symbols('E D T L R V',positive=True)
# DCM charge balance
Ip=E*D*T/L; t2=L*Ip/(V-E)
sol=sp.solve(sp.Eq(sp.Rational(1,2)*Ip*t2/T, V/R),V)
print('DCM V0 solutions',sol)
Vd=E/2*(1+sp.sqrt(1+2*D**2*R*T/L))
Rb=2*L/(D*(1-D)**2*T)
print('boundary continuity:',sp.simplify(Vd.subs(R,Rb)-E/(1-D)))
num={E:10,D:sp.Rational(1,2),T:20e-6,L:100e-6}
print('Rb',float(Rb.subs(num)),'V(240)',float(Vd.subs(num).subs(R,240)))
# energy ratio
V0=30; print('load energy/cycle',V0**2/240*20e-6,'0.5LIp^2',0.5*100e-6*1.0**2,'ratio',V0/(V0-10))
# time-domain simulation (ideal switch & diode), C=20uF
def sim(Rl,C=20e-6,cycles=3000,dt=2e-8):
    Ee,Dd,Tt,Ll=10,0.5,20e-6,100e-6
    I=0.0;Vc=10.0; steps=int(round(Tt/dt)); on=int(round(Dd*steps))
    rec=[]
    for c in range(cycles):
        for s in range(steps):
            if s<on:
                I+=Ee/Ll*dt; Vc-=Vc/(Rl*C)*dt
            else:
                if I>0:
                    I+=(Ee-Vc)/Ll*dt
                    if I<0: I=0.0
                    Vc+=(I-Vc/Rl)/C*dt
                else:
                    Vc-=Vc/(Rl*C)*dt
            if c==cycles-1: rec.append((I,Vc))
    rec=np.array(rec); return rec
for Rl in [40,240]:
    r=sim(Rl); print('R',Rl,'Vavg',r[:,1].mean(),'Imin',r[:,0].min(),'Imax',r[:,0].max(),'zero-current time us',(r[:,0]==0).sum()*2e-8*1e6)
