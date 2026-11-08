# S0513-1 verification: piecewise-linear transistor amplifier
import sympy as sp, numpy as np
Vin,Vm,g,RA,RC,RE,beta,E,IB=sp.symbols('V_in V_m g R_A R_C R_E beta E I_B',positive=True)
# active region, no R_E
IBs=sp.solve(sp.Eq(Vin, RA*IB + Vm + IB/g), IB)[0]
Vout=E-RC*beta*IBs
K=sp.simplify(-sp.diff(Vout,Vin)); print('K =',K)
V2=sp.solve(sp.Eq(Vout,0),Vin)[0]; print('V2 =',sp.simplify(V2), ' check Vm+E/K:',sp.simplify(V2-(Vm+E/K)))
# (3) max undistorted amplitude: Vm+vs<=VB<=V2-vs
vsmax=sp.simplify((V2-Vm)/2); print('vs_max =',vsmax, ' VB =',sp.simplify(Vm+vsmax))
# (3b) VB = Vm + E/(4K), vs = E/(2K): fraction of period with Vin<Vm
th=np.linspace(0,2*np.pi,2000001)
Kn=5.0;En=10.0;Vmn=0.6
VB=Vmn+En/(4*Kn); vs=En/(2*Kn)
vin=VB+vs*np.sin(th)
def vout(v):
    return np.clip(En-Kn*(v-Vmn),0,En)
vo=vout(vin)
print('fraction stuck at E:',np.mean(vin<Vmn),' fraction at 0:',np.mean(vin>Vmn+En/Kn))
# brute force: max vs without clipping over VB grid
best=0;bestVB=None
for VBt in np.linspace(Vmn,Vmn+En/Kn,2001):
    vsm=min(VBt-Vmn,Vmn+En/Kn-VBt)
    if vsm>best: best=vsm;bestVB=VBt
print('numeric vs_max',best,' expected',En/(2*Kn),' VB',bestVB,' expected',Vmn+En/(2*Kn))
# (4) with emitter resistor
IC=beta*IB
IB4=sp.solve(sp.Eq(Vin, RA*IB + Vm + IB/g + RE*(IB+IC)), IB)[0]
VC=E-RC*beta*IB4
K4=sp.simplify(-sp.diff(VC,Vin)); print('K4 =',K4)
print('beta->oo :',sp.limit(K4,beta,sp.oo))
# direct numeric check of K4 with numbers by solving full nonlinear (active) equations
vals={g:0.02,RA:1000,RC:2000,RE:100,beta:150,E:12,Vm:0.6}
f=sp.lambdify(Vin,VC.subs(vals))
print('K4 numeric (finite diff):',-(f(1.01)-f(1.0))/0.01,' formula:',float(K4.subs(vals)))
