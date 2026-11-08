# M50-1: Otto-type cycle without knowing gamma; max work at T_B = sqrt(T_A T_C)
import sympy as sp, numpy as np
TA,TB,TC,n,Cv,g,r=sp.symbols('T_A T_B T_C n C_V gamma r',positive=True)
TD=TA*TC/TB
# check adiabatic consistency: TB/TA = r^(g-1), TC/TD = r^(g-1)
print('TC/TD - TB/TA =',sp.simplify(TC/TD-TB/TA))
Qin=n*Cv*(TC-TB); Qout=n*Cv*(TD-TA)
W=sp.simplify(Qin-Qout); e=sp.simplify(W/Qin); print('W',W,' e',sp.factor(e))
# work by legs: adiabatic AB on gas nCv(TB-TA), CD by gas nCv(TC-TD)
Wlegs=n*Cv*(TC-TD)-n*Cv*(TB-TA); print('legs check',sp.simplify(Wlegs-W))
TBopt=sp.solve(sp.diff(W,TB),TB); print('TB opt',TBopt)
Wm=sp.simplify(W.subs(TB,sp.sqrt(TA*TC))); print('Wmax',sp.factor(Wm), sp.simplify(Wm-n*Cv*(sp.sqrt(TC)-sp.sqrt(TA))**2))
print('e at opt',sp.simplify(e.subs(TB,sp.sqrt(TA*TC))))
print('TD at opt',sp.simplify(TD.subs(TB,sp.sqrt(TA*TC))))
print('second derivative',sp.simplify(sp.diff(W,TB,2)))
# numeric check via explicit PV integration, gamma=7/5, n=1, R=8.31
R=8.31; gm=1.4; cv=R/(gm-1); TAn=300.; TCn=1200.
for TBn in [450.,600.,800.]:
    VA=1.0; VB=VA*(TAn/TBn)**(1/(gm-1))
    PA=R*TAn/VA; PB=PA*(VA/VB)**gm
    PC=PB*TCn/TBn; PD=PC*(VB/VA)**gm
    TDn=PD*VA/R
    # work by numeric integration
    Vs=np.linspace(VB,VA,200001)
    Wexp=np.trapezoid(PC*(VB/Vs)**gm,Vs); Wcomp=np.trapezoid(PA*(VA/Vs)**gm,Vs)
    Wnet=Wexp-Wcomp
    print(TBn,'TD num',TDn,'formula',TAn*TCn/TBn,'W num',Wnet,'formula',cv*(TCn-TBn-TAn*TCn/TBn+TAn),'r',VA/VB)
print('Wmax formula',cv*(np.sqrt(TCn)-np.sqrt(TAn))**2)
