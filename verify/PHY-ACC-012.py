# O20212-1 : line loss with distributed capacitance; reactor compensation at the load end
import sympy as sp, numpy as np
P,V,w,C,R,L,j,N,k=sp.symbols('P V omega C R L j N k',positive=True)
Ir=2*P/V; Ic=w*C*V
lump=R*(Ir**2+Ic**2)
# distributed (N sections): reactive current in section k (counted from load) = Ic*k/N - IL
dist=sp.limit(R/N*sp.summation(Ir**2+(Ic*k/N)**2,(k,1,N)),N,sp.oo); print('distributed loss:',sp.simplify(dist))
Vopt=sp.solve(sp.diff(dist,V),V); print('Vopt:',Vopt, ' Vopt^2:',sp.simplify(Vopt[0]**2),' min loss:',sp.simplify(dist.subs(V,Vopt[0])))
print('lumped Vopt^2:',sp.simplify(sp.solve(sp.diff(lump,V),V)[0]**2),' min',sp.simplify(lump.subs(V,sp.solve(sp.diff(lump,V),V)[0])))
IL=V/(w*L)
distL=sp.limit(R/N*sp.summation(Ir**2+(Ic*k/N-IL)**2,(k,1,N)),N,sp.oo); distL=sp.simplify(distL)
Lopt=sp.solve(sp.diff(distL,L),L); print('distributed with reactor:',distL,'  Lopt',Lopt,'  loss at Lopt',sp.simplify(distL.subs(L,Lopt[0])))
print('distributed loss at lumped-optimal L=1/(w^2C):',sp.simplify(distL.subs(L,1/(w**2*C))))
# numeric telegraph-ladder check (complex phasors, amplitude convention, avg power = |I|^2 R/2 per resistor)
def ladder(Pv,Vv,wv,Cv,Rv,Lv=None,Ns=2000):
    r=2*Pv/Vv**0  # placeholder
    rload=Vv**2/(2*Pv)
    Vn=Vv+0j; I=Vn/rload
    if Lv: I+=Vn/(1j*wv*Lv)
    loss=0.0
    for i in range(Ns):
        I+=1j*wv*(Cv/Ns)*Vn           # shunt capacitance of this section
        loss+=2*(Rv/Ns)*abs(I)**2/2     # two wires
        Vn+=2*(Rv/Ns)*I                 # voltage drop over both wires
    return loss,abs(Vn)
Pv,Vv,wv,Cv,Rv=1.0e9,5.0e5,2*np.pi*60,1.0e-5,1.0
lo,_=ladder(Pv,Vv,wv,Cv,Rv); print('ladder loss',lo,' formula',float(dist.subs({P:Pv,V:Vv,w:wv,C:Cv,R:Rv})),' lumped',float(lump.subs({P:Pv,V:Vv,w:wv,C:Cv,R:Rv})))
for Lv in [1/(wv**2*Cv),2/(wv**2*Cv),1.5/(wv**2*Cv),3/(wv**2*Cv)]:
    lo,_=ladder(Pv,Vv,wv,Cv,Rv,Lv); print('L=%.3f/(w^2C): ladder'%(Lv*wv**2*Cv),lo,' formula',float(distL.subs({P:Pv,V:Vv,w:wv,C:Cv,R:Rv,L:Lv})))
# original-style numbers, 10 ohm per wire: reactive fraction
Rv=10.0
print('10ohm: lumped',float(lump.subs({P:Pv,V:Vv,w:wv,C:Cv,R:Rv})),' distributed',float(dist.subs({P:Pv,V:Vv,w:wv,C:Cv,R:Rv})))
