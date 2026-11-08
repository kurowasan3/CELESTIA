# N32-1 (revised): potentiometer, wire PQ 100.0cm 20.0ohm, driver E0=7.50V r0 unknown (+ at P).
# standard 1.50V balances at 25.0cm, X balances at 40.0cm; then S at Q, ammeter (1.00 ohm) reads 0.600A.
import sympy as sp
E0=sp.Rational(15,2); Rw=20; L=100; Es=sp.Rational(3,2)
k=Es/25; drop=k*L; I0=drop/Rw; r0=(E0-drop)/I0
Ex=k*40
print('k=',k,'wire drop=',drop,'I0=',I0,'r0=',r0,'Ex=',Ex)
# S at Q: Kirchhoff with unknown r and I
I,r=sp.symbols('I r'); i=sp.Rational(3,5); rA=1
sol=sp.solve([sp.Eq(E0, r0*I+Rw*(I-i)), sp.Eq(Rw*(I-i), Ex+(r+rA)*i)],[I,r],dict=True)[0]
print('I=',sol[I],' r=',sol[r])
# consistency: Thevenin  i=(drop-Ex)/(Rw||r0 + rA + r)
Rth=Rw*r0/(Rw+r0); print('Thevenin i=',(drop-Ex)/(Rth+rA+sol[r]))
# energy
Pchem=Ex*i; PJx=sol[r]*i**2; Psup=E0*sol[I]
print('Pchem=',Pchem,float(Pchem),' Joule_X=',PJx,' E0 supplies=',Psup,float(Psup),' ratio=',float(Pchem/Psup))
# total balance
Ptot=r0*sol[I]**2+Rw*(sol[I]-i)**2+(sol[r]+rA)*i**2+Pchem
print('energy balance check:',sp.simplify(Ptot-Psup))
# direction before contact: potentials (P=0): Q at -drop, X's negative terminal (via A) at -Ex  -> A-side higher -> current from X side to Q
print('V_Q=',-drop,' V_Xminus=',-Ex)
