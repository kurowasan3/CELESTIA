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

# ---- (5) rewritten: Thevenin equivalent vs. real power bookkeeping ----
E0s,r0s,Rws,isym=sp.symbols('E_0 r_0 R_w i',positive=True)
Is=(E0s+Rws*isym)/(r0s+Rws)
assert sp.simplify(E0s-(r0s*Is+Rws*(Is-isym)))==0
P0=E0s**2/(r0s+Rws); VT=Rws*E0s/(r0s+Rws); rT=r0s*Rws/(r0s+Rws)
print('(5a) E0*I - (P0 + VT*i) =',sp.simplify(E0s*Is-(P0+VT*isym)))
assert sp.simplify(E0s*Is-(P0+VT*isym))==0
# external power is the same in both descriptions: VPQ*i
VPQ_real=Rws*(Is-isym); VPQ_th=VT-rT*isym
assert sp.simplify(VPQ_real-VPQ_th)==0
internal_real=sp.simplify(r0s*Is**2+Rws*(Is-isym)**2); internal_th=rT*isym**2
print('     real internal heat - thevenin internal heat =',sp.factor(sp.simplify(internal_real-internal_th)),' (= P0, independent of i)')
assert sp.simplify(internal_real-internal_th-P0)==0
nv={E0s:E0,r0s:r0,Rws:Rw,isym:i}
print('     numbers: P0',P0.subs(nv),' VT i',(VT*isym).subs(nv),' E0 I',(E0s*Is).subs(nv),' VPQ i',(VPQ_th*isym).subs(nv),
      ' real internal',internal_real.subs(nv),' thevenin internal',internal_th.subs(nv))
print('DCC-007 (5) OK')
