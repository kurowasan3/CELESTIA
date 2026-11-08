# P69-1 rectangular frame with current near long straight wire
import sympy as sp
mu,I1,I2,L,w,r,a,F0,x=sp.symbols('mu0 I1 I2 L w r a F0 x',positive=True)
B=lambda s: mu*I1/(2*sp.pi*s)
Fnet=I2*L*(B(r)-B(r+w))           # attraction
print('F(r)=',sp.simplify(Fnet))
# (2) F(a)/F(2a)=3 -> w
sol=sp.solve(sp.Eq(Fnet.subs(r,a)/Fnet.subs(r,2*a),3),w); print('w=',sol)
wv=sol[0]
mu_sol=sp.solve(sp.Eq(Fnet.subs({r:a,w:wv}),F0),mu); print('mu0=',mu_sol)
# flux via integral
Phi=sp.integrate(B(x)*L,(x,r,r+w)); print('Phi=',sp.simplify(Phi))
print('F == -I2 dPhi/dr ?',sp.simplify(Fnet+I2*sp.diff(Phi,r))==0)
# (3) external work to move from r=a to infinity (quasi-static) = integral of F dr
W=sp.integrate(Fnet.subs(w,wv),(r,a,sp.oo)); print('W=',sp.simplify(W), ' I2*Phi(a)=',sp.simplify(I2*Phi.subs({r:a,w:wv})))
print('W in F0:',sp.simplify(W.subs(mu,mu_sol[0])))
# (4) battery-2 energy saving = int eps*I2 dt with eps=-dPhi/dt -> I2*(Phi(a)-0)
# (5) mutual inductance M=Phi/I1 ; battery-1 saving = I1*Delta(M I2) = same
M=Phi/I1; print('battery1 saving', sp.simplify(I1*(M*I2).subs({r:a,w:wv})))
# energy: Wext + decrease of M I1 I2 = savings of both batteries
print('balance', sp.simplify(W + (M*I1*I2).subs({r:a,w:wv}) - 2*I2*Phi.subs({r:a,w:wv})))
