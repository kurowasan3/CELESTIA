# J080-1: string tension via hanging mass in accelerating / coasting train; beats with tuning fork
import sympy as sp, numpy as np
M,g,a,L,rho,f0,n,beta=sp.symbols('M g a L rho f0 n beta',positive=True)
fs=sp.sqrt(M*g/rho)/(2*L)
# equilibrium of hanging mass in train frame: T vector balances M g_eff, g_eff=(−a, −g) -> T=M sqrt(g^2+a^2)
T=M*sp.sqrt(g**2+a**2); fa=sp.sqrt(T/rho)/(2*L)
ratio=sp.simplify(fa/fs); print('fa/fs=',ratio)
asol=sp.solve(sp.Eq(ratio,(f0+n)/(f0-n)),a); print('a=',[sp.simplify(x) for x in asol])
# incline coasting frictionless: train accel g sin(beta) down slope, g_eff = g - a_train -> magnitude g cos(beta)
gv=sp.Matrix([0,-g]); acc=g*sp.sin(beta)*sp.Matrix([sp.cos(beta),-sp.sin(beta)])
geff=gv-acc; print('|g_eff|=',sp.simplify(sp.sqrt(geff.dot(geff))))
fb=(f0-n)*sp.sqrt(sp.cos(beta)); print('beats on slope=',f0-fb)
# numeric example
vals={f0:440,n:4,g:9.8}
print('a numeric',[x.subs(vals).evalf() for x in asol])
# check that beats go through zero when a increases from 0 (frequency passes f0)
