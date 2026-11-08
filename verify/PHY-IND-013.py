# P72-1 rigid ladder moving in B=b x (z dir); rails resistance 0; rungs at x, x+a, x+3a with R, 2R, R
import sympy as sp
b,d,v,a,R,x,C,beta,B0=sp.symbols('b d v a R x C beta B0',positive=True)
pos=[x,x+a,x+3*a]; Rs=[R,2*R,R]
U=sp.Symbol('U')   # V_S - V_N ; rung emf pushes from N to S, magnitude b*x_i*d*v (v x B = -y)
eps=[b*p*d*v for p in pos]
I=[(e-U)/Ri for e,Ri in zip(eps,Rs)]   # current N->S through rung
Us=sp.solve(sp.Eq(sum(I),0),U)[0]
I=[sp.simplify(i.subs(U,Us)) for i in I]; print('rung currents N->S:',I)
# force on rungs: current along -y, F = I*d*(-y) x (B z) = -I d B x_hat
Fx=sp.simplify(sum(-i*d*b*p for i,p in zip(I,pos))); print('magnetic Fx on frame =',Fx,' -> external f =',sp.simplify(-Fx))
P=sp.simplify(sum(i**2*Ri for i,Ri in zip(I,Rs))); print('Joule =',P,' f*v =',sp.simplify(-Fx*v), sp.simplify(P+Fx*v)==0)
# (3) stationary frame in uniform B(t)=B0+beta t : loop emfs by Faraday, mesh analysis
# mesh currents: J1 loop(rung1-rung2), J2 loop(rung2-rung3); flux increasing in +z -> emf clockwise? use generic: emf_k = beta*area_k, direction such that induced current in rung_left flows N->S ... solve numerically equivalently:
J1,J2=sp.symbols('J1 J2')
# moving-case loop emf (rung j minus rung i) = b d v (xj-xi) ; static case = beta * d*(xj-xi). Kirchhoff identical => same currents iff beta=b v
e1=beta*d*a; e2=beta*d*2*a
# define loop orientation: rung left current S->N = -J1 ... use KVL around loop1 with mesh J1 flowing N->S in right rung, S->N in left rung
eqs=[sp.Eq(J1*R + (J1-J2)*2*R, e1), sp.Eq((J2-J1)*2*R + J2*R, e2)]
s=sp.solve(eqs,[J1,J2]); print('static mesh:',s)
# rung currents N->S: rung1 = -J1, rung2 = J1-J2, rung3 = J2
stat=[-s[J1], s[J1]-s[J2], s[J2]]
print('static rung currents with beta=b v:',[sp.simplify(c.subs(beta,b*v)) for c in stat])
print('match:',all(sp.simplify(c.subs(beta,b*v)-i)==0 for c,i in zip(stat,I)))
print('static net force (uniform B):',sp.simplify(sum(-c*d*(B0) for c in stat)))
# (5) rung 2 replaced by capacitor C: currents only via rungs 1,3
U2=(eps[0]/R+eps[2]/R)/(2/R); Vc=sp.simplify(eps[1]-U2); print('cap voltage (emf2 - U) =',Vc,' Q=',sp.simplify(C*Vc))
I13=sp.simplify((eps[2]-eps[0])/(2*R)); print('I13=',I13,' f=',sp.simplify(b*d*3*a*I13))
