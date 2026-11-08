# S0515-1 verification: MHD generator driven by constant pressure difference with wall drag kv
import sympy as sp
dp,w,d,l,B,k,rho,R,I,v,Eb=sp.symbols('Delta_p w d ell B k rho R I v E_b',positive=True)
r=rho*d/(l*w)
# (1) open circuit
v0=dp*w*d/k; V0=v0*B*d; print('v0=',v0,' V0=',sp.simplify(V0))
# (2) with R: unknowns v, I
sol=sp.solve([sp.Eq(dp*w*d, k*v+I*B*d), sp.Eq(v*B*d, I*(R+r))],[v,I],dict=True)[0]
print('v=',sp.simplify(sol[v]),'\nI=',sp.simplify(sol[I]))
V=sp.simplify(R*sol[I])
eps0=V0; rp=r+B**2*d**2/k
print('check V = eps0 - r\' I :',sp.simplify(V-(eps0-rp*sol[I])))
# (3) maximise P = R I^2
P=sp.simplify(R*sol[I]**2)
Rs=sp.solve(sp.diff(P,R),R); print('R*=',Rs)
Rst=Rs[0]; print('R*-r\' =',sp.simplify(Rst-rp))
Pmax=sp.simplify(P.subs(R,Rst)); print('Pmax=',Pmax,' check eps0^2/(4r\')',sp.simplify(Pmax-eps0**2/(4*rp)))
# (4) efficiency at R*
Pin=dp*w*d*sol[v].subs(R,Rst)
eta=sp.simplify(Pmax/Pin); print('eta=',eta, ' expected a/(2(2r+a)):',sp.simplify(eta-(B**2*d**2/k)/(2*(2*r+B**2*d**2/k))))
# energy balance
Ib=sol[I]; vb=sol[v]
print('power balance:',sp.simplify(dp*w*d*vb-(R+r)*Ib**2-k*vb**2))
# (5) pump: dp=0, battery Eb opposing generator EMF (upper electrode +)
Ir=sp.Symbol('I_r',real=True); vr=sp.Symbol('v_r',real=True)
s5=sp.solve([sp.Eq(0, k*vr+Ir*B*d), sp.Eq(vr*B*d-Eb, Ir*r)],[vr,Ir],dict=True)
print('pump:',s5)
vp=Eb*B*d/(B**2*d**2+k*r)
print('check vp:',sp.simplify(s5[0][vr]-vp) if s5 else None)
