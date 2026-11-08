# T20223-1: semipermeable membrane; Y = heat ballast; effective exponent; membrane rupture
import sympy as sp
n,nu,R,T,T0,V,V0,M,g,S,dh,dT,dV=sp.symbols('n nu R T T0 V V0 M g S dh dT dV',positive=True)
P=M*g/S
# (1) X partial pressure equal on both sides; piston feels only X
Veq=sp.solve(sp.Eq(P*(V+V0),n*R*T0),V)[0]; print('V =',Veq)
print('bottom force =',sp.simplify(P*S+nu*n*R*T0/V0*S))
# (2) isobaric heating
dTsol=sp.solve(sp.Eq(P*S*dh,n*R*dT),dT)[0]
Q=sp.Rational(3,2)*(1+nu)*n*R*dTsol+P*S*dh
print('Q =',sp.factor(Q)); print('nu for Q=4Mg dh:',sp.solve(sp.Eq(Q,4*M*g*dh),nu))
# (3) adiabatic: (3/2)(1+nu) n R dT = - p dVx, p=nRT/Vx
Vx=sp.symbols('V_x',positive=True)
Tf=sp.Function('T')
ode=sp.Eq(sp.Rational(3,2)*(1+nu)*Tf(Vx).diff(Vx),-Tf(Vx)/Vx)
s=sp.dsolve(ode); print(s)
kappa=1+sp.Rational(2,3)/(1+nu); print('kappa=',sp.simplify(kappa))
print('T ratio when Vx halves:',sp.simplify(2**(sp.Rational(2,3)/(1+nu))),' nu=1 ->',2**sp.Rational(1,3))
# (4) rupture: energy U+Mgh conserved, final eq (1+nu)nRTf = P*Wf
Tf_,Wf,Wi=sp.symbols('T_f W_f W_i',positive=True)
sol=sp.solve([sp.Eq(sp.Rational(3,2)*(1+nu)*n*R*(Tf_-T0),-P*(Wf-Wi)),sp.Eq((1+nu)*n*R*Tf_,P*Wf)],[Tf_,Wf],dict=True)[0]
sol={k:sp.simplify(v.subs(Wi,n*R*T0/P)) for k,v in sol.items()}
print('T_f/T0 =',sp.factor(sol[Tf_]/T0),'  W_f/W_i =',sp.factor(sol[Wf]/(n*R*T0/P)))
print('nu=1:',sol[Tf_].subs(nu,1)/T0, sol[Wf].subs(nu,1)/(n*R*T0/P))
