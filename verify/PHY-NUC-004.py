# N63-1 verification: photon hypothesis vs neutral particle (Chadwick-type argument)
import sympy as sp
E,Ep,K,M,c,p=sp.symbols('E Ep K M c p',positive=True)
# photon backscatter head-on: energy & momentum conservation (nonrel recoil)
sol=sp.solve([sp.Eq(E,Ep+K), sp.Eq(E/c,-Ep/c+sp.sqrt(2*M*K))],[E,Ep],dict=True)[0]
print('E =',sp.simplify(sol[E]))   # (K+sqrt(2 M c^2 K))/2
mp=938.0; Kp=5.7; KN=1.4
Efun=lambda K,Mc2:(K+(2*Mc2*K)**0.5)/2
Ep_=Efun(Kp,mp); EN_=Efun(KN,14*mp); print('E needed (p):',Ep_,' (N):',EN_)
# predicted N recoil if E=Ep_
k=sp.symbols('k',positive=True)
kN=sp.nsolve(k+sp.sqrt(2*14*mp*k)-2*Ep_,k,0.5); print('predicted K_N:',kN)
# exact relativistic check for proton (Compton-like backscatter on proton)
Ekin=lambda E,Mc2: 2*E**2/(Mc2+2*E)   # max KE transfer, relativistic
Er=sp.nsolve(2*k**2/(mp+2*k)-Kp,k,50); print('relativistic E(p):',Er)
# neutral particle: elastic head-on, max KE to target X: 4 M X/(M+X)^2 K0
m,X,V,V1,VX=sp.symbols('m X V V1 VX',positive=True)
s=sp.solve([sp.Eq(m*V,m*V1+X*VX),sp.Eq(m*V**2,m*V1**2+X*VX**2)],[V1,VX],dict=True)
s=[t for t in s if t[VX]!=0][0]
KX=sp.simplify(X*s[VX]**2/(m*V**2)); print('K_X/K0 =',sp.factor(KX))
x=sp.symbols('x',positive=True)
R=sp.simplify((KX.subs(X,1)/KX.subs(X,14)).subs(m,x))
xs=sp.solve(sp.Eq(R,Kp/KN),x); print('M/mp candidates',xs)
xv=[float(t) for t in xs if t>0][0]
K0=Kp*(xv+1)**2/(4*xv); print('M/mp=',xv,' K0=',K0)
# closed form
r=(14*Kp/KN)**0.5; print('closed form (14-r)/(r-1)=',(14-r)/(r-1))
