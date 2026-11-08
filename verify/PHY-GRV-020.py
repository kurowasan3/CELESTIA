# Q112-1: force between hemispheres of a rotating uniform asteroid
import numpy as np, sympy as sp
rng=np.random.default_rng(1)
G=M=R=1.0
def sample(n,sign):
    out=[]
    while len(out)<n:
        p=rng.uniform(-1,1,(3*n,3)); p=p[(np.linalg.norm(p,axis=1)<1)&(sign*p[:,0]>0)]
        out.extend(p.tolist())
    return np.array(out[:n])
# MC pair force (x-component on hemisphere x>0 from x<0)
est=[]
for rep in range(20):
    A=sample(3000,1); B=sample(3000,-1)
    d=A[:,None,:]-B[None,:,:]; r=np.linalg.norm(d,axis=2)
    Fx=-(G*(M/2/3000)**2*d[:,:,0]/r**3).sum()
    est.append(Fx)
print('MC attraction',-np.mean(est),'+-',np.std(est)/np.sqrt(20),' theory 3/16=',3/16)
# sympy: centroid of hemisphere and force with rotation
x,r,th,ph,w=sp.symbols('x r theta phi omega',positive=True)
Gs,Ms,Rs=sp.symbols('G M R',positive=True)
rho=Ms/(sp.Rational(4,3)*sp.pi*Rs**3)
# x = r sin th cos ph with ph in (-pi/2,pi/2)
xc=sp.integrate(r*sp.sin(th)*sp.cos(ph)*r**2*sp.sin(th),(r,0,Rs),(th,0,sp.pi),(ph,-sp.pi/2,sp.pi/2))/(sp.Rational(2,3)*sp.pi*Rs**3)
print('x_cm',sp.simplify(xc))
F=(Ms/2)*xc*(Gs*Ms/Rs**3-w**2); print('pillar force',sp.factor(F))
print('omega_c',sp.solve(F,w))
# cap check at omega_c: normal force density -(GM/R^3)x+w^2 x =0 pointwise -> any cap zero
print('pointwise',sp.simplify((-Gs*Ms/Rs**3*x+w**2*x).subs(w,sp.sqrt(Gs*Ms/Rs**3))))
# numeric example: rho=3000 kg/m^3 -> critical period
Gv=6.67e-11; rv=3000.0
Tc=2*np.pi/np.sqrt(4*np.pi*Gv*rv/3); print('critical period [h]',Tc/3600)
