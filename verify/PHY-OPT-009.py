# N08-1: Michelson with two wavelengths -> visibility beats
import sympy as sp, numpy as np
x,l1,l2=sp.symbols('x lambda1 lambda2',positive=True)
I=(1+sp.cos(4*sp.pi*x/l1))+(1+sp.cos(4*sp.pi*x/l2))
I2=2*(1+sp.cos(2*sp.pi*x*(1/l1+1/l2))*sp.cos(2*sp.pi*x*(1/l1-1/l2)))
print("identity:",sp.simplify(sp.expand_trig(sp.expand(I-I2))) , [float((I-I2).subs({x:v,l1:0.589,l2:0.5896})) for v in [0.1,7.3,123.4]])
L=289.4e3; N=982  # nm
lam1=4*L/(2*N+1); lam2=4*L/(2*N-1); print("lam1,lam2,dlam:",lam1,lam2,lam2-lam1, "Delta x check", lam1*lam2/(2*(lam2-lam1)))
print("first vanishing for Na 589.0/589.6 (mm):",589.0*589.6/(4*0.6)*1e-6)
# numeric: count maxima between successive vanishing points
xs=np.linspace(0,0.7e6,14000001)  # nm
Iv=np.cos(4*np.pi*xs/lam1)+np.cos(4*np.pi*xs/lam2)
env=np.abs(np.cos(2*np.pi*xs*(1/lam1-1/lam2)))
z=xs[1:-1][(env[1:-1]<env[:-2])&(env[1:-1]<env[2:])]; print("vanishing x(nm):",z[:3])
seg=(xs>z[0])&(xs<z[1]); s=Iv[seg]; pk=np.sum((s[1:-1]>s[:-2])&(s[1:-1]>s[2:])); print("maxima between:",pk)
# unequal intensities 9:1 near envelope minimum
I1,I2_=9,1
xx=np.linspace(z[0]-2000,z[0]+2000,400001)
Iu=I1*(1+np.cos(4*np.pi*xx/lam1))+I2_*(1+np.cos(4*np.pi*xx/lam2))
print("max/min near vanishing:",Iu.max()/Iu.min())
