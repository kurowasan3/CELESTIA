# PHY-MAG-025 交流電圧で2回加減速する周回型質量分析と同位体の除去率
import numpy as np, sympy as sp
q=1.60e-19; u=1.66e-27; B=0.50; v0=2.0e5; V0=100.0; d=1.0e-3; eps=0.01
m=100*u
T=2*np.pi*m/(q*B)
def top_dev(mass,f,t0):
    """exact: enter X center at t0 with speed v0; return top-point deviation after 2 crossings"""
    Tm=2*np.pi*mass/(q*B)
    r0=mass*v0/(q*B)
    K=0.5*mass*v0**2
    t=t0+Tm/2           # first crossing (half turn)
    K+=q*V0*np.cos(2*np.pi*f*t)
    t+=Tm
    K+=q*V0*np.cos(2*np.pi*f*t)
    v=np.sqrt(2*K/mass); r=mass*v/(q*B)
    return 2*(r-r0)
# (2) for species m, f_n = (n+1/2)/T : all phases return
for n in [0,3,8]:
    f=(n+0.5)/T
    devs=[abs(top_dev(m,f,t0)) for t0 in np.random.rand(2000)*10*T]
    print('n',n,'f=%.4e'%f,'max dev for m =',max(devs))
# (4) fraction for m'=m(1+eps)
def frac_num(n,N=200000):
    f=(n+0.5)/T
    t0=np.random.rand(N)*1/f*50
    devs=np.abs(np.array([top_dev(m*(1+eps),f,t) for t in t0[:20000]]))
    return np.mean(devs<=d/2)
def frac_th(n):
    a=np.pi*(2*n+1)*eps
    s=d*B*v0/(8*V0*np.sin(a/2))
    return 1.0 if s>=1 else 2/np.pi*np.arcsin(s)
for n in [0,4,7,8,10]:
    print('n',n,'num',frac_num(n),'theory',frac_th(n))
print('threshold 2n+1 >=',2/np.pi/eps*np.arcsin(d*B*v0/(4*V0)))
print('f_8 =',(8.5)/T,' K0(eV)=',0.5*m*v0**2/q,' r0=',m*v0/(q*B))
# symbolic check of exact sum
ph,al=sp.symbols('phi alpha')
print(sp.simplify(sp.cos(ph)-sp.cos(ph+al) - 2*sp.sin(al/2)*sp.sin(ph+al/2)))
