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
print(' K0(eV)=',0.5*m*v0**2/q,' (>> qV0 = 100 eV)')
# symbolic check of exact sum
ph,al=sp.symbols('phi alpha')
print(sp.simplify(sp.cos(ph)-sp.cos(ph+al) - 2*sp.sin(al/2)*sp.sin(ph+al/2)))

# (6) new: P = (2/pi) asin(min(1, sigma/|sin(alpha/2)|)), sigma = d B v0/(8 V0)
sigma=d*B*v0/(8*V0); print('sigma=',sigma)
rng=np.random.default_rng(0)
def frac_num_e(n,e_,N=20000):
    f=(n+0.5)/T
    t0=rng.random(N)*50/f
    devs=np.abs(np.array([top_dev(m*(1+e_),f,t) for t in t0]))
    return np.mean(devs<=d/2)
def frac_th_e(n,e_):
    a=np.pi*(2*n+1)*e_
    return 2/np.pi*np.arcsin(min(1.0,sigma/abs(np.sin(a/2))))
n=4
for e_,label in [(1/9,'(2n+1)eps=1 (odd): minimum'),(3/9,'(2n+1)eps=3 (odd): minimum'),(2/9,'(2n+1)eps=2 (even): P=1'),(0.05,'generic'),(0.2,'generic')]:
    pn=frac_num_e(n,e_); pt=frac_th_e(n,e_)
    print(f'eps={e_:.4f} {label}: num P={pn:.4f} theory P={pt:.4f}  Pmin=(2/pi)asin(sigma)={2/np.pi*np.arcsin(sigma):.4f}')
    assert abs(pn-pt)<0.02
# P not monotone in n (fixed eps)
e_=0.03
Ps=[frac_th_e(k,e_) for k in range(0,60)]
print('P(n) for eps=0.03, n=0..59 (theory):',np.round(Ps,3))
assert not np.all(np.diff(Ps)<=0)
# coincidence: (2n+1)eps=2k  => f_n(m) = f_{n+k}(m')
nn,kk=4,1; e_=2*kk/(2*nn+1)
fn_m=(2*nn+1)/(2*T); Tm2=T*(1+e_); fiso=(2*(nn+kk)+1)/(2*Tm2)
print('f_n(m)=',fn_m,' f_{n+k}(m\')=',fiso); assert abs(fn_m/fiso-1)<1e-12
print('max |dev| of isotope at that f:',max(abs(top_dev(m*(1+e_),fn_m,t)) for t in rng.random(2000)*50/fn_m))
