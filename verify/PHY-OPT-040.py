# S0409-1 complementary gratings (Babinet): numerical sum over many point sources
import numpy as np, sympy as sp
d=1.0; N=40; M=400
def amp(wf,s,start=0.0):
    """amplitude in direction with sin(theta)/lambda = s; open part [start, start+wf*d] in each period"""
    tot=0
    for k in range(N):
        x=k*d+start+(np.arange(M)+0.5)/M*wf*d
        tot+=np.sum(np.exp(2j*np.pi*s*x))*(wf*d/M)
    return tot
w=1/3
G=np.array([amp(w,m) for m in range(5)]); Gc=np.array([amp(1-w,m,start=w) for m in range(5)])
print('G : I_m/I_0',np.round(abs(G)**2/abs(G[0])**2,5))
print('G\': I_m/I_0',np.round(abs(Gc)**2/abs(Gc[0])**2,5))
print('|G_m|,|G\'_m|,|G_m+G\'_m| :',[(round(abs(a),4),round(abs(b),4),round(abs(a+b),8)) for a,b in zip(G,Gc)])
for s in [0.37,1.21]:
    print('generic s=',s,' |G|=',abs(amp(w,s)),' |G\'|=',abs(amp(1-w,s,start=w)),' open plate=',abs(amp(1.0,s)))
x=sp.symbols('x')
print('I1/I0 w=d/3:',sp.simplify((sp.sin(sp.pi/3)/(sp.pi/3))**2),'  w=2d/3:',sp.simplify((sp.sin(2*sp.pi/3)/(2*sp.pi/3))**2))
print('I2/I1 :',sp.simplify((sp.sin(2*x)/(2*x))**2/(sp.sin(x)/x)**2))
print('cos^2 x = 1/4 :',sp.solve(sp.Eq(sp.cos(x)**2,sp.Rational(1,4)),x))
print('numeric 27/(4pi^2)=',27/(4*np.pi**2),' 27/(16pi^2)=',27/(16*np.pi**2))
