# M53-1: rectangle in P-T diagram (two isobars, two isotherms) -> orientation in P-V
import sympy as sp, numpy as np
R,T0,P0=sp.symbols('R T0 P0',positive=True)
ln2=sp.log(2)
st={'A':(T0,P0),'B':(2*T0,P0),'C':(2*T0,2*P0),'D':(T0,2*P0)}
V={k:R*T/P for k,(T,P) in st.items()}; print({k:sp.simplify(v/(R*T0/P0)) for k,v in V.items()})
def Q(a,b,kind):
    Ta,Pa=st[a];Tb,Pb=st[b]
    if kind=='p': return sp.Rational(5,2)*R*(Tb-Ta)
    if kind=='T': return R*Ta*sp.log(V[b]/V[a])
legs=[('A','B','p'),('B','C','T'),('C','D','p'),('D','A','T')]
Qs=[sp.simplify(Q(*l)) for l in legs]; print('Q legs',Qs)
Wnet=sp.simplify(sum(Qs)); print('W net',Wnet)
# signed area in PV (shoelace on numerical curve)
Rn,T0n,P0n=8.31,300.,1e5
def path():
    pts=[]
    T=np.linspace(T0n,2*T0n,2000); pts+= [(Rn*t/P0n,P0n) for t in T]
    P=np.linspace(P0n,2*P0n,2000); pts+= [(Rn*2*T0n/p,p) for p in P]
    T=np.linspace(2*T0n,T0n,2000); pts+= [(Rn*t/(2*P0n),2*P0n) for t in T]
    P=np.linspace(2*P0n,P0n,2000); pts+= [(Rn*T0n/p,p) for p in P]
    return np.array(pts)
p=path(); x,y=p[:,0],p[:,1]
Wn=np.sum((x[1:]-x[:-1])*(y[1:]+y[:-1])/2)  # integral P dV along path
print('numeric W',Wn,' formula',float(Wnet.subs({R:Rn,T0:T0n})))
# reversed cycle
Qr=[-q for q in Qs]; Qin=sum(q for q in Qr if q.subs({R:1,T0:1}).evalf()>0); W=-Wnet
e=sp.simplify(W/Qin); print('Qin',sp.simplify(Qin),' e',e,float(e),' Carnot',0.5)
