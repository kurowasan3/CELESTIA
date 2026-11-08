# P02-1 quasi-static stick-slip of a rod on two supports with different materials
import sympy as sp, numpy as np
# --- quasi-static simulation: positions on rod measured from end A; G at l
def simulate(l,p0,q0,mua,muap,mub,mubp,steps=200000,dl=None):
    p,q=p0,q0; dl=dl or (q0-p0)/steps*0.5
    slipping=None; stops=[]
    while q-p>1e-4*(q0-p0):
        x=l-p; y=q-l           # distances from G
        Na=y/(x+y); Nb=x/(x+y)  # in units Mg
        # decide which support slides relative to rod
        if slipping is None:
            slipping='b' if mub*Nb<mua*Na else 'a'
        if slipping=='b':
            # b sliding: continue while kinetic at b < static max at a
            if mubp*Nb<mua*Na: q-=dl
            else: stops.append(('b',q)); slipping='a'
        else:
            if muap*Na<mub*Nb: p+=dl
            else: stops.append(('a',p)); slipping='b'
        if len(stops)>=4: break
    return stops
l=0.37; mua,muap,mub,mubp=0.8,0.6,0.3,0.25
st=simulate(l,0.0,1.0,mua,muap,mub,mubp)
print(st)
# formula: x0=l, y0=1-l ; first slip b iff mub*x0<mua*y0
x0,y0=l,1-l; print('b first?',mub*x0<mua*y0)
y1=mubp/mua*x0; x1=muap/mub*y1; y2=mubp/mua*x1
print('pred q1,p2,q3',l+y1,l-x1,l+y2)
q1,p2,q3=st[0][1],st[1][1],st[2][1]
print('l from formula',q1*p2/(q1+p2-q3))
# symbolic check of inversion
L,rb,ra,lsym=sp.symbols('L r_b r_a l',positive=True)
Q1=lsym+rb*lsym; P2=lsym-ra*rb*lsym; Q3=lsym+rb*ra*rb*lsym
print(sp.simplify(Q1*P2/(Q1+P2-Q3)-lsym))
print('rb',sp.simplify((Q1-lsym)/lsym),'ra',sp.simplify((lsym-P2)/(Q1-lsym)))
# example: nearer support slides first
print('case near slips: x0=0.6,y0=0.4, mua=0.8, mub=0.3 ->',0.3*0.6<0.8*0.4)
print('--- problem numbers: mua=1, mub=0.8, mubp=0.5, muap=0.4, l=0.4L')
st=simulate(0.4,0.0,1.0,1.0,0.4,0.8,0.5,steps=400000)
print(st, 'expected q1=.6 p2=.3 q3=.45 p4=.375')
q1,p2,q3=0.6,0.3,0.45
print('l',q1*p2/(q1+p2-q3),'rb',(q1-q3)/p2,'ra',p2*(q3-p2)/(q1*(q1-q3)))
q1s,p2s,q3s=sp.symbols('q1 p2 q3',positive=True)
ll=q1s*p2s/(q1s+p2s-q3s)
print('rb',sp.simplify((q1s-ll)/ll-(q1s-q3s)/p2s),'ra',sp.simplify((ll-p2s)/(q1s-ll)-p2s*(q3s-p2s)/(q1s*(q1s-q3s))))
