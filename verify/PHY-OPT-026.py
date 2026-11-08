# J086-1: isosceles prism (apex alpha, index n) on liquid (n_L); horizontal ray into AB, reflection at base, exit AC
import numpy as np, sympy as sp
def refract(dvec,nrm,n1,n2):
    # nrm: unit normal pointing into incident side (against dvec)
    c=-np.dot(nrm,dvec); r=n1/n2; k=1-r*r*(1-c*c)
    if k<0: return None
    return r*dvec+(r*c-np.sqrt(k))*nrm
def trace(alpha,n,nL,yin=0.25):
    H=1.0; half=H*np.tan(alpha/2)  # apex at (0,H), base from (-half,0) to (half,0)
    B=np.array([-half,0]);A=np.array([0,H]);C=np.array([half,0])
    # entry point on AB at height yin*H
    P=B+(A-B)*yin; d=np.array([1.0,0])
    tAB=(A-B)/np.linalg.norm(A-B); nAB=np.array([-tAB[1],tAB[0]])  # outward normal (left-up)
    inc=np.degrees(np.arccos(abs(np.dot(d,nAB))))
    d1=refract(d,nAB,1.0,n)
    th1=np.degrees(np.arccos(abs(np.dot(d1,nAB))))
    # hit base y=0
    t=-P[1]/d1[1]; Q=P+t*d1
    th2=np.degrees(np.arccos(abs(d1[1])))
    tir = n*np.sin(np.radians(th2))>nL
    d2=np.array([d1[0],-d1[1]])
    # hit AC: line through C and A
    tAC=(A-C)/np.linalg.norm(A-C); nAC=np.array([tAC[1],-tAC[0]])  # outward right-up
    # solve Q+s d2 = C + w tAC
    M=np.array([d2,-tAC]).T; s,w=np.linalg.solve(M,C-Q)
    d3=refract(d2,-nAC,n,1.0)
    return inc,th1,th2,tir,Q,d3,w
alpha=np.radians(90); n=1.50
inc,th1,th2,tir,Q,d3,w=trace(alpha,n,1.33,0.2)
print('inc',inc,' th1',th1,' formula',np.degrees(np.arcsin(np.sin(alpha/2)/n)),' th2',th2, 'TIR',tir,' exit dir',d3, 'Q',Q,'w on AC',w)
N=lambda n,a: np.cos(a/2)*np.sqrt(n*n-np.sin(a/2)**2)+np.sin(a/2)**2
print('n sin th2 =',n*np.sin(np.radians(th2)),' formula N',N(n,alpha))
# inverse: n* = 1.47 -> n
ns=1.42; a=alpha
nrec=np.sqrt(((ns-np.sin(a/2)**2)/np.cos(a/2))**2+np.sin(a/2)**2); print('n from n*=1.42:',nrec, ' check N',N(nrec,a))
# monotonic in n
print([round(N(x,a),4) for x in (1.3,1.4,1.5,1.6)])
# symbolic identity n cos(a/2 - th1)
nn,aa=sp.symbols('n a',positive=True); th=sp.asin(sp.sin(aa/2)/nn)
print(sp.simplify(sp.expand_trig(nn*sp.cos(aa/2-th)) - (sp.cos(aa/2)*sp.sqrt(nn**2-sp.sin(aa/2)**2)+sp.sin(aa/2)**2)))
# red vs violet example
for nc in (1.51,1.53):
    print(nc, N(nc,a))
# alpha=90 glass 1.5
print('alpha90 n1.5',N(1.5,np.pi/2))
# range of entry heights that hit base before AC
for yin in (0.1,0.2,0.3,0.4,0.45,0.5):
    r=trace(alpha,n,1.33,yin); print(yin,'Q x',r[4][0], ' within base?',abs(r[4][0])<np.tan(alpha/2),' exit dir',r[5],' w on AC (0..len)',r[6])
