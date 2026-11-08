# Q162-1 半無限導体の4探針法と側面による鏡像
import sympy as sp, numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl
a,b,rho,I=sp.symbols('a b rho I',positive=True)
c=rho*I/(2*sp.pi)
def Vpt(x,y):  # potential at surface point (x,y) ; source at (0,b0) sink at (3a,b0) with face at y=0, electrodes at y=b
    pass
# electrodes on line y=b, face at y=0 -> images at y=-b
src=[(0,b,1),(0,-b,1),(3*a,b,-1),(3*a,-b,-1)]
def V(x,y,imgs=True):
    s=0
    for (X,Y,q) in (src if imgs else src[0::2]):
        s+=q/sp.sqrt((x-X)**2+(y-Y)**2)
    return c*s
V23=sp.simplify(V(a,b)-V(2*a,b)); V23n=sp.simplify(V(a,b,False)-V(2*a,b,False))
print('no face: V23 =',V23n,' -> rho =',sp.solve(sp.Eq(V23n,sp.Symbol('V')),rho))
ratio=sp.simplify(V23/V23n); print('ratio =',ratio)
print('b=a:',sp.nsimplify(ratio.subs(b,a)),float(ratio.subs(b,a)),' b->0:',sp.limit(ratio,b,0),' b->oo:',sp.limit(ratio,b,sp.oo))
# normal derivative on face y=0 vanishes (3D potential incl. depth z)
x,y,z=sp.symbols('x y z',real=True)
phi=sum(q/sp.sqrt((x-X)**2+(y-Y)**2+z**2) for (X,Y,q) in src)
print('dphi/dy at y=0:',sp.simplify(sp.diff(phi,y).subs(y,0)))
# lattice simulation (resistor cubic lattice, unit resistors -> rho = 1*h, h=1)
def sim(Nx,Ny,Nz,bb,aa):
    idx=lambda i,j,k:(i*Ny+j)*Nz+k; N=Nx*Ny*Nz
    I_,J_=[],[]
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                n=idx(i,j,k)
                if i+1<Nx: I_.append(n);J_.append(idx(i+1,j,k))
                if j+1<Ny: I_.append(n);J_.append(idx(i,j+1,k))
                if k+1<Nz: I_.append(n);J_.append(idx(i,j,k+1))
    I_=np.array(I_);J_=np.array(J_)
    Adj=sps.coo_matrix((np.ones(len(I_)),(I_,J_)),shape=(N,N)); Adj=Adj+Adj.T
    Lp=(sps.diags(np.array(Adj.sum(1)).ravel())-Adj).tocsr()
    x0=(Nx-3*aa)//2
    s=idx(x0,bb,0); t=idx(x0+3*aa,bb,0); p2=idx(x0+aa,bb,0); p3=idx(x0+2*aa,bb,0)
    rhs=np.zeros(N); rhs[s]=1; rhs[t]=-1
    keep=np.arange(1,N)
    v=spl.cg(Lp[keep][:,keep],rhs[keep],rtol=1e-10,maxiter=20000)[0]
    Vv=np.zeros(N); Vv[keep]=v
    return Vv[p2]-Vv[p3]
aa=8
v_face=sim(90,70,45,aa,aa); v_far=sim(90,70,45,35,aa)
th=lambda bb: 1+2*aa/np.sqrt(aa**2+4*bb**2)-aa/np.sqrt(aa**2+bb**2)
print('lattice ratio V(b=a)/V(b=35) =',v_face/v_far,' theory =',th(aa)/th(35))
print('lattice V(b=35)*2*pi*a (should ~1 x th(35) for rho=1):',v_far*2*np.pi*aa, th(35))
# 格子の自由境界は実効的に半格子間隔外側にあるので b_eff = a + 0.5 で比較
print('theory with b_eff=a+0.5:',th(aa+0.5)/th(35.5))
