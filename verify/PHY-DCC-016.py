# Q159-1 無限格子の隣接点間の合成抵抗（重ね合わせ）を大きな有限格子の数値計算で確認
import numpy as np, scipy.sparse as sps, scipy.sparse.linalg as spl
def R_between(edges,N,a,b):
    i,j=np.array(edges).T
    L=sps.coo_matrix((np.ones(len(i)),(i,j)),shape=(N,N)); L=L+L.T
    deg=np.array(L.sum(1)).ravel(); Lap=(sps.diags(deg)-L).tocsr()
    # ground node 0 (far corner)
    keep=np.arange(1,N); A=Lap[keep][:,keep]
    rhs=np.zeros(N); rhs[a]=1; rhs[b]=-1
    v=spl.spsolve(A.tocsc(),rhs[keep]); V=np.zeros(N); V[keep]=v
    return V[a]-V[b], V
def square(n):
    idx=lambda x,y:x*n+y; E=[]
    for x in range(n):
        for y in range(n):
            if x+1<n:E.append((idx(x,y),idx(x+1,y)))
            if y+1<n:E.append((idx(x,y),idx(x,y+1)))
    return E,n*n,idx(n//2,n//2),idx(n//2+1,n//2)
def triangular(n):
    idx=lambda x,y:x*n+y; E=[]
    for x in range(n):
        for y in range(n):
            if x+1<n:E.append((idx(x,y),idx(x+1,y)))
            if y+1<n:E.append((idx(x,y),idx(x,y+1)))
            if x+1<n and y+1<n:E.append((idx(x,y),idx(x+1,y+1)))
    return E,n*n,idx(n//2,n//2),idx(n//2+1,n//2)
def honeycomb(n):
    # brick-wall representation: square lattice with every other vertical bond removed
    idx=lambda x,y:x*n+y; E=[]
    for x in range(n):
        for y in range(n):
            if x+1<n:E.append((idx(x,y),idx(x+1,y)))
            if y+1<n and (x+y)%2==0:E.append((idx(x,y),idx(x,y+1)))
    return E,n*n,idx(n//2,n//2),idx(n//2+1,n//2)
def cubic(n):
    idx=lambda x,y,z:(x*n+y)*n+z; E=[]
    for x in range(n):
        for y in range(n):
            for z in range(n):
                if x+1<n:E.append((idx(x,y,z),idx(x+1,y,z)))
                if y+1<n:E.append((idx(x,y,z),idx(x,y+1,z)))
                if z+1<n:E.append((idx(x,y,z),idx(x,y,z+1)))
    c=n//2; return E,n**3,idx(c,c,c),idx(c+1,c,c)
for name,f,n,exact in [('square',square,301,1/2),('triangular',triangular,301,1/3),('honeycomb',honeycomb,301,2/3),('cubic',cubic,41,1/3)]:
    E,N,a,b=f(n); r,V=R_between(E,N,a,b)
    print(f'{name:10s} R_AB/R = {r:.5f}  (expected {exact:.5f})')
# power fraction in edge AB (square): current in AB = (V_a-V_b)/R
E,N,a,b=square(301); r,V=R_between(E,N,a,b)
print('square: current in AB =',V[a]-V[b],' power fraction AB =',(V[a]-V[b])**2/r)
# capacitor lattice: same math with admittance; energy fraction = (1/2)C V^2 / ((1/2)(2C) V^2) = 1/2
