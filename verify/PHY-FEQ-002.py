# T20171-1: Jenga-like tower, min-cut brute force over moving sets (max-flow min-cut for static friction network)
import itertools, sympy as sp
def analyse(mu1,mu2,w=0.0,M=1.0,g=1.0):
    # blocks: L1:a0,a1,a2 (along x), L2:b0,b1,b2 (along x, on top of a_i), L3:c0,c1,c2 (along y, crossing all b)
    blocks=['a0','a1','a2','b0','b1','b2','c0','c1','c2']
    # normal forces
    edges=[]  # (u,v,N,mu) v may be 'floor'
    Nc=(M+w/3)*g/3   # each crossing c_j-b_i
    for i in range(3):
        for j in range(3): edges.append((f'b{i}',f'c{j}',Nc,mu1))
        Nb=M*g+3*Nc          # b_i on a_i
        edges.append((f'a{i}',f'b{i}',Nb,mu1))
        edges.append((f'a{i}','floor',M*g+Nb,mu2))
    best=None
    others=[b for b in blocks if b!='b1']
    for r in range(len(others)+1):
        for sub in itertools.combinations(others,r):
            S=set(sub)|{'b1'}
            cap=sum(N*mu for u,v,N,mu in edges if (u in S)!=(v in S))
            if best is None or cap<best[0]-1e-12: best=(cap,S)
    return best
for mu1,mu2 in [(0.6,0.3),(0.6,0.5),(0.6,0.39),(0.6,0.41)]:
    c,S=analyse(mu1,mu2); print(mu1,mu2,'F=',round(c,6),sorted(S),' theory',(mu1+3*mu2) if 3*mu2<2*mu1 else 3*mu1)
# weight on top, mu2=0.8mu1: threshold w=6M
for w in [5.0,5.9,6.1,8.0]:
    mu1=0.5; mu2=0.4
    c,S=analyse(mu1,mu2,w); 
    FA=mu1*(3+2*w/3); FB=mu1*(1+w/3)+mu2*(3+w/3)
    print('w=',w,'F=',round(c,6),sorted(S),' FA',FA,' FB',FB)
m1,m2,w,M=sp.symbols('mu1 mu2 w M',positive=True)
print('threshold w:',sp.solve(sp.Eq(m2*(3*M+w/3),m1*(2*M+w/3)),w))
