# J088-1: two facing concave mirrors (focal f each), vertex distance D; object at lower vertex; image after upper->lower reflections at upper vertex (hole)
import sympy as sp, numpy as np
D,f=sp.symbols('D f',positive=True)
b1=f*D/(D-f); a2=D-b1; b2=sp.simplify(f*a2/(a2-f))
print('b2(D)=',sp.factor(b2))
print('solutions b2=D:',sp.solve(sp.Eq(b2,D),D), ' limit D->f of b2:',sp.limit(b2,D,f))
m=sp.simplify((-b1/D)*(-b2/a2)); print('m(D)=',sp.factor(m),' at 3f',m.subs(D,3*f),' limit D->f',sp.limit(m,D,f))
print('intermediate image at D=3f: b1=',b1.subs(D,3*f))
# unfolded paraxial matrix check (mirror = lens f), propagate D between
def sysmat(Dv,fv):
    T=lambda L: np.array([[1,L],[0,1]]); Lm=lambda F: np.array([[1,0],[-1/F,1]])
    return T(Dv)@Lm(fv)@T(Dv)@Lm(fv)@T(Dv)   # object -> upper -> lower -> back to upper vertex plane
for Dv in (1.0,3.0,2.0):
    M=sysmat(Dv,1.0); print('D/f=',Dv,' B element (0 means imaging at hole)',M[0,1],' A (magnification)',M[0,0])
# exact ray tracing with spherical mirrors R=2f, small object height and finite apertures
def trace_exact(Dv,R,h0,slope):
    # axis z upward; lower mirror vertex z=0 center at z=R ; upper mirror vertex z=Dv center at z=Dv-R
    p=np.array([h0,0.0]); dvec=np.array([slope,1.0]); dvec/=np.linalg.norm(dvec)
    for (zc,sign) in ((Dv-R,+1),(R,-1)):
        # intersect sphere |p - (0,zc)|=R, choose the cap near vertex
        c=np.array([0,zc]); oc=p-c; bq=np.dot(oc,dvec); cq=np.dot(oc,oc)-R*R
        disc=np.sqrt(bq*bq-cq); ts=[-bq-disc,-bq+disc]; ts=[t for t in ts if t>1e-9]
        # choose point with z nearest the vertex region
        pts=[p+t*dvec for t in ts]; vz=Dv if sign>0 else 0.0
        q=min(pts,key=lambda P: abs(P[1]-vz)); nrm=(q-c)/R
        dvec=dvec-2*np.dot(dvec,nrm)*nrm; p=q
    # propagate to plane z=Dv
    t=(Dv-p[1])/dvec[1]; return (p+t*dvec)[0]
R=2.0
for Dv in (1.0,3.0):
    xs=[trace_exact(Dv,R,1e-3,s) for s in (-0.02,-0.01,0.01,0.02)]
    print('D=',Dv,' ray heights at hole plane for object h=1e-3:',xs)
