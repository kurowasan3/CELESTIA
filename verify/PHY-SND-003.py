# M69-1: S passes O at t=0 moving toward fixed reflector R (distance D), emits f0 during [0,T]
import sympy as sp, numpy as np
V,v,T,D,f0,te=sp.symbols('V v T D f0 t_e',positive=True)
# arrival time at O (x=0) of the wave emitted at te (S at x=v te)
ta_dir=te+v*te/V
ta_ref=te+(D-v*te)/V+D/V
t1=sp.simplify(ta_dir.subs(te,T)-ta_dir.subs(te,0)); t2=sp.simplify(ta_ref.subs(te,T)-ta_ref.subs(te,0))
print("t1",t1,"  t2",t2)
# frequencies = f0 / (d ta/d te)
fd=sp.simplify(f0/sp.diff(ta_dir,te)); fr=sp.simplify(f0/sp.diff(ta_ref,te))
print("f_dir",fd,"  f_ref",fr)
start_ref=ta_ref.subs(te,0); end_dir=ta_dir.subs(te,T); end_ref=ta_ref.subs(te,T)
print("end_ref-end_dir",sp.simplify(end_ref-end_dir))
overlap=sp.simplify(end_dir-start_ref); print("overlap",overlap)
N=sp.simplify((fr-fd)*overlap); print("beats N",sp.factor(N))
# alt: wave counting -> number of reflected waves arriving in [start_ref,end_dir] minus direct waves in same interval
te_ref_end=sp.solve(sp.Eq(ta_ref,end_dir),te)[0]
te_dir_start=sp.solve(sp.Eq(ta_dir,start_ref),te)[0]
Nref=f0*te_ref_end; Ndir=f0*(T-te_dir_start)
print("wave-count check:",sp.simplify(Nref-Ndir-N))
# numeric sanity
vals={V:340,v:3.4,T:2.0,D:300,f0:600}
print({k:float(e.subs(vals)) for k,e in [('t1',t1),('t2',t2),('overlap',overlap),('N',N)]})
