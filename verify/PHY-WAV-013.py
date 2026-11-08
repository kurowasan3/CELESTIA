# K20023-1 verification: radio reflection from moving sound-wave layers (RASS)
import sympy as sp
c,vs,fs,f,eps=sp.symbols('c v_s f_s f epsilon',positive=True)
Lam=vs/fs
# moving-frame: f' = f(1-vs/c); lambda'=c/f'; Bragg 2Lam = lambda' (n=1)
f0_exact=sp.solve(sp.Eq(2*Lam, c/(f*(1-vs/c))),f)[0]
print('f0 exact (layer frame):',sp.simplify(f0_exact))
print('f0 - (c fs/(2vs)+fs/2) =',sp.simplify(f0_exact-(c*fs/(2*vs)+fs/2)),' relative order:',sp.simplify((f0_exact-(c*fs/(2*vs)+fs/2))/f0_exact))
fr=f*(1-vs/c)/(1+vs/c)
df=sp.series((f-fr).subs(vs,eps*c),eps,0,2).removeO().subs(eps,vs/c)
print('Delta f first order:',sp.simplify(df))
print('Delta f at f0 (first order expr):',sp.expand(df.subs(f,c*fs/(2*vs)+fs/2)))
# photon-phonon
fr2=sp.symbols('f_r')
sol=sp.solve([sp.Eq(f,fr2+fs),sp.Eq(f/c,-fr2/c+fs/vs)],[f,fr2],dict=True)[0]
print('photon-phonon f =',sp.simplify(sol[f]),' f_r =',sp.simplify(sol[fr2]))
# numbers
cn=3.00e8; fsn=1.000e3; f0n=461.5e6
vsn=cn*fsn/(2*f0n-fsn); print('v_s =',vsn)
T=273*(vsn/331.5)**2; print('T =',T,'K',T-273,'C')
t1=3.000; t2=6.80e-6; h=cn*t2/2; print('h=',h,' avg v =',h/(t1+t2/2))
