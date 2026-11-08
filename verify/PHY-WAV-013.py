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
# (6) 新設問：下降する音波（反射面が近づく）
import numpy as np
f0d_exact=sp.solve(sp.Eq(2*Lam, c/(f*(1+vs/c))),f)[0]
print('f0 descending exact:',sp.simplify(f0d_exact))
rel=sp.simplify((f0d_exact-(c*fs/(2*vs)-fs/2))/f0d_exact); print('  relative diff:',rel)
assert sp.simplify(rel-vs**2/c**2)==0   # 差は (v_s/c)^2 の相対誤差＝2次
frd=f*(1+vs/c)/(1-vs/c)
dfd=sp.series((frd-f).subs(vs,eps*c),eps,0,2).removeO().subs(eps,vs/c)
print('f_r - f first order:',sp.simplify(dfd),' at f0:',sp.expand(dfd.subs(f,c*fs/(2*vs)-fs/2)))
assert sp.simplify(sp.expand(dfd.subs(f,c*fs/(2*vs)-fs/2))-fs+fs*vs/c)==0   # = fs - fs*vs/c（第2項は f0 自体の1次補正との積で2次の量）
sold=sp.solve([sp.Eq(f+fs,fr2),sp.Eq(f/c-fs/vs,-fr2/c)],[f,fr2],dict=True)[0]
print('photon absorbs phonon: f =',sp.simplify(sold[f]),' f_r =',sp.simplify(sold[fr2]))
assert sp.simplify(sold[f]-(c*fs/(2*vs)-fs/2))==0 and sp.simplify(sold[fr2]-sold[f]-fs)==0
# 数値：多数の反射面からの反射波を直接足し合わせ，強め合いの振動数を求める（下降・上昇の両方）
cn_,vsn_,fsn_=3.0e8,330.0,1.0e3
def refl_power(fr,sign,N=400,T=1.0e-5):
    # 反射面 m の高さ z_m(t)=z0+m*Lam+sign*vs*t（sign=+1:上昇，-1:下降）。
    # 地上で時刻 T に受ける反射波：反射時刻 t_r は z_m(t_r)=c(T-t_r)，発射時刻 t_e=2 t_r - T
    Lam_=vsn_/fsn_; z0=1000.0
    m=np.arange(N)
    tr=(cn_*T-z0-m*Lam_)/(cn_+sign*vsn_)
    te=2*tr-T
    return abs(np.sum(np.exp(2j*np.pi*fr*te)))
for sign,f0th in [(+1,cn_*fsn_/(2*vsn_)+fsn_/2),(-1,cn_*fsn_/(2*vsn_)-fsn_/2)]:
    fsc=np.linspace(f0th-3*fsn_,f0th+3*fsn_,6001)
    pw=[refl_power(x,sign) for x in fsc]
    fbest=fsc[int(np.argmax(pw))]
    print('sign',sign,' f0 numeric',fbest,' theory',f0th,' diff',fbest-f0th)
    assert abs(fbest-f0th)<3.0
print('PHY-WAV-013 (6) OK')
