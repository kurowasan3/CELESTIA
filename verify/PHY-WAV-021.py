# PHY-WAV-021 検証：動く振動波源がつくる深水波（分散性）と，ストロボによる見かけの運動
import sympy as sp, numpy as np
g,f,u,lam=sp.symbols('g f u lambda',positive=True)
v=sp.sqrt(g*lam/(2*sp.pi))
# (2) 前方: f*lam = v - u   -> s=sqrt(lam): f s^2 - c s + u = 0
s=sp.symbols('s',positive=True); c=sp.sqrt(g/(2*sp.pi))
roots=sp.solve(sp.Eq(f*s**2-c*s+u,0),s); print(roots)
disc=sp.simplify(c**2-4*f*u); ucrit=sp.solve(sp.Eq(disc,0),u)[0]; print('ucrit=',ucrit)
lam_plus=sp.simplify((g/(8*sp.pi*f**2))*(1+sp.sqrt(1-8*sp.pi*f*u/g))**2)
lam_minus=sp.simplify((g/(8*sp.pi*f**2))*(1-sp.sqrt(1-8*sp.pi*f*u/g))**2)
for L in (lam_plus,lam_minus):
    print('residual',sp.simplify((f*L-v.subs(lam,L)+u).subs({g:9.8,f:3,u:0.2})))
# (3) group velocity v/2 > u selects lam_plus
import random
for _ in range(200):
    gn=9.8; fn=random.uniform(0.5,5); un=random.uniform(0,1)*gn/(8*np.pi*fn)*0.999
    lp=float(lam_plus.subs({g:gn,f:fn,u:un})); lm=float(lam_minus.subs({g:gn,f:fn,u:un}))
    vg=lambda L: 0.5*np.sqrt(gn*L/(2*np.pi))
    assert vg(lp)>un and vg(lm)<un
# at critical: v_g = u, lam = g/(8 pi f^2), lab freq 2f
lc=lam_plus.subs(u,ucrit); print('lam crit',sp.simplify(lc))
print('vg-u at crit',sp.simplify(v.subs(lam,lc)/2-ucrit))
print('f_lab crit',sp.simplify(v.subs(lam,lc)/lc))
# (4) strobe at f_s=f : apparent displacement per flash = u/f (mod lambda) for any pattern
#   front: lab freq f_F = v/lam = f + u/lam ; back: f_B = f - u/lam
# 見かけの静止条件 u/f = n*lam_F （1回の閃光間の移動が波長の整数倍）
#   前方波 lam_F (=lam_plus) は v_g>u より f*lam_F = v-u > u つまり lam_F > u/f → n>=1 は不可（臨界 n=1 のみ等号）
for _ in range(300):
    gn=9.8; fn=random.uniform(0.5,5); un=random.uniform(0.01,0.999)*gn/(8*np.pi*fn)
    lp=float(lam_plus.subs({g:gn,f:fn,u:un}))
    assert lp>un/fn
print('lam_F > u/f confirmed; at critical lam_F*f/u =',sp.simplify(lc*f/ucrit))
# 数値シミュレーション: 山の位置を追跡し，f_s=f の閃光ごとの見かけの移動を確認
gn,fn,un=9.8,2.0,0.15
lF=float(lam_plus.subs({g:gn,f:fn,u:un})); vF=np.sqrt(gn*lF/(2*np.pi))
crest=lambda T,m: un*(m/fn)+vF*(T-m/fn)
T0=20.0; T1=T0+1/fn
p0=np.array([crest(T0,m) for m in range(0,60)]); p1=np.array([crest(T1,m) for m in range(0,62)])
a0=np.sort(p0[p0>un*T0]); a1=np.sort(p1[p1>un*T1])
print('spacing',np.diff(a0)[:3],'lamF',lF)
# a1 の山のうち a0 の山から最も近い前方位置への移動量
d=[(x-a0[a0<=x].max()) for x in a1[3:8]]
print('apparent shift per flash',d,' u/f=',un/fn)
assert np.allclose(d,un/fn)
fF=vF/lF; assert abs(fF-(fn+un/lF))<1e-12   # 静止した観測者の見る振動数 f_F=f+u/λ_F
# ストロボ f_s=f_F なら前方の縞は水面に対して静止
p2=np.array([crest(T0+1/fF,m) for m in range(0,62)]); a2=np.sort(p2[p2>un*(T0+1/fF)])
print('shift with fs=fF (mod lam):',[round(((x-a0[a0<=x].max())%lF),9) for x in a2[3:8]])
# 後方の波：apparent velocity with f_s=f equals u (moves with source)
lamB=sp.symbols('lamB',positive=True)
fB=sp.sqrt(g*lamB/(2*sp.pi))/lamB
# f lamB = v + u -> fB = f - u/lamB, apparent speed lamB*(f-fB) = u
expr=sp.simplify(lamB*(f-(f-u/lamB))); print('apparent back speed',expr)
# (5) 新設問：lambda_F = (lam0/4)(1+sqrt(1-u/uc))^2，範囲 [lam0/4, lam0]，f_F in [f, 2f]
lam0=g/(2*sp.pi*f**2)
lamF_alt=(lam0/4)*(1+sp.sqrt(1-u/ucrit))**2
assert sp.simplify(lamF_alt-lam_plus)==0
assert sp.simplify(lam_plus.subs(u,0)-lam0)==0 and sp.simplify(lam_plus.subs(u,ucrit)-lam0/4)==0
dl=sp.diff(lam_plus,u); 
fF=f+u/lam_plus
for _ in range(300):
    fn_=random.uniform(0.5,5); uc_=9.8/(8*np.pi*fn_); un_=random.uniform(1e-4,0.9999)*uc_
    sub={g:9.8,f:fn_,u:un_}
    assert float(dl.subs(sub))<0
    ff=float(fF.subs(sub)); assert fn_<ff<2*fn_
    assert abs(ff-float((v.subs(lam,lam_plus)/lam_plus).subs(sub)))<1e-9
print('f_F at u=uc:',sp.simplify(fF.subs(u,ucrit)),' at u=0:',sp.simplify(fF.subs(u,0)))
# 非分散の波との比較：f*lam=c0-u は u->c0 で 0。深水波の前方波長の下限 lam0/4 は v_g=u の点
print('min lambda_F / lambda_0 =',sp.simplify(lam_plus.subs(u,ucrit)/lam0))
print('OK')
