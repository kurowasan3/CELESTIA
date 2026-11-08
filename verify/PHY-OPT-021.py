# T20263-1 望遠レンズ系の合成焦点距離と回折によるボケ・焦点深度
import sympy as sp, numpy as np
a,w,phi=sp.symbols('a w phi',positive=True)
def lens(f): return sp.Matrix([[1,0],[-1/f,1]])
def gap(d): return sp.Matrix([[1,d],[0,1]])
# 光線ベクトル (高さ, 角度)
M2=lens(-3*a)*gap(4*a)*lens(6*a)
# スクリーン位置 s: 平行光線(高さw,角0)が高さ0になる
s=sp.symbols('s')
ray=M2*sp.Matrix([w,0])
ss=sp.solve(sp.Eq(ray[0]+s*ray[1],0),s)[0]; print('screen behind L2 =',ss)
print('exit angle',sp.simplify(ray[1]),' => f_eff =',sp.simplify(-w/ray[1]))
# 斜め平行光（角phi）の像の高さ：主光線を任意に，2本の光線が同じ点に集まることを確認
for h0 in [0,w]:
    r=gap(ss)*M2*sp.Matrix([h0,phi]); print('image height (h0=',h0,') =',sp.simplify(r[0]))
# 単一レンズ L1 だけ
print('L1 only image height =',sp.simplify((gap(6*a)*lens(6*a)*sp.Matrix([0,phi]))[0]))
# 縁の光線の角度：L1の縁(高さR)
R=sp.symbols('R',positive=True)
print('tan beta =',sp.simplify(-(M2*sp.Matrix([R,0]))[1]))
# 光路差の近似の数値確認
r=0.06; beta=np.arctan(0.015/0.18)
for (x,y) in [(0,2e-6),(3e-5,0),(2e-5,1e-6)]:
    W=np.array([r*np.cos(beta),r*np.sin(beta)]); W0=np.array([r,0]); Q=np.array([x,y])
    ex=np.linalg.norm(W-Q)-np.linalg.norm(W0-Q); ap=x*(1-np.cos(beta))-y*np.sin(beta)
    print('dL exact',ex,' approx',ap)
# 数値
lam=5.5e-7; A=0.010; f=18*A; eps=2.0e-4; Dd=1.0e3; Rr=0.015
print('blur delta (R=1.5cm) =',lam/(2*np.sin(np.arctan(Rr/f))),' ~9a lam/R=',9*A*lam/Rr)
print('object-side =',lam/(2*np.sin(np.arctan(Rr/f)))*Dd/f,' ~ lam D/(2R)=',lam*Dd/(2*Rr))
print('depth exact',lam/(2*(1-np.cos(np.arctan(Rr/f)))),' approx',lam*f**2/Rr**2)
rs=f*np.sqrt(lam/eps); print('r_s max =',rs,' object res =',lam*Dd/(2*rs))
b=np.arctan(rs/f); print('check depth at rs exact',lam/(2*(1-np.cos(b))))
