# P40-1 顕微鏡のピント合わせの許容幅（眼を接眼レンズの後焦点に置く）
import sympy as sp, numpy as np
f1,f2,l,Dn,y1,s=sp.symbols('f1 f2 l D y1 s',positive=True)
# 接眼レンズ：中間像は L2 の前方 s(<f2)，虚像は前方 s' = s f2/(f2-s)
sp_=s*f2/(f2-s)
# 眼（L2 の後方 f2）から見た虚像の視角：像の高さ y1*sp_/s，距離 sp_+f2
ang=(y1*sp_/s)/(sp_+f2)
print('visual angle simplified:',sp.simplify(ang))      # -> y1/f2
# 虚像までの距離（眼から）= D となる s
s_D=sp.solve(sp.Eq(sp_+f2,Dn),s)[0]; print('s for D:',sp.simplify(s_D))
# 数値
F1,F2,LL,DD=4.0,25.0,165.0,250.0
def obj(sv):
    b=LL-sv; return b*F1/(b-F1), b
a_inf,b_inf=obj(F2); a_D,b_D=obj(float(s_D.subs({f2:F2,Dn:DD})))
print('a_inf=%.5f a_D=%.5f delta=%.3e mm'%(a_inf,a_D,a_inf-a_D))
m_inf=b_inf/a_inf; m_D=b_D/a_D
print('m1',m_inf,m_D,' M',m_inf*DD/F2,m_D*DD/F2)
print('approx delta f2^2/(D m^2)=',F2**2/DD/m_inf**2)
# 独立検証：行列（光線行列）法で物体位置から最終虚像位置・視角を計算
def system(a):
    T=lambda d: np.array([[1,d],[0,1]]); Lm=lambda f: np.array([[1,0],[-1/f,1]])
    M=T(F2)@Lm(F2)@T(LL)@Lm(F1)@T(a)   # 物体 -> 眼の位置
    # 物体の先端から出る光線 (y=1, u)。眼の位置で y=0 となる u を選ぶ（主光線）-> 視角 = 出射角
    u=-M[0,0]/M[0,1]; out=M@np.array([1.0,u]); ang=out[1]
    # 最終像位置：眼の位置から後ろ向きに像距離 = -B/D 型の計算（物体点からの光束の収束点）
    Mo=Lm(F2)@T(LL)@Lm(F1)@T(a)  # L2 直後
    # 像：y_out independent of u -> 位置 z で M00+..: 求める z: (T(z)Mo)[0,1]=0
    z=-Mo[0,1]/Mo[1,1]
    return ang, z
for a in [a_inf-1e-9, a_D]:
    ang,z=system(a); print('a=%.5f visual angle per unit height=%.4f (expected %.4f)  image z from L2=%.2f -> dist from eye %.2f'%(a,ang,(LL-F2- (0))*0+ (1/F2)*( (LL-(F2 if a==a_inf-1e-9 else float(s_D.subs({f2:F2,Dn:DD}))))/a ),z,F2-z))
