# 2004[2] 類題：入れ子の平行板（相互の結合係数の対称性と，スイッチを操作しない電池の充放電）
import sympy as sp
e0,S,D,d,Q,q,V0,V1,V2=sp.symbols('epsilon_0 S D d Q q V_0 V_1 V_2',real=True)
# 面電荷の重ね合わせ：K-N 間は K と N の電荷による Q/(e0 S)，L-M 間はさらに q/(e0 S) が加わる（L,M の外側では打ち消す）
z1=sp.symbols('z1',positive=True)  # M の高さ（結果によらない）
EKL=Q/(e0*S); ELM=(Q+q)/(e0*S); EMN=Q/(e0*S)   # 下向きの大きさ
VKN=sp.simplify(EKL*(D-z1-d)+ELM*d+EMN*z1); VLM=sp.simplify(ELM*d)
print('VKN',VKN,' VLM',VLM)
assert sp.simplify(VKN-(Q*D+q*d)/(e0*S))==0
# (2) KN に V0，LM 短絡
s2=sp.solve([sp.Eq(VKN,V0),sp.Eq(VLM,0)],[Q,q]); print('(2)',s2)
# (3) LM に V0，KN 短絡
s3=sp.solve([sp.Eq(VKN,0),sp.Eq(VLM,V0)],[Q,q]); print('(3)',s3)
assert sp.simplify(s2[q]-s3[Q])==0
# (4)
s4=sp.solve([sp.Eq(VKN,V2),sp.Eq(VLM,V1)],[Q,q]); print('(4)',{k:sp.factor(v) for k,v in s4.items()})
# (5) 初め Q=0, LM=V1
q0=sp.solve(sp.Eq(VLM.subs(Q,0),V1),q)[0]
Qf,qf=s4[Q],s4[q]
dq=sp.simplify(qf-q0); print('dq',dq, ' = -Q?',sp.simplify(dq+Qf))
W1=V1*dq; W2=V2*Qf
U=lambda QQ,qq: (QQ**2*(D-d)+(QQ+qq)**2*d)/(2*e0*S)     # 電場のエネルギー
U2=lambda QQ,qq: sp.Rational(1,2)*(QQ*VKN.subs({Q:QQ,q:qq})+qq*VLM.subs({Q:QQ,q:qq}))
assert sp.simplify(U(Qf,qf)-U2(Qf,qf))==0
H=sp.factor(sp.simplify(W1+W2-(U(Qf,qf)-U(0,q0)))); print('heat',H)
assert sp.simplify(H-e0*S*(V2-V1)**2/(2*(D-d)))==0
# 数値の時間発展（各回路に抵抗 R1,R2）でも熱の総和を確認
import numpy as np
from scipy.integrate import solve_ivp
vals={e0:1.0,S:1.0,D:3.0,d:1.0}
A=np.array([[3.0,1.0],[1.0,1.0]])  # [VKN,VLM]=A[Q,q]/(e0 S)
R1,R2,v1,v2=2.0,0.7,1.0,4.0
def rhs(t,y):
    Qn,qn,H=y
    vkn,vlm=A@np.array([Qn,qn])
    I2=(v2-vkn)/R2; I1=(v1-vlm)/R1
    return [I2,I1,I2**2*R2+I1**2*R1]
sol=solve_ivp(rhs,[0,200],[0.0,1.0,0.0],rtol=1e-10,atol=1e-12)
Hn=sol.y[2,-1]; print('numeric heat',Hn,' theory',(v2-v1)**2/(2*2.0))
assert abs(Hn-(v2-v1)**2/4)<1e-6
print('q final',sol.y[1,-1],' theory',float(qf.subs({**vals,V1:v1,V2:v2})))
