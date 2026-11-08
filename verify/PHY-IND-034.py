# J132-1 検証：結合した2つのコイルの直列接続（和動・差動）と共振周波数による相互インダクタンスの測定
import sympy as sp, numpy as np
mu,S,l,N1,N2,I1,I2,L1,L2,M,C,fp,fm=sp.symbols('mu S l N1 N2 I1 I2 L1 L2 M C f_p f_m',positive=True)
# (1) 鉄心：磁束 Phi = mu S (N1 I1 + N2 I2)/l
Phi=mu*S*(N1*I1+N2*I2)/l
L1e=sp.diff(N1*Phi,I1); L2e=sp.diff(N2*Phi,I2); Me=sp.diff(N1*Phi,I2)
print('(1)',L1e,L2e,Me,' M^2-L1L2=',sp.simplify(Me**2-L1e*L2e))
Lp=sp.factor(L1e+L2e+2*Me); Lm=sp.factor(L1e+L2e-2*Me); print('   L+=',Lp,' L-=',Lm)
# 磁気エネルギー B^2/(2mu)*S l と一致するか（和動，I1=I2=I）
I=sp.symbols('I'); Bf=mu*(N1+N2)*I/l
print('   energy check:',sp.simplify(sp.Rational(1,2)*Lp*I**2-Bf**2/(2*mu)*S*l))
# (2) U>=0 条件：判別式
U=sp.Rational(1,2)*L1*I1**2+sp.Rational(1,2)*L2*I2**2+M*I1*I2
x=sp.symbols('x',real=True)
q=sp.expand(U.subs({I1:x*I2})/I2**2)
print('(2) discriminant:',sp.factor(sp.discriminant(q,x)))
# (3)
Lpv=1/(4*sp.pi**2*fp**2*C); Lmv=1/(4*sp.pi**2*fm**2*C)
Mv=sp.simplify((Lpv-Lmv)/4); print('(3) M=',Mv)
kk=sp.simplify(((Lpv-Lmv)/4)/((Lpv+Lmv)/4)); print('    k (L1=L2)=',kk)
num={C:1e-6,fp:500,fm:1000}
print('    L+=',float(Lpv.subs(num)),' L-=',float(Lmv.subs(num)),' M=',float(Mv.subs(num)),' L1=',float(((Lpv+Lmv)/4).subs(num)),' k=',float(kk.subs(num)))
# (4) コイル2短絡：定常交流で jw(L2 I2 + M I1)=0 -> Leff = L1 - M^2/L2
L1v=(Lpv+Lmv)/4; Leff=sp.simplify(L1v-Mv**2/L1v)
fs=sp.simplify(1/(2*sp.pi*sp.sqrt(Leff*C)))
print('(4) fs=',fs,' num=',float(fs.subs(num)))
# 数値での別法：2x2の回路方程式の行列から共振（虚部0）を探す
Ln=float(L1v.subs(num)); Mn=float(Mv.subs(num)); Cn=1e-6
for f in [float(fs.subs(num))]:
    wv=2*np.pi*f
    # 1次側から見たインピーダンス Z = jwL1 + 1/(jwC) - (jwM)^2/(jwL2)
    Z=1j*wv*Ln+1/(1j*wv*Cn)-(1j*wv*Mn)**2/(1j*wv*Ln)
    print('   Z at fs =',Z)
# 和動の共振の数値確認
wv=2*np.pi*500; print('   Z+ at 500Hz=',1j*wv*(2*Ln+2*Mn)+1/(1j*wv*Cn))
wv=2*np.pi*1000; print('   Z- at 1000Hz=',1j*wv*(2*Ln-2*Mn)+1/(1j*wv*Cn))
