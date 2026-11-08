# P62-1 浮いた金属板Dを含むコンデンサーの繰り返し操作
from fractions import Fraction as Fr
import sympy as sp
# 単位 C0=1, V=1. C1=3 (上の隙間 d/3), C2=3/2 (下の隙間 2d/3), K2: C=1
C1,C2,C=Fr(3),Fr(3,2),Fr(1)
V=Fr(1)
# 状態：上極板電荷 Qt, D 正味電荷 qD, K2 上極板電荷 qK
Qt=Fr(0); qD=Fr(0); qK=Fr(0)
def energy(phiTop,phiD,phiK):
    return (C1*(phiTop-phiD)**2+C2*phiD**2+C*phiK**2)/2
heat=Fr(0); Wbat=Fr(0)
Vn=[Fr(0)]
phiTop=Fr(0);phiD=Fr(0);phiK=Fr(0)
for n in range(1,40):
    # 操作A: S1 閉 (上=V), S2 開. D 孤立(qD), K2 孤立
    E0=energy(phiTop,phiD,phiK)
    phiD_new=(C1*V+qD)/(C1+C2)
    Qt_new=C1*(V-phiD_new)
    Wb=V*(Qt_new-Qt); Wbat+=Wb
    phiTop,phiD=V,phiD_new; Qt=Qt_new
    heat+=E0+Wb-energy(phiTop,phiD,phiK)
    if n==1: print('after first A: D upper face',-Qt,' lower face',C2*phiD)
    # 操作B: S1 開 (上極板孤立 Qt), S2 閉: D と K2 上極板が同電位
    E0=energy(phiTop,phiD,phiK)
    tot=qD+qK
    # D upper face = -Qt (固定), lower C2*phi, K2 C*phi: -Qt + (C2+C)phi = tot
    phi=(tot+Qt)/(C2+C)
    phiD=phi; phiK=phi; phiTop=phi+Qt/C1
    qD=-Qt+C2*phi; qK=C*phi
    heat+=E0-energy(phiTop,phiD,phiK)
    Vn.append(phi)
    if n<=3: print('n',n,'V_n',phi,' qD',qD,' -C V_n',-C*phi, ' recurrence check',Fr(2,5)+Fr(4,15)*Vn[-2])
Vinf=C1*V/(C1+C2+C)
print('V_40 approx',float(Vn[-1]),' Vinf',Vinf, float(Vinf))
print('closed form n=3:',Fr(6,11)*(1-Fr(4,15)**3),' sim',Vn[3])
print('total heat ~',float(heat),' battery work ~',float(Wbat))
Sig=C1+C2+C
print('theory heat',C1*(C2+C)/Sig/2,float(C1*(C2+C)/Sig/2),' battery',C1*(C2+C)/Sig)
# 一度に両スイッチを閉じた場合: 直列 C1 と (C2+C) を電池で充電 → 熱 = 1/2 Ceq V^2
print('direct heat',C1*(C2+C)/Sig/2)
