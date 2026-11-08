# P88-1: 原子の質量と原子核の質量の換算（β-, β+, EC）と完全電離による崩壊の入れ替わり
import sympy as sp
# 記号で条件を変換: 原子質量 Ma(Z)=MN(Z)+Z me (束縛無視)
MN0,MN1,MNm,me=sp.symbols('M_N(Z) M_N(Z+1) M_N(Z-1) m_e')
Z=sp.symbols('Z')
Ma=lambda MN,z: MN+z*me
# β-: MN(Z) > MN(Z+1)+me  → Ma(Z)-Ma(Z+1) = MN0-MN1-me
print('beta- :',sp.simplify(Ma(MN0,Z)-Ma(MN1,Z+1)-(MN0-MN1-me)))
# β+: MN(Z) > MN(Z-1)+me → Ma(Z)-Ma(Z-1)-2me = MN0-MNm-me
print('beta+ :',sp.simplify(Ma(MN0,Z)-Ma(MNm,Z-1)-2*me-(MN0-MNm-me)))
# EC: MN(Z)+me > MN(Z-1) → Ma(Z)-Ma(Z-1)
print('EC    :',sp.simplify(Ma(MN0,Z)-Ma(MNm,Z-1)-(MN0+me-MNm)))
u=931.5; mec2=0.511
Ni,Cu,Zn=63.92797,63.92976,63.92914
print('Cu->Zn beta- Q',(Cu-Zn)*u)
print('Cu->Ni Q_EC',(Cu-Ni)*u,' beta+ Kmax',(Cu-Ni)*u-2*mec2)
Be,Li=7.01693,7.01600
print('Be->Li Q_EC',(Be-Li)*u,' beta+?',(Be-Li)*u-2*mec2)
# 163Dy/Ho (keV)
dMa=-2.83   # Ma(Dy)-Ma(Ho) c^2
dB=12.3     # Btot(Ho)-Btot(Dy)
# 中性原子: Ma = MN + Z me - Btot/c^2 → MN(Dy)-MN(Ho)-me = Ma(Dy)-Ma(Ho) + (Btot(Dy)-Btot(Ho))
q_cont=dMa-dB; print('bare Dy beta- to continuum Q',q_cont)
BK=13.6e-3*67**2; print('BK(Ho66+) Bohr',BK)
print('bound-state Q',q_cont+BK)
# 中性 Ho の EC: Q = Ma(Ho)-Ma(Dy) = 2.83 keV > 0, 完全電離 Ho67+ は EC 不可
# 検算: 別の帳簿（裸の原子核のエネルギー収支を直接）
MN_Dy_minus_MN_Ho = dMa + (0) - (-dB)  # placeholder
# 直接: MN(Dy)=Ma(Dy)-66me+B(Dy); MN(Ho)=Ma(Ho)-67me+B(Ho)
# MN(Dy)-MN(Ho)-me = (Ma(Dy)-Ma(Ho)) + (B(Dy)-B(Ho)) = dMa - dB
print('direct',dMa-dB)
