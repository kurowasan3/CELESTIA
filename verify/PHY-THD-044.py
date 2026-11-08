# Q140-1: 密閉した剛体容器内で加熱する氷（飽和蒸気量による上限と定積での蒸発熱）
import sympy as sp
R,T,M,ps,V,Lv,Lf,ci,cw,m=sp.symbols('R T M p_s V L_v L_f c_i c_w m',positive=True)
ms=ps*V*M/(R*T)
lv=Lv-R*T/M        # 定積での単位質量あたりの蒸発熱
sens=lambda m: m*(ci*10+Lf+cw*100)
Q=lambda m,mv: sens(m)+mv*lv
num={R:8.31,T:373,M:0.018,ps:1.01e5,V:1.0,Lv:2.26e6,Lf:3.3e5,ci:2.1e3,cw:4.2e3}
msv=float(ms.subs(num)); print('(1) m_s =',msv,'kg')
print('(2) L_v - RT/M =',float(lv.subs(num)))
for mv in (0.40,0.80):
    q=float(Q(mv,min(mv,msv)).subs(num)); print('(3) m=',mv,' Q=',q)
p40=float((0.40*R*T/(M*V)).subs(num)); print('(4) p(0.40)=',p40)
# 独立検証: エンタルピーで計算  Q = ΔU = ΔH - Δ(pV),  ΔH(蒸発)=m_v L_v, 終状態の pV = m_v R T / M （初めは真空で pV≈0, 液体の体積は無視）
for mv in (0.40,0.80):
    mvap=min(mv,msv); dH=sens(mv)+mvap*Lv; dpV=mvap*R*T/M
    q2=float((dH-dpV).subs(num)); q1=float(Q(mv,mvap).subs(num)); assert abs(q1-q2)<1e-6
    # 終状態の圧力が飽和蒸気圧を超えないこと
    pf=float((mvap*R*T/(M*V)).subs(num)); assert pf<=1.01e5*(1+1e-12)
    print('  check m=',mv,' Q(エンタルピー)=',q2,' p_f=',pf)
print('OK')
# ---- 新しい (5)：最小容積，Q が V' によらないこと，ピストンの場合との差 p_s V_min
Vp,q_=sp.symbols("V' q",positive=True)
Vmin=m*R*T/(M*ps); print('(5) V_min=',Vmin,' m=0.80 kg で',float(Vmin.subs(num).subs(m,0.80)),'m^3')
# 経路 (i)(ii)(iii) と，独立に「終状態の内部エネルギー」から：U = 液体基準 + 蒸気 m(Lv - RT/M) ， (iii) 等温膨張で ΔU=0
U_vap_rel=m*(Lv-R*T/M)       # 100℃ の液体を基準とした 100℃ の水蒸気の内部エネルギー（理想気体：体積によらない）
Qrigid=m*q_+U_vap_rel
assert sp.diff(Qrigid,Vp)==0
# 別の確かめ：エンタルピー ΔH = m(q+Lv) から Δ(pV)=p_f V' = m R T/M（V' によらない）を引く
pf=m*R*T/(M*Vp); assert sp.simplify(m*(q_+Lv)-pf*Vp-Qrigid)==0
Qpiston=m*(q_+Lv)
assert sp.simplify(Qpiston-Qrigid-ps*Vmin)==0
print("    Q_rigid = m(q+Lv-RT/M)（V' によらない），Q_piston - Q_rigid = p_s V_min を確認")
# 数値：m=0.40 kg（V'=1.0 m^3 ≥ V_min）では (3) の値と一致
qv=float((ci*10+Lf+cw*100).subs(num))
assert abs(float(Qrigid.subs(num).subs({m:0.40,q_:qv}))-float(Q(0.40,0.40).subs(num)))<1e-6
print('OK (5)')
