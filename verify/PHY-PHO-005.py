# P81-1: 接触電位差を含む光電管．エネルギー準位モデルから独立に導出し，数値データを照合
import sympy as sp
e,h,nu,V,WA,WC=sp.symbols('e h nu V W_A W_C',positive=True)
# 基準準位: C の基準準位 0, A の基準準位 -eV．表面直外の静止電子のエネルギー
EvC=WC; EvA=-e*V+WA
# 電位 φ = -(エネルギー)/e → φA-φC
dphi=sp.simplify(-(EvA-EvC)/e); print('phiA-phiC =',dphi)
# 最大エネルギー電子: エネルギー hnu (基準0からの) が EvA 以上なら到達
Vstop=sp.solve(sp.Eq(h*nu,EvA),V)[0]; print('V at cutoff =',sp.simplify(Vstop))  # = (WA-h nu)/e → V0=(h nu-WA)/e
Vsat=sp.solve(sp.Eq(EvA,EvC),V)[0]; print('V sat =',Vsat)
print('width =',sp.simplify(Vsat-Vstop))
# 数値: ν1=1.20e15: 開始 -0.45 V, ν2=0.90e15: 開始 +0.79 V, 飽和 2.20 V
ee=1.6e-19
h_num=ee*(0.79-(-0.45))/(1.20e15-0.90e15); print('h=',h_num)
WA_num=(h_num*1.20e15-ee*0.45)/ee; print('WA[eV]=',WA_num)
WC_num=WA_num-2.20; print('WC[eV]=',WC_num)
print('nu_min=',WC_num*ee/h_num)
# 別法: ν2 側から WA
print('WA from nu2 =',(h_num*0.90e15+ee*0.79)/ee)
# 素朴な解釈の限界振動数
print('naive threshold',WA_num*ee/h_num)
# モンテカルロ: 平行平板, 等方放出, 運動エネルギー一様 0..Kmax の電子 → 到達条件 法線成分 ≥ e(φC-φA)
import numpy as np
rng=np.random.default_rng(0); hnu=6.6e-34*1.2e15/ee; Kmax=hnu-2.3; WAv=4.5; WCv=2.3
K=rng.uniform(0,Kmax,200000); ct=np.sqrt(rng.uniform(0,1,200000)) # 余弦則で方向
Kn=K*ct**2
for Vv in [-0.46,-0.44,0.0,2.19,2.21]:
    barrier=-(Vv+(WCv-WAv))  # e(φC-φA) in eV
    print(Vv, np.mean(Kn>=barrier))
