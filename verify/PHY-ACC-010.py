# J139-1 検証: 極板間隔を周期的に切り替えるLC回路（パラメトリック増幅）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
L,C,Q0=sp.symbols('L C Q0',positive=True)
# (1) 仕事 = エネルギー増加 = 力×距離
W1=Q0**2/(2*(C/2))-Q0**2/(2*C)
eps0,S,d=sp.symbols('epsilon0 S d',positive=True)
F=Q0**2/(2*eps0*S); Wforce=(F*d).subs(eps0,C*d/S)
print('W1',sp.simplify(W1),sp.simplify(Wforce))
# (2) t1, I1
wp=1/sp.sqrt(L*C/2); t1=sp.pi/(2*wp); I1=Q0*wp
print('t1',sp.simplify(t1),'I1',sp.simplify(I1))
# (3) Q1
Q1=I1*sp.sqrt(L*C); print('Q1',sp.simplify(Q1))
# 数値: L=1,C=1,Q0=1
Lv,Cv=1.0,1.0
def run(n_cycles,mode):
    # mode 'pump': |Q|max で C->C/2, Q=0 で C/2->C;  'reverse': 逆
    t=0.0; y=[1.0,0.0]; Ccur = Cv/2 if mode=='pump' else Cv
    Emax=[]; times=[]
    for k in range(2*n_cycles):
        target_zero_Q = (abs(y[1])<1e-9)  # 今I=0なら次はQ=0を待つ
        f=lambda t,y,Cc=Ccur:[y[1],-y[0]/(Lv*Cc)]
        if target_zero_Q:
            ev=lambda t,y:y[0]
        else:
            ev=lambda t,y:y[1]
        ev.terminal=True; ev.direction=0
        sol=solve_ivp(f,[t,t+100],y,events=ev,rtol=1e-11,atol=1e-13,max_step=0.01)
        t=sol.t_events[0][0]; y=list(sol.y_events[0][0])
        if target_zero_Q: y[0]=0.0
        else: y[1]=0.0
        # 切り替え
        if mode=='pump':
            Ccur = Cv if target_zero_Q else Cv/2
        else:
            Ccur = Cv/2 if target_zero_Q else Cv
        if not target_zero_Q:
            times.append(t); Emax.append(y[0]**2/(2*Cv))
    return times,Emax
tp,Ep=run(5,'pump')
tau=(np.pi/2)*np.sqrt(Lv*Cv)*(1+1/np.sqrt(2))
print('pump times',np.round(tp,6),'pred',[round(n*tau,6) for n in range(1,6)])
print('pump E/E0',np.round(np.array(Ep)/0.5,6))
tr,Er=run(4,'reverse')
print('reverse E/E0',np.round(np.array(Er)/0.5,6))
# 1000倍
print('n for 2^n>=1000:',next(n for n in range(30) if 2**n>=1000))
