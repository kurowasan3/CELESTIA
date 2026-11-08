# O20083-1 : すきま容積のあるポンプでタンクを等温で加圧する
import sympy as sp, numpy as np
P0,V2,VT,a,n=sp.symbols('P0 V2 V_T a n',positive=True)
V1=a*V2
p=sp.Function('p')
# 漸化式 p_n (VT+V2) = p_{n-1} VT + P0 V1
r=VT/(VT+V2); ps=a*P0
pn=ps+(P0-ps)*r**n
print('check recurrence:',sp.simplify(pn.subs(n,n+1)*(VT+V2)-(pn*VT+P0*V1)))
# 逐次シミュレーション（モル数で追跡, 等温, RT=1）
def simulate(av,VTv,strokes,P0v=1.0,V2v=1.0):
    V1v=av*V2v; ptank=P0v; out=[]
    for i in range(strokes):
        ncyl=P0v*V1v            # 吸気後のシリンダー
        # 圧縮: 体積 V1->V2 を細かく
        Vs=np.linspace(V1v,V2v,200001); opened=False; ntank=ptank*VTv
        for Vc in Vs:
            pc=ncyl/Vc
            if not opened and pc>=ptank: opened=True; Vopen=Vc
            if opened:
                ptot=(ncyl+ntank)/(Vc+VTv); ncyl=ptot*Vc; ntank=ptot*VTv
        ptank=ntank/VTv; pdead=ncyl/V2v
        Vin=V2v*pdead/P0v       # 戻り行程で吸気弁が開く体積
        out.append((ptank,Vopen if opened else None,Vin))
    return out
res=simulate(5,9,20)
for i,(pt,vo,vi) in enumerate(res[:4]+res[13:15]):
    pass
for i in [0,1,2,12,13,14]:
    print('stroke',i+1,' tank p',round(res[i][0],6),' formula',float(pn.subs({a:5,P0:1,V2:1,VT:9,n:i+1})),' outlet opens V',res[i][1],' inlet opens V',round(res[i][2],5))
nmin=sp.log(4)/sp.log(sp.Rational(10,9)); print('n>=',float(nmin))
print('first exceed 4P0 at stroke', next(i+1 for i,x in enumerate(res) if x[0]>=4))
# 出口弁が開く体積 (n回目の圧縮): V1 P0/p_{n-1}
print('outlet open vol formula stroke2:',float((5*1*1/pn).subs({a:5,P0:1,V2:1,VT:9,n:1})))
