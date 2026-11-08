# T20251-1 片足だけ摩擦のあるそりの斜面滑走（前後の向きで加速・減速が入れかわる）
import sympy as sp
M,g,al,b,h,mu,mus,a,NP,NQ,f=sp.symbols('M g alpha b h mu mu_s a N_P N_Q f',real=True)
k=mu*h/b
res={}
for case,xP in [('P_down',b),('P_up',-b)]:
    # x: 斜面下向き, y: 斜面垂直上向き. 重心 G 原点, P=(xP,-h), Q=(-xP,-h)
    eqs=[sp.Eq(M*a, M*g*sp.sin(al)-mu*NP), sp.Eq(NP+NQ, M*g*sp.cos(al)),
         sp.Eq(xP*NP-(-h)*(-mu*NP) + (-xP)*NQ, 0)]  # 重心まわり
    s=sp.solve(eqs,[a,NP,NQ],dict=True)[0]
    # 別解：P 点まわり（慣性力 -M a を重心に加える）
    rG=sp.Matrix([-xP,h]); Fg=sp.Matrix([M*g*sp.sin(al)-M*a,-M*g*sp.cos(al)])
    rQ=sp.Matrix([-2*xP,0]); FQ=sp.Matrix([0,NQ])
    tau=(rG[0]*Fg[1]-rG[1]*Fg[0])+(rQ[0]*FQ[1]-rQ[1]*FQ[0])
    s2=sp.solve([eqs[0],eqs[1],sp.Eq(tau,0)],[a,NP,NQ],dict=True)[0]
    print(case,'a=',sp.simplify(s[a]),' agree:',sp.simplify(s[a]-s2[a]),sp.simplify(s[NP]-s2[NP]))
    print('   NP=',sp.simplify(s[NP]),' NQ=',sp.simplify(s[NQ]))
    res[case]=s
aD=g*sp.sin(al)-mu*g*sp.cos(al)*b/(2*b-mu*h)
aU=g*sp.sin(al)-mu*g*sp.cos(al)*b/(2*b+mu*h)
print('formula check',sp.simplify(res['P_down'][a]-aD),sp.simplify(res['P_up'][a]-aU))
# 静止（P下向き）での静止摩擦と垂直抗力
s=sp.solve([sp.Eq(f,M*g*sp.sin(al)),sp.Eq(NP+NQ,M*g*sp.cos(al)),sp.Eq(b*NP-h*f-b*NQ,0)],[f,NP,NQ],dict=True)[0]
print('static NP=',sp.simplify(s[NP]),' NQ=',sp.simplify(s[NQ]))
musmin=sp.simplify(s[f]/s[NP]); print('mu_s min =',musmin)
# 数値例
v={g:9.8,b:0.3,h:0.2,mu:0.6}
import math
tl=float((mu*b/(2*b+mu*h)).subs(v)); tu=float((mu*b/(2*b-mu*h)).subs(v))
alv=math.atan((tl+tu)/2)
print('tan range',tl,tu,' alpha',alv,' aDown',float(aD.subs(v).subs(al,alv)),' aUp',float(aU.subs(v).subs(al,alv)))
print('NQ(down)>0?',float(res['P_down'][NQ].subs(v).subs({al:alv,M:1})))
print('mu_s min at alpha',float(musmin.subs(v).subs(al,alv)),' < mu?')
