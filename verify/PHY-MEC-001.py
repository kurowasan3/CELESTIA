# T20101-1: elastic collision at loop bottom, A not removed; h=4R. ranges of r=m2/m1
import sympy as sp, numpy as np
r,R,g=sp.symbols('r R g',positive=True)
h=4*R; v2=2*g*h
vB=2/(1+r); vA=(1-r)/(1+r)   # factors of v
# B loops: vB^2 v^2 >= 5 g R
rB=sp.solve(sp.Eq(vB**2*v2,5*g*R),r); print('B boundary',rB)
# A safe: hA=vA^2 h <=R or >=5R/2
print(sp.solve(sp.Eq(vA**2*h,R),r), [sp.radsimp(x) for x in sp.solve(sp.Eq(vA**2*h,sp.Rational(5,2)*R),r)])
print('4sqrt10/5-1=',float(4*sp.sqrt(10)/5-1),' (13-4sqrt10)/3=',float((13-4*sp.sqrt(10))/3))
# brute force scan
def ok(rv):
    b=(2/(1+rv))**2*4; a=((1-rv)/(1+rv))**2*4  # heights in units of R
    return b>=2.5 and (a<=1 or a>=2.5)
rs=np.linspace(0.001,3,300000); good=rs[[ok(x) for x in rs]]
# find intervals
br=np.where(np.diff(good)>1e-4)[0]; print('intervals:',good[0],good[br],good[br+1],good[-1])
# part 4: r=1/2: B height on friction slope
rv=sp.Rational(1,2); vB2=(2/(1+rv))**2*2*g*h
h4=R+(vB2-2*g*R)/(4*g); print('h4',sp.simplify(h4/R),'A height',sp.simplify(((1-rv)/(1+rv))**2*h/R))
# numeric check of slope deceleration via integration: along slope a=-g(sin+mu cos)=-2g sin
th=0.5; s=float((vB2-2*g*R).subs({g:9.8,R:1}))/(2*2*9.8*np.sin(th)); print('rise',s*np.sin(th),'expected',float(((vB2-2*g*R)/(4*g)).subs({g:9.8,R:1})))
# also verify normal force inside loop for B at top >=0 for r=1/2
print('B top v^2/(gR)=',sp.simplify((vB2-4*g*R)/(g*R)))
