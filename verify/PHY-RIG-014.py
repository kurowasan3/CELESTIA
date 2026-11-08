# P14 類題：最低点のペダルを後ろ向きに押す自転車（仮想仕事と力のモーメントの2通り）
import sympy as sp
F,c,R,r1,r2,M,dx=sp.symbols('F c R r1 r2 M dx',positive=True)
# kinematics: bike forward dx -> wheel angle dx/R (clockwise, forward rolling), crank angle = (dx/R)*(r2/r1) same sense
dphi=dx/R*r2/r1
pedal=dx-c*dphi   # bottom pedal ground displacement (forward +)
print('pedal disp =',sp.simplify(pedal))
# work by F (backward): W = -F*pedal ; effective force on bike (forward +) = -F*(1 - c r2/(R r1))
Feff=sp.simplify(-F*pedal/dx); print('Feff=',Feff)
# force method: crank torque balance: F*c = T*r1 (chain tension T on upper run pulling...), rear: T*r2 = f*R
T=F*c/r1; f=T*r2/R
net=-F+f
print('net external horizontal =',sp.simplify(net)); assert sp.simplify(net-Feff)==0
vals={F:30,c:0.17,R:0.34,r1:2,r2:1,M:15}
print('a =',float((Feff/M).subs(vals)),' f=',float(f.subs(vals)))
print('balance gear ratio r1/r2 = c/R =',0.17/0.34)
# 符号の確認（力のモーメント）：ペダル(0,-c)に力(-F,0)→ z成分 = x*Fy - y*Fx
x,y=0,-c; Fx,Fy=-F,0
print('torque z on crank =',x*Fy-y*Fx,' (negative = clockwise = forward pedaling when bike faces +x)')
