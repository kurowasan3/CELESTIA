# T20102-1: rod as cross arm of a bridge; Thevenin from two terminal speeds
import sympy as sp
E,R,Rx,eps,Bl,VX,VM=sp.symbols('E R R_x epsilon Bl V_X V_M',real=True)
Rrod=R
def solve_nodes(Ev,Rv,Rxv,epsv):
    # node L at E, Y at 0. Arms: L-X R, X-Y 2R, L-M Rx, M-Y R. Rod X->M: (VX-VM-eps)/Rrod
    Irod=(VX-VM-epsv)/Rv
    eqs=[(Ev-VX)/Rv - VX/(2*Rv) - Irod, (Ev-VM)/Rxv - VM/Rv + Irod]
    s=sp.solve(eqs,[VX,VM],dict=True)[0]
    Ir=sp.simplify(Irod.subs(s))
    Ibat=sp.simplify(((Ev-s[VX])/Rv + (Ev-s[VM])/Rxv))  # current leaving + terminal
    return s,Ir,Ibat
s,Ir,Ib=solve_nodes(E,R,Rx,eps)
Eth=sp.simplify(sp.solve(sp.Eq(Ir,0),eps)[0]); print('E_th=',sp.factor(Eth))
Rtot=sp.simplify(-1/sp.diff(Ir,eps)); print('R_rod+R_th=',sp.factor(Rtot),' R_th=',sp.factor(Rtot-R))
# numeric: R=3, Rx=6, E=1.5, Bl=0.1
num={E:1.5,R:3,Rx:6}
print('E_th',Eth.subs(num),'Rtot',Rtot.subs(num))
v0=Eth.subs(num)/0.1; F=0.02; I=F/0.1; v1=(Eth.subs(num)+I*Rtot.subs(num))/0.1
print('v0',v0,'v1',v1)
# inverse: from v0=5, v1=19 recover Rx,E
Rxs=sp.symbols('Rxs',positive=True); Es=sp.symbols('Es',positive=True)
RtotX=Rtot.subs({R:3,Rx:Rxs}); EthX=Eth.subs({R:3,Rx:Rxs,E:Es})
rx=sp.solve(sp.Eq(RtotX,(19-5)*0.1/0.2),Rxs); print('Rx from data',rx)
es=sp.solve(sp.Eq(EthX.subs(Rxs,rx[0]),0.5),Es); print('E from data',es)
# terminal with F: eps=Bl*v1
s1,Ir1,Ib1=solve_nodes(1.5,3,6,0.1*v1)
print('rod current (X->M)',Ir1,'battery current',Ib1, 'VX,VM',s1)
P_F=F*v1
J=( (1.5-s1[VX])**2/3 + s1[VX]**2/6 + (1.5-s1[VM])**2/6 + s1[VM]**2/3 + Ir1**2*3 )
print('P_F',P_F,'Joule',J,'battery power out',1.5*Ib1,' check P_F + Pbat = J:',P_F+1.5*Ib1, J)
# speed for zero battery current
vz=sp.solve(sp.Eq(solve_nodes(1.5,3,6,eps)[2],0),eps); print('eps for Ibat=0',vz,[x/0.1 for x in vz])
