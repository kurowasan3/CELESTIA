# P65-1 cube network: nodal analysis
import sympy as sp
R,V,C=sp.symbols('R V C',positive=True)
nodes='ABCDEFGH'
edges=['AB','BC','CD','DA','EF','FG','GH','HE','AE','BF','CG','DH']
def solve(edge_list, fixed, inject=None):
    # fixed: dict node->potential ; returns potentials
    P={n:(fixed[n] if n in fixed else sp.Symbol('p'+n)) for n in nodes}
    eqs=[]
    for n in nodes:
        if n in fixed: continue
        s=0
        for e in edge_list:
            if n in e:
                o=e.replace(n,'')
                s+=(P[n]-P[o])/R
        if inject and n in inject: s-=inject[n]
        eqs.append(s)
    sol=sp.solve(eqs,[P[n] for n in nodes if n not in fixed],dict=True)[0]
    return {n:sp.simplify(sp.sympify(P[n]).subs(sol)) for n in nodes}
def cur_out(P,edge_list,n):
    return sp.simplify(sum((P[n]-P[e.replace(n,'')])/R for e in edge_list if n in e))
# (1) full cube A-G
P=solve(edges,{'A':V,'G':0}); I=cur_out(P,edges,'A'); print('R_AG',sp.simplify(V/I), 'VA-VB',P['A']-P['B'])
# R_AB full
P=solve(edges,{'A':V,'B':0}); print('R_AB full',sp.simplify(V/cur_out(P,edges,'A')))
e2=[e for e in edges if e!='AB']
P=solve(e2,{'A':V,'B':0}); Rp=sp.simplify(V/cur_out(P,e2,'A')); print("R' (AB removed)",Rp)
# open-circuit voltage A-B with AB removed, battery A-G
P=solve(e2,{'A':V,'G':0}); Vth=sp.simplify(P['A']-P['B']); print('Vth',Vth, 'I_batt final',cur_out(P,e2,'A'))
# Thevenin R seen from AB with A,G shorted (both 0), inject at B
P=solve(e2,{'A':0,'G':0},inject={'B':1}); Rth=sp.simplify(P['B']); print('Rth (battery shorted)',Rth)
# initial battery current (cap = short, i.e. full cube) 
P=solve(edges,{'A':V,'G':0}); print('I0',cur_out(P,edges,'A'))
# check Thevenin consistency: V_AB with edge = Vth*R/(R+Rth)
print('check', sp.simplify(Vth*R/(R+Rth) - sp.Rational(2,5)*V))
print('Q',C*Vth,'tau_ch',C*Rth,'tau_dis',C*Rp,'heat',sp.simplify((C*Vth)**2/(2*C)))

# ---- final problem version: edge AB = resistor R in series with capacitor C
print('--- final version')
print('(1) I_batt0 =', sp.Rational(6,5), 'V/R ; I_AB0 = (VA-VB)/R =', sp.Rational(2,5),'V/R')
print('(2) Q =',Vth*C,' I_batt_inf =',sp.Rational(10,11),'V/R')
P=solve(e2,{'A':V,'G':0}); print('    zero-current edges:',[e for e in e2 if sp.simplify(P[e[0]]-P[e[1]])==0])
i_open=sp.simplify(Vth/(R+Rp)); heatAB=sp.simplify((C*Vth)**2/(2*C)*R/(R+Rp))
print('(3) i just after opening =',i_open,' heat in AB resistor =',heatAB)
# numerical transient check: MNA with cap state, R=C=V=1
import numpy as np
from scipy.integrate import solve_ivp
idx={n:i for i,n in enumerate(nodes)}
def node_pots(q,battery):
    # unknown potentials of nodes; AB branch: R in series with cap(voltage q/C): current from A to B = (pA-pB-q)/R
    Rn,Cn,Vn=1.0,1.0,1.0
    N=8; Gm=np.zeros((N,N)); b=np.zeros(N)
    for e in e2:
        i,j=idx[e[0]],idx[e[1]]
        Gm[i,i]+=1;Gm[j,j]+=1;Gm[i,j]-=1;Gm[j,i]-=1
    i,j=idx['A'],idx['B']
    Gm[i,i]+=1;Gm[j,j]+=1;Gm[i,j]-=1;Gm[j,i]-=1
    b[i]+=q; b[j]-=q   # source term from capacitor voltage
    if battery:
        for n,val in (('A',Vn),('G',0.0)):
            k=idx[n]; Gm[k,:]=0; Gm[k,k]=1; b[k]=val
    else:
        k=idx['G']; Gm[k,:]=0; Gm[k,k]=1; b[k]=0
    p=np.linalg.solve(Gm,b)
    return p
def rhs(t,y,battery):
    q=y[0]; p=node_pots(q,battery)
    iAB=p[idx['A']]-p[idx['B']]-q
    return [iAB, iAB**2]
s1=solve_ivp(rhs,[0,40],[0,0],args=(True,),rtol=1e-10,atol=1e-12)
q_inf=s1.y[0,-1]; print('num Q',q_inf,'vs',8/11)
p0=node_pots(0,True); print('num I_AB0',p0[0]-p0[1],'vs 0.4')
s2=solve_ivp(rhs,[0,80],[q_inf,0],args=(False,),rtol=1e-10,atol=1e-12)
p=node_pots(q_inf,False); print('num i open',-(p[0]-p[1]-q_inf),'vs',10/33)
print('num heat AB',s2.y[1,-1],'vs',40/363)
