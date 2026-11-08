# M09-1: person on scale in gondola with movable pulley (3 rope segments support system)
import sympy as sp
m,M,g,al,T,N,h,w = sp.symbols('m M g alpha T N h w', positive=True)
# person: T + N - m g = m al ; whole (person+G): 3T - (m+M) g = (m+M) al
sol = sp.solve([sp.Eq(T+N-m*g, m*al), sp.Eq(3*T-(m+M)*g, (m+M)*al)], [T,N], dict=True)[0]
# alternative: G alone: 2T - M g - N = M al (rope at pulley B: 2 segments; scale pushed down by N)
chk = sp.simplify(2*sol[T]-M*g-sol[N]-M*al)
print('G-alone eq residual:', chk)
T0 = sol[T].subs(al,0); w0 = sp.simplify(sol[N].subs(al,0)/g)
print('(1) T0=',sp.factor(T0),' w0=',sp.factor(w0))
assert sp.simplify(T0-(m+M)*g/3)==0 and sp.simplify(w0-(2*m-M)/3)==0
print('(2) N>=0 <=> 2m-M>=0')
alpha = sp.solve(sp.Eq(sol[N]/g, w), al)[0]
print('(3) alpha=',sp.simplify(alpha))
assert sp.simplify(alpha-(3*w/(2*m-M)-1)*g)==0
# (4) rope geometry: lengths C->B, B->A, A->hand. Ceiling at 0, G top (B and hand level) at -y
y=sp.symbols('y')
Ltot_const = sp.Symbol('Lr')
# rope through hand = change of (CB+BA+A-hand) = 3*dy
seg = 3*y  # each segment length = y + const
L = sp.diff(seg,y)*h
TL = sp.simplify(sol[T]*L)
dE = (m+M)*al*h + (m+M)*g*h   # KE = 1/2 (m+M) v^2 , v^2=2 al h
print('(4) L=',L,' TL=',sp.factor(TL),' dE-TL=',sp.simplify(dE-TL))
assert sp.simplify(dE-TL)==0
# numeric sanity, m=60,M=40 -> w0=80/3 kg
print('numeric w0 (m=60,M=40):', w0.subs({m:60,M:40}))
print('ALL OK')
