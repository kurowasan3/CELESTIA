import sympy as sp
a,b,f=sp.symbols('a b f',positive=True)
h1,h2,d=sp.Rational(12),sp.Rational(3),sp.Rational(30)
h0=sp.sqrt(h1*h2); mag=sp.sqrt(h1/h2)
sol=sp.solve([sp.Eq(b/a,mag),sp.Eq(b-a,d)],[a,b],dict=True)[0]
D=sol[a]+sol[b]; F=sol[a]*sol[b]/D
print('h0',h0,'D',D,'f',F,'check lens eq',sp.simplify(1/sol[a]+1/sol[b]-1/F))
# the other position: a'=b, b'=a -> magnification a/b
print('m2 =',sol[a]/sol[b],'image2',h0*sol[a]/sol[b])
# merging: quadratic a^2 - D a + D f = 0 discriminant
Dv=sp.symbols('D',positive=True); print('disc zero at',sp.solve(sp.Eq(Dv**2-4*Dv*F,0),Dv))
