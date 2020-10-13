syms x y t
f1(x,y)= x^3-1
f2(x,y) = y^3-1
f=[f1;f2]
g1(x,y)= 1 + x * y^2 + x^(2)  + y^3
g2(x,y) = 2 - x*y^2 + 17*x^(2) * y + 4*x^3
g=[g1;g2]
gamma=exp(2*pi*i*rand())

H(x,y,t) = g(x,y)*t + f(x,y)*(1-t)*gamma; %Homotopy

inp=symvar(H)
Hx = jacobian(H,[x,y])
Ht = jacobian(H,[t])

xCur = [1 1] %% start sol, t=0
tCur = 0
dt=0.0001
xPred=predict(inp, Hx, Ht, xCur, tCur, dt)
xCor = correct(inp, Hx, H, xPred', tCur, dt)
double(subs(H, inp, [xCor' tCur+dt]))
