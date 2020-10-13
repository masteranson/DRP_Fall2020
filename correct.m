function xCor=correct(inp, Hx, H, xPred, tCur, dt)
    Hx0=subs(Hx, inp, [xPred tCur+dt])
    H0=subs(H, inp, [xPred tCur+dt])
    xCor=xPred' - (Hx0 \ H0)
end
