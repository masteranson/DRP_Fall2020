function xPred=predict(inp, Hx, Ht, xCur, tCur, dt)
    Hx0=subs(Hx, inp, [xCur tCur])
    Ht0=subs(Ht, inp, [xCur tCur])
    cx0t0 = Hx0 \ Ht0 % instead of inverse
    xPred = xCur' + dt * cx0t0
end
