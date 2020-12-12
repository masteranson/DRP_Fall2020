# Homotopy Continuation

Homotopy continuation is a classical technique in numerical algebraic geometry that is used to solve for all possible solutions in a non-linear system of algebraic equations. For my final project in the [Directed Reading Program](https://sites.gatech.edu/gtatgt/sample-page/directed-reading-program/), I implemented a very simple type of a homotopy continuation method, called a ***Total Degree Parameter Homotopy*** in MATLAB. My primary reference for this work is based on the textbook written by [Sommese et al.](https://www.worldscientific.com/worldscibooks/10.1142/5763).


To test my numerical solver, I attempted to solve for the concentrations of a saturated silver chloride system in solution:

Reaction:
$AgCl_{(s)}^{2} \leftrightarrow Ag^{+} + Cl^{-}$
$Ag^{+} + OH^- \leftrightarrow AgOH$
$Ag^{+} + Cl^{-} \leftrightarrow AgCl$
$AgCl + Cl^{-} \leftrightarrow AgCl_{2}^{-}$
$AgCl + Ag^{+} \leftrightarrow Ag_{2}Cl^{+}$
$AgCl_{2}^{+} + Cl^{-} \leftrightarrow AgCl_3^{--}$
$AgCl_{2}^{+} + Ag^{+} \leftrightarrow Ag_{3}Cl^{++}$
$H_{2}O_{(l)} \leftrightarrow H_{+} + OH^{-}$

Saturation States:
$[AgCl_{(s)}] = 1$
$[H_{2}O_{(l)}] = 1$

Total Conservation of Species:

$Ag^{+} + AgOH + AgCl +  AgCl_{2}^{-} + 2Ag_{2}Cl^{+} + AgCl_3^{--} + 3Ag_{3}Cl^{++} = S_{Ag}$
$Cl^- + AgCl + 2AgCl_2^{-} + Ag_2Cl^+ 3AgCl_3^{--} + Ag_3Cl^{++} = S_{Cl}$
$OH^- + AgOH = S_{OH}$
$H^{+} = S_{H}$

Equilibrium Constants:

$K_{1} = \frac{[Ag^{+}][Cl^{-}]}{[AgCl_{(s)}]}$
$K_{2} = \frac{[AgOH]}{[Ag^+][OH^-]}$
$K_{3} = \frac{[AgCl]}{[Ag^+][Cl^-]}$
$K_{4} = \frac{[AgCl_2^-]}{[AgCl][Cl^-]}$
$K_{5} = \frac{[Ag_2Cl^-]}{[AgCl][Ag]}$
$K_{6} = \frac{[AgCl_3^{--}]}{[AgCl_2^-][Cl^-]}$
$K_{7} = \frac{[Ag_3Cl^{++}]}{[Ag^+][AgCl_3^{--}]}$
$K_{8} =\frac{[H^+][OH-]}{[H_2O_{(l)}]}$


Using a heuristic method of reduction introduced by [Keith Meintjes](https://www.sciencedirect.com/science/article/abs/pii/0096300387900762), this is the resulting polynomial system:

${\alpha}_1(Ag^+)^4+{\alpha}_2(Ag^+)^3Cl^- +{\alpha}_3(Ag^+)^3 +{\alpha}_4(Ag^+)+{\alpha}_5=0$
${\beta}_1(Ag^+)(Cl^-)^2+{\beta}_2(Cl^-)^2+{\beta}_3$=0

Polynomial Coefficients:
${\alpha}_1= 2K_1K_3K_5K_7$
${\alpha}_2=K_2$
${\alpha}_3=K_1K_3K_5+1$
${\alpha}_4=-K_1(K_1K_3K_4+1)$
${\alpha}_5=-2K_1^3K_3K_4K_6$
${\beta} = \frac{K_2}{K_8}$
${\beta} = \frac{1}{K_8}$
${\beta}=-1$
