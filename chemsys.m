%% Saturation of Silver Chloride

syms f1(x1,x3) f2(x1,x3) 
syms f(x1,x3) r(x1,x3)
syms fz(x1,x3,z) fs(x1,x3,z)
syms z %homogenizing variable

%Cofficients
a1 = 1.069e4;
a2 = 2e14;
a3 = 1.0;
a4 = -1.8e10;
a5 = -1.283e24;
b1 = 2e16;
b2 = 1e14;
b3 = -1;

%% System and homogenization
%TODO: Fix

f1(x1,x3) = a1*x1^4 + a2*x1^3*x3 + a3*x1^3 + a4*x1 +a5; 
f2(x1,x3) = b1*x1*x3^2 + b2*x3^2 + b3; 
f(x1,x3)  = [f1; f2]; %Target System


degree  = polynomialDegree(f);
fz(x1,x3,z) = f(x1/z, x3/z);
fz = [z^degree(1); z^degree(2)].* fz;

%% Homogenized Total Degree Starting System

% Random patch
r = (rand + rand*1i)*x1  + (rand + rand*1i)*x3 + (rand + rand*1i)*z - 1; 

%Starting system
fs =[(x1^degree(1) - z^degree(1)); (x3^degree(2) - z^degree(2)); r];
%% Starting Solutions
%Bezuit bound 12 starting solutions

bezuit_bound = degree(1)*degree(2);
z_starting_solutions = ones(bezuit_bound,1);

x1_starting_solutions = roots([1 zeros(1,degree(1) - 1) -1]);
x3_starting_solutions = roots([1 zeros(1,degree(2) - 1) -1]);
[x1_starting_solutions,x3_starting_solutions] = ndgrid(x1_starting_solutions, x3_starting_solutions);
starting_solutions =[x1_starting_solutions(:), x3_starting_solutions(:), z_starting_solutions];

