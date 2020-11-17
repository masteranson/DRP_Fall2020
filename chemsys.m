%% Saturation of Silver Chloride

syms f1(x1,x3) f2(x1,x3) 
syms f(x1,x3,z) r(x1,x3,z)
syms fz(x1,x3,z) fs(x1,x3,z)
syms z %homogenizing variable
syms t %time

syms a1 a2 a3 a4 a5 b1 b2 b3
coefficients = [a1, a2, a3, a4, a5, b1, b2, b3];
random_coefficients = [(rand + rand*1i), (rand + rand*1i), (rand + rand*1i),(rand + rand*1i), (rand + rand*1i),(rand + rand*1i),(rand + rand*1i),(rand + rand*1i)];
target_coefficients =  [(rand + rand*1i), (rand + rand*1i), (rand + rand*1i),(rand + rand*1i), (rand + rand*1i),(rand + rand*1i),(rand + rand*1i),(rand + rand*1i)]; % [1, 1, 1, 1, 1, 1, 1, 1];

%% System and homogenization

f1(x1,x3) = a1*x1^4 + a2*x1^3*x3 + a3*x1^3 + a4*x1 +a5; 
f2(x1,x3) = b1*x1*x3^2 + b2*x3^2 + b3; 
f_t  = [f1; f2]; %Target System
f = subs(f_t, coefficients, random_coefficients);

degree  = polynomialDegree(f);
fz(x1,x3,z) = f(x1/z, x3/z);
fz = [z^degree(1); z^degree(2)].* fz; %Homogenized system

%% Homogenized Total Degree Starting System

% Random patch
a = (rand + rand*1i);
b = (rand + rand*1i);
c = (rand + rand*1i);
r = a*x1  + b*x3 + c*z - 1; 

%Starting system
fs =[(x1^degree(1) - z^degree(1)); (x3^degree(2) - z^degree(2)); r];
%% Starting Solutions
%Bezuit bound 12 starting solutions

bezuit_bound = degree(1)*degree(2);
z_starting_solutions = ones(bezuit_bound,1);

x1_starting_solutions = roots([1 zeros(1,degree(1) - 1) -1]);
x3_starting_solutions = roots([1 zeros(1,degree(2) - 1) -1]);
[x1_starting_solutions,x3_starting_solutions] = ndgrid(x1_starting_solutions, x3_starting_solutions);
%x1_starting_solutions(:) = x1_starting_solutions(:) .* 1./(a*x1_starting_solutions(:) + b*x3_starting_solutions(:) + c*1);
%x3_starting_solutions(:) = x3_starting_solutions(:) * 1./(a*x1_starting_solutions(:) + b*x3_starting_solutions(:) + c*1);
%z_starting_solutions = z_starting_solutions * 1./(a*x1_starting_solutions(:) + b*x3_starting_solutions(:) + c*1);

starting_solutions =[x1_starting_solutions(:), x3_starting_solutions(:), z_starting_solutions];
starting_solutions = double(starting_solutions);
starting_solutions = starting_solutions ./ (starting_solutions * [a; b; c]);
