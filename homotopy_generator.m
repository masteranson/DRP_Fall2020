function [predictor, corrector,jac,h] = homotopy_generator(starting,target,patch, variables)
syms t %time
%variables = [variables, t];
%syms  h(variables) p(variables)

target = [target; patch]; %Appended Patch to target solution
gamma = rand + rand*1i; %Gamma Trick
h_variables = [variables, t];
h(h_variables) =  gamma*starting*(1-t) + target*t; %Homotopy
%h = formula(h);
dht = jacobian(h,t);
jac = jacobian(h,variables);
temp = jac^(-1);
p = jac^(-1)*dht; %Predictor
%p(x1,x3,z,t) = jac^(-1)*dht; %Predictor
j = jac^(-1)*h; %Newton's Corrector f(x)/f'(x)

%Convert to numerical solutions
predictor  = matlabFunction(p);
corrector = matlabFunction(j);
jac_eval = matlabFunction(jac);

end