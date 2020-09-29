%% Initial Variables
t0 = 0; %Initial time
n = 1000; %steps
dt  = 1/n; %Timestep
gamma = 1 + 1i; %Randomly selected Gamma


%% System Setup
syms f1(x) f2(x) h(x,t) p(x,t)

f1(x) = x^2 - 1; %Starting System
f2(x) = x^2 + 1; %Target System

h(x,t) = f1(x)*t + f2(x)*(1-t)*gamma; %Homotopy

initial_solution(t) = solve(h, x); %x(t) solution of h(x,t)

starting_point = initial_solution(t0);

dhx = diff(h,x); %dh/dx
dht = diff(h,t); %dh/dt

p(x,t) = (dhx)^(-1)*dht; %Predictor

%follow up
%j = det(jacobian(h(x,t0),[x,t])); %Corrector 
j = (det(diff(h,x)))^(-1); %Corrector 


%% Predictor Corrector Loop

xval = double(starting_point); %tracking variable
tval = t0;

for counter = 1:n
    
xval = double(xval + p(xval,tval)*dt); % Predicted Value
xval = xval - h(xval,tval).*j(xval,tval);% Corrected Value

tval = t0 + dt;
counter = counter + 1;

if mod(counter,200)==0
    fprintf('Iteration %d\n',counter)
end

end
tracked_solutions = length(xval);

fprintf('Tracked %d solutions:\n',tracked_solutions);

for counter = 1:tracked_solutions
    fprintf('Solution %d: %f%+fi\n', counter, [real(xval(counter)), imag(xval(counter))]);
end
