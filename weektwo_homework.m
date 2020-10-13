%% Initial Variables
t0 = 0; %Initial time
n = 1000; %steps
dt  = 1/n; %Timestep
h_t = 0.001; %timestep half criterion

%% Target System Setup
syms f1(x,y) f2(x,y) f(x,y) 

f1(x,y) = x + y; 
f2(x,y) = x + y + 1; 
f(x,y) = f1 + f2; %Target System

%% Polynomial Homogenization
syms z

degree  = polynomialDegree(f);
f(x,y,z) = (z^degree).*f(x/z,y/z);

%% Starting System
syms fs(x,y)

%Projective Transformation Euclidean Patch
syms r(x,y,z)

r(x,y,z) = (rand + rand*1i)*x +(rand + rand*1i)*y + (rand + rand*1i)*z - 1;

%fs = gamma*(zn^d1 - z(n-1)^d1;z(n-(n-1))^d1-z(0)^d1)
fs(x,y,z) = [(x^degree - y^degree); (y^degree- z^degree)];

%% Homotopy
syms  h(x,y,t) p(x,y,t)

gamma = 1/r;
h(x,y,z,t) =  gamma*fs(x,y,z)*t + f(x,y,z)*(1-t); %Homotopy

%initial_solution(t) = solve(h,x,y,z); %x(t) solution of h(x,t)
%starting_point = initial_solution(t0);

dht = diff(h,t);
j = jacobian(h(x,y,z,t),[x,y,z,t]);

p(x,y,z,t) = j^(-1)*dht; %Predictor
j = jacobian(h(x,y,z,t0),[x,y,z,t])^(-1)*h(x,y,z,t0); %Corrector


%% Predictor Corrector Loop
%TODO: Convert everything to numerical calculation for speedup
%TODO: Parallelize Code

xval = double(starting_point); %tracking variable
tval = t0;

for counter = 1:n
    
    xval_temp = xval;
    tval_temp = tval;
    xval = double(xval + p(xval,tval)*dt); % Predicted Value
    xval = xval - j(xval,tval)*h(xval,tval);% Corrected Value
    tval = tval+ dt;   

    if mod(counter,200) == 0
        fprintf('Iteration %d\n',counter)
        half_dt = dt/2;
        for counterr = 1:2 %Double half step
            xval_temp = double(xval + p(xval_temp,tval)*half_dt); % Predicted Value
            xval_temp = xval_temp - j(xval_temp,tval)*h(xval_temp,tval);% Corrected Value
            tval_temp = tval+ half_dt;
        end          
        if abs(xval-xval_temp)/xval > h_t %Comparison
            dt = half_dt;
            fprintf('Timestep halved to %f\n',dt);
            xval = xval_temp;
        end   
    end
end

%% Results

tracked_solutions = length(xval);
fprintf('Tracked %d solutions:\n',tracked_solutions);

for counter = 1:tracked_solutions
    fprintf('Solution %d: %f%+fi\n', counter, [real(xval(counter)), imag(xval(counter))]);
end
