%% Initial Variables
t0 = 0; %Initial time
n = 199; %steps
dt  = 1/n; %Timestep
h_t = 0.001; %timestep half criterion

%% Target System Setup
syms f1(x,y) f2(x,y) f(x,y) 

f1(x,y) = x + y; 
f2(x,y) = x + y + 1; 
f(x,y) = [f1 ; f2]; %Target System

%% Polynomial Homogenization
syms z %homogenizing variable

degree  = polynomialDegree(f);
f(x,y,z) = (z^degree(1)).*f(x/z,y/z);

%% Starting Solution Generation

%Solutions based on roots of unity [e^i2pi/degree; 1]
%Solutions before homogenization (canonical chart z0 = 1)
starting_solutions =[exp(1i*2*pi)/(degree(1)); exp(1i*2*pi)/(degree(2)); 1];

%Projective Transformation Euclidean Patch (eliminates scaling)
syms r(x,y,z) l(x,y,z)

% Generic Linear Form
l(x,y,z) =(rand + rand*1i)*x +(rand + rand*1i)*y + (rand + rand*1i)*z; 

r(x,y,z) = l - 1; %Generic Patch
starting_solutions = (1/l(1,1,1))*starting_solutions;

%% Homogenized Total Degree Starting System
syms fs(x,y) 

%fs = gamma*(zn^d1 - z(n-1)^d1;z(n-(n-1))^d1-z(0)^d1)
fs(x,y,z) =[(x^degree(1) - z^degree(1)); (y^degree(2) - z^degree(2)); r];

%% Homotopy
syms  h(x,y,t) p(x,y,t)

f = [f; r]; %Appended Patch to target solution
gamma = rand + rand*1i; %Gamma Trick
h(x,y,z,t) =  gamma*fs(x,y,z)*(1-t) + f(x,y,z)*t; %Homotopy

%% Predictor Corrector

dht = jacobian(h,t);
jac = jacobian(h,[x,y,z]);

p(x,y,z,t) = jac^(-1)*dht; %Predictor
j = jac^(-1)*h; %Newton's Corrector f(x)/f'(x)

%% Predictor Corrector Loop
%TODO: Convert everything to numerical calculation for speedup
%TODO: Parallelize Code

xval = double(starting_solutions); %tracking variable
tval = t0;

for counter = 1:n
    
    xval_temp = xval;
    tval_temp = tval;
    xval = double(xval + p(xval(1),xval(2),xval(3),tval)*dt); % Predicted Value
    xval = xval - j(xval(1),xval(2),xval(3),tval+dt);% Corrected Value
    tval = tval+ dt;   

    %TODO: Visualization R^3 [Re,Ir, T]
    fprintf('Solution step %d: %f %f %f\n', counter, xval(1), xval(2),xval(3));
    fprintf('Residual: %f\n',norm(h(xval(1),xval(2),xval(3),tval)));

    %Adaptive step size
    %TODO: Implement step size bigger criterion
    if mod(counter,200) == 0
        fprintf('Iteration %d\n',counter)
        half_dt = dt/2;
        for counterr = 1:2 %Double half step
            xval_temp = double(xval + p(xval_temp,tval)*half_dt); % Predicted Value
            xval_temp = xval_temp - j(xval_temp,tval)*h(xval_temp,tval);% Corrected Value
            tval_temp = tval+ half_dt;
        end          
        if norm(xval-xval_temp)/norm(xval) > h_t %Comparison
            dt = half_dt;
            fprintf('Timestep halved to %f\n',dt);
            xval = xval_temp;
        end   
    end
end

%% Results
%TODO: Fix it
% 
% tracked_solutions = length(xval);
% fprintf('Tracked %d solutions:\n',tracked_solutions);
% 
% for counter = 1:tracked_solutions
%     fprintf('Solution %d: %f%+fi\n', counter, [real(xval(counter)), imag(xval(counter))]);
% end
