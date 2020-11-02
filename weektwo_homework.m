%% Initial Variables
clear
t0 = 0; %Initial time
n = 5000; %steps
dt = 1/n; %Starting Timestep
h_t = 0.001; %timestep half criterion upper bound
s_t = 1e-7;
max_correction = 3;
correction_criterion = 0.0001;

%% Target System Setup
% syms f1(x,y) f2(x,y) f(x,y) 
% 
% f1(x,y) = x + y; 
% f2(x,y) = x + y + 1; 
% f(x,y) = [f1 ; f2]; %Target System
chemsys

%% Polynomial Homogenization
% syms z %homogenizing variable
% 
% degree  = polynomialDegree(f);
% f(x,y,z) = (z^degree(1)).*f(x/z,y/z);

%% Starting Solution Generation

%Solutions based on roots of unity [e^i2pi/degree; 1]
%Solutions before homogenization (canonical chart z0 = 1)
% starting_solutions =[exp(1i*2*pi)/(degree(1)); exp(1i*2*pi)/(degree(2)); 1];
% 
% %Projective Transformation Euclidean Patch (eliminates scaling)
% syms r(x,y,z) l(x,y,z)
% 
% % Generic Linear Form
% l(x,y,z) =(rand + rand*1i)*x +(rand + rand*1i)*y + (rand + rand*1i)*z; 
% 
% r(x,y,z) = l - 1; %Generic Patch
% starting_solutions = (1/l(1,1,1))*starting_solutions;
% 
% %% Homogenized Total Degree Starting System
% syms fs(x,y) 
% 
% %fs = gamma*(zn^d1 - z(n-1)^d1;z(n-(n-1))^d1-z(0)^d1)
% fs(x,y,z) =[(x^degree(1) - z^degree(1)); (y^degree(2) - z^degree(2)); r];

%% Homotopy
syms  h(x1,x3,z,t) p(x,y,t)

f = [f; r]; %Appended Patch to target solution
gamma = rand + rand*1i; %Gamma Trick
h =  gamma*fs*(1-t) + f*t; %Homotopy

%% Predictor Corrector

dht = jacobian(h,t);
jac = jacobian(h,[x1,x3,z]);
p(x1,x3,z,t) = jac^(-1)*dht; %Predictor
j = jac^(-1)*h; %Newton's Corrector f(x)/f'(x)

%Convert to numerical solutions
p  = matlabFunction(p);
j = matlabFunction(j);

%% Predictor Corrector Loop

xval = double(starting_solutions); %tracking variable vector [n,1]
tval = t0;

tracking_values = zeros(length(starting_solutions),length(starting_solutions(1,:)),n); % For plotting
%(solutions, variables, h)
time = zeros(1,n); % For plotting

xval = double(xval); %Numerical speedup
counter = 1;


while tval < 1

    xval = rk45(p,xval,tval,dt,bezuit_bound);
    %temp_p = p(xval(:,1),xval(:,2),xval(:,3),tval + dt*dt);
    %temp_p = reshape(temp_p,bezuit_bound,3);
    %xval = xval + temp_p ; % Euler Predicted Value
    
    temp = j(xval(:,1),xval(:,2),xval(:,3),tval+dt);
    temp = reshape(temp,bezuit_bound,3);
    xval_c = xval - temp ;% Corrected Value
    
    corrector_counter = 1;
    while corrector_counter <= 3
        check = (abs(norm(xval)-norm(xval_c)) >= correction_criterion);
        if check 
            xval_c = xval;
            temp = j(xval(:,1),xval(:,2),xval(:,3),tval+dt);
            temp = reshape(temp,bezuit_bound,3);
            xval = xval - temp;
            corrector_counter = corrector_counter + 1;
        elseif ~check && corrector_counter == 1
            xval = xval_c;
            corrector_counter = 69;
        end
    end 
    
    
     tval = tval + dt;   
     tracking_values(:,:,counter) = xval;
     time(counter) = tval;
   
%    Adaptive step size
    if mod(counter,200) == 0
        fprintf('Iteration %d\n',counter)
        %dt = adaptive_stepsize_check(j,p,xval,tval,dt,h_t,s_t);
        %fprintf('Residual: %f\n',norm(h(xval(1),xval(2),xval(3),tval)));
    end
    counter = counter + 1;
end


%% Result Visualization



tracked_solutions = length(xval);
fprintf('Tracked %d solutions\n',tracked_solutions);

 for counter = 1:tracked_solutions
     for solutions = 1:3
        plot3(real(tracking_values(counter,solutions,:)),imag(tracking_values(counter,solutions,:)),time,'LineWidth',2);
     end
     
    view(3)

    figure(tracked_solutions)
    title('Homotopy Tracking Path');
    xlabel('\Re');
    ylabel('\Im');
    zlabel('T');
    hold on
    max_val = max(norm(tracking_values(:)));
    [x, y] = meshgrid(-max_val*2:1:max_val*2);
    z = zeros(length(x),length(y));
    surf(x,y,z);
    z = ones(length(x),length(y));
    surf(x,y,z);

     fprintf('Solution %d: ', counter);
     for counterr = 1: length(xval)
        fprintf('%f%+fi ', [real(xval(counterr,counter)), imag(xval(counterr,counter))]);
     end
     fprintf('\n');
     hold off;
 end