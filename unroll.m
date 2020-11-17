tracking_values = zeros(length(starting_solutions),length(starting_solutions(1,:)),n); % For plotting
%(solutions, variables, h)
%dt = 1/n;
dt = 1e-12;
time = zeros(1,n); % For plotting
counter = 1;
tval = 0;
xval = double(starting_solutions(1:1,:));

Hx = matlabFunction(jac_eval)
Ht = matlabFunction(jacobian(h,t))

%xval = rk45(p,xval,tval,dt,bezuit_bound);
temp_p = predictor(xval(:,1),xval(:,2),xval(:,3),tval + dt);
temp_p = reshape(temp_p,1,3);
xval_p = xval + temp_p ; % Euler Predicted Value

Hx0 = Hx(xval_p(:,1),xval_p(:,2),xval_p(:,3),tval + dt)
Ht0 = Ht(xval_p(:,1),xval_p(:,2),xval_p(:,3),tval + dt)
%% WANT: x(dt) ~ x(0) + dx/dt (0) * dt = x(0) + dxdt0 * dt
dxdt0 = Hx0 \ (-Ht0) %% vs x_pred
xval_p_new = xval + dt * dxdt0'

%% test predictor
hval=hEval(xval(:,1),xval(:,2),xval(:,3),0);
hval_p=hEval(xval_p(:,1),xval_p(:,2),xval_p(:,3),dt);
hval_p_new=hEval(xval_p_new(:,1),xval_p_new(:,2),xval_p_new(:,3),dt);

%fprintf('Before Corrector %d %d %d\n',xval(1,1),xval(1,2),xval(1,3));

temp = corrector(xval(:,1),xval(:,2),xval(:,3),tval+dt);
temp = reshape(temp,bezuit_bound,3);
xval_c = xval - temp ;% Corrected Value

corrector_counter = 1;

while corrector_counter < 3
    check = (norm(h(xval_c(:,1),xval_c(:,2),xval_c(:,3),tval + dt)) >= correction_criterion);
    if check
        %xval_c = xval;
        temp = corrector(xval_c(:,1),xval_c(:,2),xval_c(:,3),tval + dt);
        temp = reshape(temp,bezuit_bound,3);
        xval_c = xval_c - temp;
        corrector_counter = corrector_counter + 1;
    elseif ~check && corrector_counter < 3
        xval = xval_c;
    else
        dt = dt/2;
        corrector_counter = counter+counter + 1;
    end
end

%     diagonal_vecs = diag((1./vecnorm(xval')));
%     xval = diagonal_vecs*xval;

tval = tval + dt;
tracking_values(:,:,counter) = xval;
time(counter) = tval;