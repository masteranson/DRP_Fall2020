function [xval, tracking_values, time] = linear_homotopy(predictor, corrector, starting_solutions, n, bezuit_bound, correction_criterion, h_t, jac_eval,h)

    tracking_values = zeros(length(starting_solutions),length(starting_solutions(1,:)),n); % For plotting
    %(solutions, variables, h)
    %dt = 1/n;
    dt = 1e-2;
    time = zeros(1,n); % For plotting
    counter = 1;
    tval = 0;
    xval = double(starting_solutions);
    successes = 0; % consecutive
    while tval < .98
        fprintf('at time: %f\n',tval);

        %xval = rk45(p,xval,tval,dt,bezuit_bound);
        temp_p = predictor(xval(:,1),xval(:,2),xval(:,3),tval + dt);
        temp_p = reshape(temp_p,bezuit_bound,3);
        % WANT: x(dt) ~ x(0) + dx/dt (0) * dt
        %       dx/dt = - (dH/dx)^(-1) * dH/dt
        % -->  x_pred = dx/dt (0) * dt = - (dH/dx)^(-1) * dH/dt * dt
        xvaltemp = xval - dt*temp_p ; % Euler Predicted Value
        
        %fprintf('Before Corrector %d %d %d\n',xval(1,1),xval(1,2),xval(1,3));
        
        corrector_counter = 0;
        correctorSuccess = false;
        while corrector_counter < 3  && ~correctorSuccess
            % do corrector
            temp = corrector(xvaltemp(:,1),xvaltemp(:,2),xvaltemp(:,3),tval + dt);
            temp = reshape(temp,bezuit_bound,3);
            xvaltemp = xvaltemp - temp;
            corrector_counter = corrector_counter + 1;
            % below: switch to rel_error(Newton)
            fprintf('Residual xvaltemp: %f\n',norm(h(xvaltemp(:,1),xvaltemp(:,2),xvaltemp(:,3),tval + dt)));
            residual_big = (norm(h(xvaltemp(:,1),xvaltemp(:,2),xvaltemp(:,3),tval + dt)) >= correction_criterion);
            if residual_big
                corrector_counter = corrector_counter + 1;
            else 
                xval = xvaltemp;
                correctorSuccess = true;
            end
        end
        
        if correctorSuccess
            fprintf('success!\n');
            tval = tval + dt;
            successes = successes + 1;
            if successes > 2 
                dt = 2 * dt;
                successes = 0;
            end
        else
            fprintf('not success...\n');
            dt = dt/2;
            successes = 0;
        end

         tracking_values(:,:,counter) = xval;
         time(counter) = tval;

    %    Adaptive step size
        if mod(counter,1) == 0
            %condition_number =  cond(jac_eval(xval(1,1),xval(1,2),xval(1,3),tval));
            %fprintf('Condition Number  %d\n', condition_number);
            %fprintf('After Corrector %d %d %d\n',xval(1,1),xval(1,2),xval(1,3));
            fprintf('Iteration %d\n',counter)
            %dt = adaptive_stepsize_check(j,p,xval,tval,dt,h_t,s_t);
            fprintf('Residual xval: %f\n',norm(h(xval(:,1),xval(:,2),xval(:,3),tval)));
            dt
        end
        counter = counter + 1;
        
    end

end