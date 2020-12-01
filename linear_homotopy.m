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
    endgame_tolerance = 0.000001;
    num_solutions = bezuit_bound;
    convergence_range = 1e-10;
    checked = 1;
    
    while (1-tval)> convergence_range
        
        if (tval > 1 - endgame_tolerance) && (checked)
            fprintf('Endgame reached at time: %f\n',tval);
            [xval,num_solutions,tracking_values] = endgame(jac_eval,  xval,  tval + dt,tracking_values);
            checked = 0;
        end
               
        %fprintf('at time: %f\n',tval);
        
        %xval = rk45(p,xval,tval,dt,bezuit_bound);
        
        temp_p = predictor(xval(:,1),xval(:,2),xval(:,3),tval + dt);
        temp_p = reshape(temp_p,num_solutions,3);
        xvaltemp = xval - dt*temp_p ; % Euler Predicted Value
        
        % WANT: x(dt) ~ x(0) + dx/dt (0) * dt
        %       dx/dt = - (dH/dx)^(-1) * dH/dt
        % -->  x_pred = dx/dt (0) * dt = - (dH/dx)^(-1) * dH/dt * dt
     
        corrector_counter = 0;
        correctorSuccess = false;
        while corrector_counter < 3  && ~correctorSuccess
            %TODO: corrector switch to rel_error(Newton)
            temp = corrector(xvaltemp(:,1),xvaltemp(:,2),xvaltemp(:,3),tval + dt);
            temp = reshape(temp,num_solutions,3);
            %xvaltemp_nc = xvaltemp;
            xvaltemp = xvaltemp - temp;
            corrector_counter = corrector_counter + 1;

%            post_correct = norm(h(xvaltemp(:,1),xvaltemp(:,2),xvaltemp(:,3),tval + dt));
%            pre_correct = norm(h(xvaltemp_nc(:,1),xvaltemp_nc(:,2),xvaltemp_nc(:,3),tval + dt));            
%            rel_error = norm((xvaltemp(:,1) - xvaltemp_nc(:,1))/(xvaltemp(:,1)))            
%            residual_big = rel_error >= correction_criterion;
%            fprintf('Residual xvaltemp: %f\n',norm(h(xvaltemp(:,1),xvaltemp(:,2),xvaltemp(:,3),tval)));
            residual_big = norm(h(xvaltemp(:,1),xvaltemp(:,2),xvaltemp(:,3),tval + dt)) >= correction_criterion;
%           fprintf('Residual xvaltemp: %f\n',norm(h(xvaltemp(:,1),xvaltemp(:,2),xvaltemp(:,3),tval + dt)));
            
            if residual_big
                corrector_counter = corrector_counter + 1;
            else 
                xval = xvaltemp;
                correctorSuccess = true;
            end
        end
        
        if correctorSuccess
            tval = tval + dt;
            successes = successes + 1;
            if successes > 2 
                dt = 2 * dt;
                successes = 0;
            end
        else
            dt = dt/2;
            successes = 0;
        end

         tracking_values(:,:,counter) = xval;
         time(counter) = tval;
         
        if mod(counter,1) == 0
            %condition_number =  cond(jac_eval(xval(1,1),xval(1,2),xval(1,3),tval));
            %fprintf('Condition Number  %d\n', condition_number);
            %fprintf('After Corrector %d %d %d\n',xval(1,1),xval(1,2),xval(1,3));
            fprintf('Iteration %d\n',counter)
            %dt = adaptive_stepsize_check(j,p,xval,tval,dt,h_t,s_t);
            fprintf('Residual xval: %f\n',norm(h(xval(:,1),xval(:,2),xval(:,3),tval)));
            tval
        end
        counter = counter + 1;
        
    end
  
end