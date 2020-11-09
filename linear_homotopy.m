function [xval, tracking_values, time] = linear_homotopy(predictor, corrector, starting_solutions, n, bezuit_bound, correction_criterion, h_t, jac_eval,h)

    tracking_values = zeros(length(starting_solutions),length(starting_solutions(1,:)),n); % For plotting
    %(solutions, variables, h)
    dt = 1/n;
    time = zeros(1,n); % For plotting
    counter = 1;
    tval = 0;
    xval = double(starting_solutions);
    while tval < 1
        
        %xval = rk45(p,xval,tval,dt,bezuit_bound);
        temp_p = predictor(xval(:,1),xval(:,2),xval(:,3),tval + dt);
        temp_p = reshape(temp_p,bezuit_bound,3);
        xval = xval + temp_p ; % Euler Predicted Value
        
        %fprintf('Before Corrector %d %d %d\n',xval(1,1),xval(1,2),xval(1,3));
        
        temp = corrector(xval(:,1),xval(:,2),xval(:,3),tval+dt);
        temp = reshape(temp,bezuit_bound,3);
        xval_c = xval - temp ;% Corrected Value

        corrector_counter = 1;
        
        while corrector_counter < 3
            check = (norm(h(xval_c(1),xval_c(2),xval_c(3),tval + dt)) >= correction_criterion);
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

    %    Adaptive step size
        if mod(counter,1) == 0
            %condition_number =  cond(jac_eval(xval(1,1),xval(1,2),xval(1,3),tval));
            %fprintf('Condition Number  %d\n', condition_number);
            %fprintf('After Corrector %d %d %d\n',xval(1,1),xval(1,2),xval(1,3));
            fprintf('Iteration %d\n',counter)
            %dt = adaptive_stepsize_check(j,p,xval,tval,dt,h_t,s_t);
            fprintf('Residual: %f\n',norm(h(xval(1),xval(2),xval(3),tval + dt)));
            dt
        end
        counter = counter + 1;
        
    end

end