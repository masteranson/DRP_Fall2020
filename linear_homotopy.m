function [xval, tracking_values, time] = linear_homotopy(predictor, corrector, starting_solutions, n, bezuit_bound, correction_criterion, h_t)

    tracking_values = zeros(length(starting_solutions),length(starting_solutions(1,:)),n); % For plotting
    %(solutions, variables, h)
    dt = 1/n;
    time = zeros(1,n); % For plotting
    counter = 1;
    tval = 0;
    xval = double(starting_solutions);
    while tval < 1
        
        %xval = rk45(p,xval,tval,dt,bezuit_bound);
        temp_p = predictor(xval(:,1),xval(:,2),xval(:,3),tval + dt*dt);
        temp_p = reshape(temp_p,bezuit_bound,3);
        xval = xval + temp_p ; % Euler Predicted Value

        temp = corrector(xval(:,1),xval(:,2),xval(:,3),tval+dt);
        temp = reshape(temp,bezuit_bound,3);
        xval_c = xval - temp ;% Corrected Value

        corrector_counter = 1;
        
        while corrector_counter <= 3
            check = (abs(norm(xval)-norm(xval_c)) >= correction_criterion);
            if check 
                xval_c = xval;
                temp = corrector(xval(:,1),xval(:,2),xval(:,3),tval+dt);
                temp = reshape(temp,bezuit_bound,3);
                xval = xval - temp;
                corrector_counter = corrector_counter + 1;
            elseif ~check && corrector_counter == 1
                xval = xval_c;
                corrector_counter = 69;
            end
        end 

    %     diagonal_vecs = diag((1./vecnorm(xval')));
    %     xval = diagonal_vecs*xval;

         tval = tval + dt;   
         tracking_values(:,:,counter) = xval;
         time(counter) = tval;

    %    Adaptive step size
        if mod(counter,20) == 0
            fprintf('Iteration %d\n',counter)
            %dt = adaptive_stepsize_check(j,p,xval,tval,dt,h_t,s_t);
            %fprintf('Residual: %f\n',norm(h(xval(1),xval(2),xval(3),tval)));
        end
        counter = counter + 1;
    end

end