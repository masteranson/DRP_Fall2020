function [dt] = adaptive_stepsize_check(corrector,predictor,xval,tval,dt,h_t,s_t)

    xval_doublecheck = xval;
    tval_small = tval;
    dt_small = dt/2;
    xval_smallercheck = xval;
    dt_big = dt*2;
    tval_big = tval;

    for counter = 1:2
        %Two Iterations of half step
        xval_smallercheck = double(xval_smallercheck + corrector(xval_smallercheck(1),xval_smallercheck(2),xval_smallercheck(3),tval)*dt_small); % Predicted Value
        xval_smallercheck = xval_smallercheck - predictor(xval_smallercheck(1),xval_smallercheck(2),xval_smallercheck(3),tval+dt_small);% Corrected Value
        tval_small = tval_small+ dt_small;   

        %Two Iterations of normal step
        xval = double(xval + predictor(xval(1),xval(2),xval(3),tval)*dt); % Predicted Value
        xval = xval - corrector(xval(1),xval(2),xval(3),tval+dt);% Corrected Value
        tval = tval+ dt;   
        
        if counter == 1
            xval_reference = xval;
        end
    end

    %One Iteration of double step
    xval_doublecheck = double(xval_doublecheck + corrector(xval_doublecheck(1),xval_doublecheck(2),xval_doublecheck(3),tval_big)*dt_big); % Predicted Value
    xval_doublecheck = xval_doublecheck - predictor(xval_doublecheck(1),xval_doublecheck(2),xval_doublecheck(3),tval_big+dt_big);% Corrected Value


    if norm(xval_reference-xval_smallercheck)/norm(xval_reference) > h_t %Comparison
        fprintf('Timestep halved to %f\n',dt);
        dt = dt_small;
    elseif norm(xval_reference-xval_doublecheck)/norm(xval_reference) < s_t
        fprintf('Timestep doubled to %f\n',dt);
        dt = dt_big;
    else
        fprintf('Timestep unchanged at %f\n',dt);
    end  
end