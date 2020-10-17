function [dt] = adaptive_stepsize_check(corrector,predictor,xval,tval,dt,h_t,s_t)

    xval_doublecheck = xval;
    tval_small = tval;
    dt_small = dt/2;
    xval_smallercheck = xval;
    dt_big = dt*2;
    tval_big = tval;

    for counter = 1:2
        %Two Iterations of half step
        xval_smallercheck = xval_smallercheck + dt_small*predictor(xval_smallercheck(1),xval_smallercheck(2),xval_smallercheck(3),tval+dt_small);% predicted Value
       % fprintf('correcter term%f\n'. corrector(xval_smallercheck(1),xval_smallercheck(2),xval_smallercheck(3),tval));
        xval_smallercheck = double(xval_smallercheck - corrector(xval_smallercheck(1),xval_smallercheck(2),xval_smallercheck(3),tval+dt_small)); % corrected Value
        tval_small = tval_small+ dt_small;   
        newtonStep = - corrector(xval_smallercheck(1),xval_smallercheck(2),xval_smallercheck(3),tval);
        rel_error_small = norm(newtonStep) / norm(xval_smallercheck);

        %Two Iterations of normal step
        xval = double(xval + dt*predictor(xval(1),xval(2),xval(3),tval+dt)); % Predicted Value
        xval = xval - corrector(xval(1),xval(2),xval(3),tval+dt);% Corrected Value
        tval = tval+ dt;   
        
        if counter == 1
            xval_reference = xval;
            tval_reference = tval;
        end
    end

    %One Iteration of double step
    % note: tval has been incremented by 2dt in previous loop
    fprintf('double pre-predictor%f\n',xval_reference); 
    xval_doublecheck = xval_doublecheck + dt_big*predictor(xval_doublecheck(1),xval_doublecheck(2),xval_doublecheck(3),tval);% predicted Value
           fprintf('double post-predictorf%f\n',xval_doublecheck);
    newtonStep = - corrector(xval_doublecheck(1),xval_doublecheck(2),xval_doublecheck(3),tval);
    rel_error_big = norm(newtonStep) / norm(xval_doublecheck);
    xval_doublecheck = double(xval_doublecheck + newtonStep); % corrected Value
    fprintf('double post-correctorf%f\n',xval_doublecheck); 
    
    fprintf('rel error big%f\n', rel_error_big);
    if rel_error_small > h_t %Comparison
        fprintf('xval ref%f\n',xval_reference); 
        fprintf('xval smaller check%f\n',xval_smallercheck);

        fprintf('Timestep halved to %f\n',dt);
        dt = dt_small;
    elseif rel_error_big < s_t
        fprintf('Timestep doubled to %f\n',dt);
        dt = dt_big;
    else
        fprintf('Timestep unchanged at %f\n',dt);
    end  
end