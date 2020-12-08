function [xval,num_solutions, tracking_values] = endgame(jac_eval,  xval,  tval, tracking_values)

condition_threshold = 10^6;
num_solutions = length(xval);
xval_new = xval;
tracking_values_new = tracking_values;

for counter = 1:length(xval) %Truncating tracking solutions based on predictor condition number
    condition_number = cond(jac_eval(xval(counter,1),xval(counter,2),xval(counter,3),tval));
    if condition_number > condition_threshold
        xval_new(counter,:) = [];
        num_solutions = num_solutions - 1;
        tracking_values_new(counter,:,:) = [];
        fprintf('solution staged for removal: %f\n', xval(counter,1:3));
        fprintf('log condition number: %f\n', log10(condition_number));
    end
end

xval = xval_new;
tracking_values = tracking_values_new;
end