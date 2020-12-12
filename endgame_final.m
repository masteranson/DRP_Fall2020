function [xval,num_solutions, tracking_values] = endgame_final(h,  xval,  tval, tracking_values)
% possible trunctation criteria: does residual change if we replace t -<
% t+dt?

num_solutions = length(xval);
xval_new = zeros(num_solutions,1);
comparison_value = 1;
condition_threshold = 1e-1;
xval_size = size(xval);

for counter = 1:length(xval) %Truncating tracking solutions based on predictor condition number
    
    tval_residual = norm(h(xval(counter,1),xval(counter,2),xval(counter,3),tval));
    reference_residual = norm(h(xval(counter,1),xval(counter,2),xval(counter,3),comparison_value));
    
    rel_error = abs(tval_residual - reference_residual);
    
    if rel_error > condition_threshold
        xval_new(counter,1) = 1;
        num_solutions = num_solutions - 1;
        fprintf('solution staged for removal: ');
        for counterr = 1:length(xval_size(2))
            fprintf('%f%+fi ', [real(xval(counter,counterr)), imag(xval(counter,counterr))]);
        end
        fprintf('\nrelative error: %f\n', rel_error);
    end
end

vals = logical(~xval_new);
xval = xval(vals,:);
tracking_values = tracking_values(vals,:,:);
end