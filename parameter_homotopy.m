function [predictor, corrector, jac, h] = parameter_homotopy(f_t, parameters, start_values, target_values, variables, patch)

syms t

starting_system = subs(f_t, parameters, start_values);
target_system = subs(f_t, parameters, target_values);
[predictor, corrector, jac, h] = homotopy_generator([starting_system; patch], target_system, patch, variables, t);

end

