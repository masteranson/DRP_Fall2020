function [predictor, corrector] = parameter_homotopy(f_t, parameters, start_values, target_values, variables, patch)

starting_system = subs(f_t, parameters, start_values);
target_system = subs(f_t, parameters, target_values);
[predictor, corrector] = homotopy_generator([starting_system; patch], target_system, patch, variables);

end

