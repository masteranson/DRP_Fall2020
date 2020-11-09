%% Initial Variables
n = 5000; %steps
dt = 1/n; %Starting Timestep
h_t = 0.001; %timestep half criterion upper bound
s_t = 1e-7;
max_correction = 2;
correction_criterion = 0.0001;
chemsys %Load config file

%% Linear Homotopy
[predictor, corrector] = homotopy_generator(fs, f, r, [x1,x3,z]);
[xval,tracking_values,time] = linear_homotopy(predictor, corrector, double(starting_solutions), n, bezuit_bound, correction_criterion, h_t);

%% Parameter Homotopy
[predictor,corrector] = parameter_homotopy(f_t,  coefficients, random_coefficients, target_coefficients,[x1,x3,z],r); 
p_xval = linear_homotopy(predictor, corrector, xval, n, bezuit_bound, correction_criterion, h_t);

%% Result Visualization

% tracked_solutions = length(xval);
% fprintf('Tracked %d solutions\n',tracked_solutions);
% 
%  for counter = 1:tracked_solutions
%      for solutions = 1:3
%         plot3(real(tracking_values(counter,solutions,:)),imag(tracking_values(counter,solutions,:)),time,'LineWidth',2);
%      end
%      
%     view(3)
% 
%     figure(tracked_solutions)
%     title('Homotopy Tracking Path');
%     xlabel('\Re');
%     ylabel('\Im');
%     zlabel('T');
%     hold on
%     max_val = max(norm(tracking_values(:)));
%     [x, y] = meshgrid(-max_val*2:1:max_val*2);
%     z = zeros(length(x),length(y));
%     surf(x,y,z);
%     z = ones(length(x),length(y));
%     surf(x,y,z);
% 
%      fprintf('Solution %d: ', counter);
%      for counterr = 1: length(xval)
%         fprintf('%f%+fi ', [real(xval(counterr,counter)), imag(xval(counterr,counter))]);
%      end
%      fprintf('\n');
%      hold off;
%  end