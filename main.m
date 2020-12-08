%% Initial Variables
n = 5; %steps
h_t = 0.001; %timestep half criterion upper bound
s_t = 1e-7;
max_correction = 3;
correction_criterion = 1e-10;
chemsys %Load config file

%% Linear Homotopy
[predictor, corrector,jac,h] = homotopy_generator(fs, f, r, [x1,x3,z], t);
hval=h(starting_solutions(:,1),starting_solutions(:,2),starting_solutions(:,3),0);
fprintf('Residual: %f\n',norm(hval));
[xval,tracking_values,time,tval] = linear_homotopy(predictor, corrector, starting_solutions, n, bezuit_bound, correction_criterion, h_t,jac,h);


% possible trunctation criteria: does residual change if we replace t -<
% t+dt?
sizexval = size(xval)
k=sizexval(1)
for i=1:k
  [norm(h(xval(i,1),xval(i,2),xval(i,3),tval)), norm(h(xval(i,1),xval(i,2),xval(i,3),1.0))]
end





%% Parameter Homotopy
[predictor,corrector,jac,h] = parameter_homotopy(f_t,  coefficients, random_coefficients, target_coefficients,[x1,x3,z],r); 
[p_xval,p_tracking_values, p_time] = linear_homotopy(predictor, corrector, xval, n, length(xval), correction_criterion, h_t, jac, h);

%% Result Visualization

tracked_solutions = length(xval);
fprintf('Tracked %d solutions\n',tracked_solutions);
figure
 for counter = 1:tracked_solutions
     
    subplot(3,3,counter)
    title(sprintf('Homotopy Tracking Path Number %d',counter));
    xlabel('\Re');
    ylabel('\Im');
    zlabel('T');
    hold on   
    
     for solutions = 1:3
        real_part = squeeze(real(tracking_values(counter,solutions,:)));
        imag_part = squeeze(imag(tracking_values(counter,solutions,:)));
        plot3(real_part,imag_part,time,'LineWidth',2);
     end
     
    hold off
     
    view(3)
    
%     hold on
%     max_val = max(norm(tracking_values(:)));
%     [x, y] = meshgrid(-max_val*2:1:max_val*2);
%     z = zeros(length(x),length(y));
%     surf(x,y,z);
%     z = ones(length(x),length(y));
%     surf(x,y,z);

%      fprintf('Solution %d: ', counter);
%      for counterr = 1: length(xval)
%         fprintf('%f%+fi ', [real(xval(counterr,counter)), imag(xval(counterr,counter))]);
%      end
%      fprintf('\n');
%      hold off;
 end