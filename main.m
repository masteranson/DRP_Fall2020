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

%% Parameter Homotopy
[predictor,corrector,jac,h] = parameter_homotopy(f_t,  coefficients, random_coefficients, target_coefficients,[x1,x3,z],r); 
[p_xval,p_tracking_values, p_time] = linear_homotopy(predictor, corrector, xval, n, length(xval), correction_criterion, h_t, jac, h);

%% Result Visualization

%Dehomogenize the variables
%xval =[xval(:,1) xval(:,2)]./(xval(:,3));
xval_size = size(xval);

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
    
     fprintf('Solution %d: ', counter);
     for counterr = 1:xval_size(2)
        fprintf('%f%+fi ', [real(xval(counter,counterr)), imag(xval(counter,counterr))]);
     end
     fprintf('\n');
     hold off;
 end