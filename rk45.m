%TODO: FIX!!!
function [xval] = rk45(fn,xval,t,dt,bezuit_bound)

k_1 = fn(xval(:,1),xval(:,2),xval(:,3),t);
k_1 = reshape(k_1,bezuit_bound,3);
k_2 = fn((xval(:,1) + k_1(:,1)*dt/2),(xval(:,2) + k_1(:,2)*dt/2),(xval(:,3) + k_1(:,3)*dt/2), t+ dt/2);
k_2 = reshape(k_2,bezuit_bound,3);
k_3 = fn((xval(:,1) + k_2(:,1)*dt/2),(xval(:,2) + k_2(:,2)*dt/2),(xval(:,3) + k_2(:,3)*dt/2), t+ dt/2);
k_3 = reshape(k_3,bezuit_bound,3);
k_4 = fn((xval(:,1) + k_3(:,1)*dt),(xval(:,2) + k_3(:,2)*dt),(xval(:,3) + k_3(:,3)*dt), t+ dt);
k_4 = reshape(k_4,bezuit_bound,3);
xval = xval + (dt/6)*(k_1 + 2*k_2 + 2*k_3 +k_4);
end