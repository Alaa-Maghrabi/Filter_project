function [x_fin,p_fin] =  kalm(m, n,x_est,p_est)
%kalman filter 

%Define variables.
dt=0.2; 
HexAccel_noise_mag= 0.092;
u=0.00; %acceleration magnitude

% Initialize state transition matrix
A=[ 1 0 1 0;...     % [x]
    0 1 0 1;...     % [y]
    0 0 1 0;...     % [Vx]
    0 0 0 1 ];     % [Vy]
B=[(dt^2/2); (dt^2/2); 1; 1];
C=[1 0 0 0; 0 1 0 0];
Ez=[0.3 0.1;0.1 0.3];
% Ex=[dt^4/4 0 dt^3/2 0;...
%     0 dt^4/4 0 dt^3/2; ...
%     dt^3/2 0 dt^2 0; ...
%     0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2;0.6
Ex=[dt^4/4 0 dt^2/2 0;...
    0 dt^4/4 0 dt^2/2; ...
    dt^2/2 0 1 0; ...
    0 dt^2/2 0 1].*HexAccel_noise_mag^2;

    z(1,1)=m;
    z(2,1)=n;
    
    
    % Predicted state and covariance
    x_prd(:,1) = A * x_est(:,1) + B*u;
    p_prd = A * p_est * A' + Ex;
    % Kalman gain
    K = p_prd*C'*inv(C*p_prd*C'+Ez);

    % Estimated state and covariance
    x_fin(:,1) = x_prd(:,1) + K * (z(:,1) - C * x_prd(:,1));
    p_fin = p_prd - K * C * p_prd;
    % Compute the estimated measurements
    y(:,1) = C * x_fin(:,1);

end