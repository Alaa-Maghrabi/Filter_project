function [x_fin,p_fin] =  kalm(m, n,x_est,p_est)
%KALM Kalman filter called in the Kalman-like Particle filter
%
%   input -----------------------------------------------------------------
%   
%       o x_est   : (1 x 1),  Previous estimation
%       o p_est   : (N x 1),  N-dimensional error covariance
%       o m   : (1 x 1),  X coordinate of the measurement
%       o n   : (1 x 1),  Y coordinate of the measurement
%
%   output ----------------------------------------------------------------
%
%       o x_fin   : (1 x 1),  coordinates of the estimation
%       o p_fin   : (1 x 1),  Y coordinates of the estimation
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%Define variables.
dt=0.2; 
HexAccel_noise_mag= 0.05;
u=0.00; %acceleration magnitude

% Initialize state transition matrix
A=[ 1 0 1 0;...     % [x]
    0 1 0 1;...     % [y]
    0 0 1 0;...     % [Vx]
    0 0 0 1 ];     % [Vy]
B=[(dt^2/2); (dt^2/2); 1; 1];
C=[1 0 0 0; 0 1 0 0];
Ez=[10 0;0 0.05];
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