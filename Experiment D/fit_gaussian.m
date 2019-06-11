function []=fit_gaussian(x_P,y_P)
%FIT_GAUSSIAN Estimate the position using a gaussian distribution 
%
%   input -----------------------------------------------------------------
%   
%       o x_P   : (1 x Nb_particles),  X coordinate of the particles
%       o y_P   : (1 x Nb_particles),  Y coordinate of the particles
%
%   output ----------------------------------------------------------------
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%fit normal distribution to data

pd = fitdist(x_P,'Normal');

x_values = 0:0.01:4;
f = pdf(pd,x_values);
plot(x_values,f,'LineWidth',2)

pd = fitdist(y_P,'Normal');

x_values = 0:0.01:6;
f = pdf(pd,x_values);
plot(f,x_values,'LineWidth',2)