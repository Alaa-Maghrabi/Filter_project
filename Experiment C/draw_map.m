function [X_real, Y_real] = draw_map(size)
%DRAW_MAP draw the map for the experiment, contains also the 4 persons
%location and the tag's location
%
%   input -----------------------------------------------------------------
%   
%       o size   : (1 x 1),  Size of the datapoint
%
%   output ----------------------------------------------------------------
% 
%       o X_real   : (size x 1),  real trajectory coordinate X
%       o Y_real   : (size x 1),  real trajectory coordinate Y
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold on
X_real=[];
Y_real=[];


%plot the four anchors
scatter(0,0,'s','filled', 'handlevisibility', 'off')
scatter(3.75,0.05,'s','filled', 'handlevisibility', 'off')
scatter(3.93,5.81,'s','filled', 'handlevisibility', 'off')
scatter(0.05,5.76,'s','filled', 'handlevisibility', 'off')

%plot second tag
scatter(2.3,2.65,'d','k','filled', 'handlevisibility', 'off')

%plot path
xx=[2 2];
yy=[4.5 1];
X_real= [X_real linspace(2,2,size)];
Y_real= [Y_real linspace(4.5,1,size)];
plot(xx,yy,'g','LineWidth', 2.5, 'handlevisibility', 'off');

end