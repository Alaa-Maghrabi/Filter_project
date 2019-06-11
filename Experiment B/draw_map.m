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

Section = size/8;

%plot the four anchors
scatter(0,0,'s','filled', 'handlevisibility', 'off')
scatter(3.75,0.05,'s','filled', 'handlevisibility', 'off')
scatter(3.93,5.81,'s','filled', 'handlevisibility', 'off')
scatter(0.05,5.76,'s','filled', 'handlevisibility', 'off')

%plot the four people
rectangle('Position',[1.15 3.7 0.3 0.25], 'handlevisibility', 'off')
rectangle('Position',[2.2 3.7 0.3 0.25], 'handlevisibility', 'off')
rectangle('Position',[1.15 1.3 0.3 0.25], 'handlevisibility', 'off')
rectangle('Position',[2.2 1.3 0.3 0.25], 'handlevisibility', 'off')

% %plot the real trajectory
xx=[3 3];
yy=[4.5 1];
X_real= [X_real linspace(3,3,floor(2*Section))];
Y_real= [Y_real linspace(4.5,1,floor(2*Section))];
plot(xx,yy,'g','LineWidth', 2.5, 'handlevisibility', 'off');

xx1=[3 1];
yy1=[1 1];
X_real= [X_real linspace(3,1,floor(2*Section))];
Y_real= [Y_real linspace(1,1,floor(2*Section))];
plot(xx1,yy1,'g','LineWidth', 2.5, 'handlevisibility', 'off');

xx2=[1 1];
yy2=[1 4.5];
X_real= [X_real linspace(1,1,floor(2*Section))];
Y_real= [Y_real linspace(1,4.5,floor(2*Section))];
plot(xx2,yy2,'g','LineWidth', 2.5, 'handlevisibility', 'off');

xx3=[1 2];
yy3=[4.5 4.5];
X_real= [X_real linspace(1,2,floor(1*Section))];
Y_real= [Y_real linspace(4.5,4.5,floor(1*Section))];
plot(xx3,yy3,'g','LineWidth', 2.5, 'handlevisibility', 'off');

xx4=[2 2];
yy4=[4.5 2.75];
X_real= [X_real linspace(2,2,floor(1*Section))];
Y_real= [Y_real linspace(4.5,2.75,floor(1*Section))];
plot(xx4,yy4,'g','LineWidth', 2.5, 'handlevisibility', 'off');

end