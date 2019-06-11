function [X_real, Y_real] = draw_map()

hold on
X_real=[];
Y_real=[];


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

%plot the real trajectory
xx=[3 3];
yy=[4.5 1];
X_real= [X_real linspace(2.75,2.75,20)];
Y_real= [Y_real linspace(4.5,1,20)];
plot(xx,yy,'g','LineWidth', 2.5, 'handlevisibility', 'off');

xx1=[3 1];
yy1=[1 1];
X_real= [X_real linspace(2.75,0.75,20)];
Y_real= [Y_real linspace(1,1,20)];
plot(xx1,yy1,'g','LineWidth', 2.5, 'handlevisibility', 'off');

xx2=[1 1];
yy2=[1 4.5];
X_real= [X_real linspace(0.75,0.75,20)];
Y_real= [Y_real linspace(1,4.5,20)];
plot(xx2,yy2,'g','LineWidth', 2.5, 'handlevisibility', 'off');

xx3=[1 2];
yy3=[4.5 4.5];
X_real= [X_real linspace(0.75,1.75,10)];
Y_real= [Y_real linspace(4.5,4.5,10)];
plot(xx3,yy3,'g','LineWidth', 2.5, 'handlevisibility', 'off');

xx4=[2 2];
yy4=[4.5 2.75];
X_real= [X_real linspace(1.75,1.75,10)];
Y_real= [Y_real linspace(4.5,2.75,10)];
plot(xx4,yy4,'g','LineWidth', 2.5, 'handlevisibility', 'off');


end

% 

% 
% xx1=[2.75 0.75];
% yy1=[1 1];
% X_real= [X_real linspace(2.75,0.75,20)];
% Y_real= [Y_real linspace(1,1,20)];
% plot(xx1,yy1,'g','LineWidth', 2.5, 'handlevisibility', 'off');
% 
% xx2=[0.75 0.75];
% yy2=[1 4.5];
% X_real= [X_real linspace(0.75,0.75,20)];
% Y_real= [Y_real linspace(1,4.5,20)];
% plot(xx2,yy2,'g','LineWidth', 2.5, 'handlevisibility', 'off');
% 
% xx3=[0.75 1.75];
% yy3=[4.5 4.5];
% X_real= [X_real linspace(0.75,1.75,10)];
% Y_real= [Y_real linspace(4.5,4.5,10)];
% plot(xx3,yy3,'g','LineWidth', 2.5, 'handlevisibility', 'off');
% 
% xx4=[1.75 1.75];
% yy4=[4.5 2.75];
% X_real= [X_real linspace(1.75,1.75,10)];
% Y_real= [Y_real linspace(4.5,2.75,10)];
% plot(xx4,yy4,'g','LineWidth', 2.5, 'handlevisibility', 'off');