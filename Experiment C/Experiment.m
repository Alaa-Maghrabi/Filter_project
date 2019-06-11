clear all
figure

%Read the data from the files
fileID2= fopen('tagontag3.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:6:length(B));
ycor=B(2:6:length(B));
xcor1=B(4:6:length(B));
ycor1=B(5:6:length(B));

%Define variables.
dt=0.2;
Nb_particles=150;
HexAccel_noise_mag=0.2;
Ex=[dt^4/4 0 dt^2/2 0;...
    0 dt^4/4 0 dt^2/2; ...
    dt^2/2 0 1 0; ...
    0 dt^2/2 0 1].*HexAccel_noise_mag^2;
P=Ex;

x_est = zeros(4, 1);
x_est1 = zeros(4, 1);

for j=1:Nb_particles
    p_est(:,:,j) = P;
    p_est1(:,:,j) = P;
end

speed = zeros(2,Nb_particles);
speed1= zeros(2,Nb_particles);
err0(1)=0;

%draw the map
[X_real,Y_real]=draw_map(length(xcor)-1);

%Initialize the measurement values
for i=1:1:floor(length(xcor))
    m(i)=xcor(i);
    n(i)=ycor(i);
    m1(i)=xcor1(i);
    n1(i)=ycor1(i);
    z(1,i)=m(i);
    z(2,i)=n(i);
end

%Define the particles
x_P= mvnrnd(m(1),0.5,Nb_particles);
y_P= mvnrnd(n(1),0.5,Nb_particles);
x_fin(1)=mean(x_P);
y_fin(1)=mean(y_P);

x_P1= mvnrnd(m1(1),0.5,Nb_particles);
y_P1= mvnrnd(n1(1),0.5,Nb_particles);
x_fin1(1)=mean(x_P1);
y_fin1(1)=mean(y_P1);

for i=1:1:length(xcor)-1
    
    %Particle filter for the moving tag
    [x_fin(i+1),y_fin(i+1),x_P,y_P,p_est,speed] = particle_fil(p_est...
        ,Nb_particles,m(i+1),n(i+1),x_P,y_P,speed);
    
    %Particle filter for the fixed tag
    [x_fin1(i+1),y_fin1(i+1),x_P1,y_P1,p_est1,speed1] = particle_fil(p_est1...
        ,Nb_particles,m1(i+1),n1(i+1),x_P1,y_P1,speed1);
        
    scatter(x_fin(i),y_fin(i),'b','filled', 'handlevisibility', 'off')
    scatter(x_fin1(i),y_fin1(i),'m','filled', 'handlevisibility', 'off')
    plot([x_fin(i) X_real(i)],[y_fin(i) Y_real(i)],'b','LineWidth',0.5,'handlevisibility','off')
    plot([x_fin1(i) 2.3],[y_fin1(i) 2.65],'m','LineWidth',0.5,'handlevisibility','off')
    drawnow
    pause(0.25)
    cla
        
end
%The fixed tag real position
scatter(2.3,2.65,'d','k','filled', 'handlevisibility', 'off')




