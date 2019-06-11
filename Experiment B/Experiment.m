%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Particle filter & Kalman filter

clear all

%Read the data from the files
fileID2= fopen('Alessio7.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:3:length(B));
ycor=B(2:3:length(B));

%Define variables.
dt=0.2;
Nb_particles=1500;
HexAccel_noise_mag=0.5;
Ex=[dt^4/4 0 dt^2/2 0;...
    0 dt^4/4 0 dt^2/2; ...
    dt^2/2 0 1 0; ...
    0 dt^2/2 0 1].*HexAccel_noise_mag^2;
P=Ex;

x_est = zeros(4, 1);
p_est = P;
err0(1)=0;

%Draw the map
[X_real,Y_real]=draw_map(length(xcor));

%Initialize the measurement values
for i=1:1:length(xcor)
    m(i)=xcor(i);
    n(i)=ycor(i);
    z(1,i)=m(i);
    z(2,i)=n(i);
    x_est(1,1)=m(1);
    x_est(2,1)=n(1);
    x_est(3,1)=0.0;
    x_est(4,1)=0.0;
    x_fin(1)=m(1);
    y_fin(1)=n(1);
end

%Define the particles
x_P= mvnrnd(x_fin,1,Nb_particles);
y_P= mvnrnd(y_fin,1,Nb_particles);


for i=1:1:length(X_real)-1
    %kalman filter
    [x_est(:,i+1),p_est] = kalman(m(i+1),n(i+1),x_est(:,i),p_est);
    
    %real points
    X_=[X_real(i+1);Y_real(i+1)];
    
    %Particle filter
    [x_fin(i+1),y_fin(i+1),x_P,y_P] = particle_filter(x_fin(i),y_fin(i),x_est(:,i+1)...
        ,Nb_particles,m(i+1),n(i+1),x_P,y_P);
    
    X_F=[x_fin(i+1);y_fin(i+1)];
    
    err0(i+1)=my_distance(X_,X_F,'L2');
    
    %scatter(X_real(i),Y_real(i),'g','filled', 'handlevisibility', 'on')
    %scatter(x_est(1,i),x_est(2,i),'r','filled', 'handlevisibility', 'off')
    %scatter(m(i),n(i),'y','filled', 'handlevisibility', 'off')
    scatter(x_fin(i),y_fin(i),'b','filled', 'handlevisibility', 'off')
    plot([x_fin(i) X_real(i)],[y_fin(i) Y_real(i)],'b','LineWidth',1,'handlevisibility','off')

    drawnow
    pause(0.25)
    cla

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kalman-Like Particle filter

clear all
figure

%Read data from the files
fileID2= fopen('onetag4.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:3:length(B));
ycor=B(2:3:length(B));

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

for j=1:Nb_particles
    p_est(:,:,j) = P;
end

speed = zeros(2,Nb_particles);

err0(1)=0;

%draw the map
[X_real,Y_real]=draw_map(length(xcor));

%Initialize the measurement values
for i=1:1:length(xcor)
    m(i)=xcor(i);
    n(i)=ycor(i);
    z(1,i)=m(i);
    z(2,i)=n(i);
end

%Define Paticles
x_P= mvnrnd(m(1),0.5,Nb_particles);
y_P= mvnrnd(n(1),0.5,Nb_particles);
x_fin(1)=mean(x_P);
y_fin(1)=mean(y_P);

for i=1:1:length(X_real)-1
    
    %real points
    X_=[X_real(i+1);Y_real(i+1)];
    
    %Particle filter
    [x_fin(i+1),y_fin(i+1),x_P,y_P,p_est,speed] = particle_fil(p_est...
        ,Nb_particles,m(i+1),n(i+1),x_P,y_P,speed);
    
    
    X_F=[x_fin(i+1);y_fin( i+1)];
    err0(i+1)=my_distance(X_,X_F,'L2');
    
    
    %scatter(X_real(i),Y_real(i),'g','filled', 'handlevisibility', 'on')
    %scatter(m(i),n(i),'r','filled', 'handlevisibility', 'off')
    scatter(x_fin(i),y_fin(i),'b','filled', 'handlevisibility', 'off')
    plot([x_fin(i) X_real(i)],[y_fin(i) Y_real(i)],'b','LineWidth',1,'handlevisibility','off')
    %plot([m(i) X_real(i)],[n(i) Y_real(i)],'r','LineWidth',1,'handlevisibility','on')
    drawnow
    pause(0.1)
    cla
    
    %legend('True trajectory','Measured data','particle filter 1','particle filter 2');
    
end


