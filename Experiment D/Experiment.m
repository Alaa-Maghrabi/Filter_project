clear all
figure

class=2;

%Read the data from the files
fileID2= fopen('alssio1.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:6:length(B));
ycor=B(2:6:length(B));
ori = B(3:6:length(B));

if(class==2)
    xcor1=B(4:6:length(B));
    ycor1=B(5:6:length(B));
end

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

%initialize the measurement values
for i=1:1:floor(length(xcor))
    m(1,i)=xcor(i);
    n(1,i)=ycor(i);
    if(class==2)
        m(2,i)=xcor1(i);
        n(2,i)=ycor1(i);
    end
    z(1,i)=m(1,i);
    z(2,i)=n(1,i);
    z1(1,i)=m(2,i);
    z1(2,i)=n(2,i);
end

%define the particles
x_P= mvnrnd(0.5*(m(1,1)+m(2,1)),0.5,Nb_particles);
y_P= mvnrnd(0.5*(n(1,1)+n(2,1)),0.5,Nb_particles);
x_fin(1)=mean(x_P);
y_fin(1)=mean(y_P);

for i=1:1:length(X_real)-1 
    
    %real points
    X_=[X_real(i+1);Y_real(i+1)];
    
    %Compare the measurements to check the NLOS
    dist_meas=my_distance(z(:,i+1),z1(:,i+1),'L2');
    if(dist_meas>0.75)
        dist_meas = my_distance(z(:,i+1),z(:,i),'L2');
        dist_meas1= my_distance(z1(:,i+1),z(:,i),'L2');
        if(dist_meas>dist_meas1)
            m(1,i+1)=m(2,i+1);
            n(1,i+1)=n(2,i+1);
        else
            m(2,i+1)=m(1,i+1);
            n(2,i+1)=n(1,i+1);            
        end
    end
    
    %Particle filter
    [x_fin(i+1),y_fin(i+1),x_P,y_P,p_est,speed] = particle_fil_d(class,p_est...
        ,Nb_particles,m(:,i+1),n(:,i+1),x_P,y_P,speed);
         
    
    X_F=[x_fin(i+1);y_fin( i+1)];
    err0(i+1)=my_distance(X_,X_F,'L2');
    
    %scatter(X_real(i),Y_real(i),'g','filled', 'handlevisibility', 'on')
    scatter(m(1,i),n(1,i),'r','filled', 'handlevisibility', 'on')
    scatter(m(2,i),n(2,i),'r','filled', 'handlevisibility', 'on')
    quiver(x_fin(i),y_fin(i),cos(ori(i)-pi/2),sin(ori(i)-pi/2),'b','filled','handlevisibility','on')
    scatter(x_fin(i),y_fin(i),'b','filled', 'handlevisibility', 'off')
    plot([x_fin(i) X_real(i)],[y_fin(i) Y_real(i)],'b','LineWidth',1,'handlevisibility','off')
    %plot([m(i) X_real(i)],[n(i) Y_real(i)],'r','LineWidth',1,'handlevisibility','on')
    drawnow
    pause(0.25)
    cla
        
end