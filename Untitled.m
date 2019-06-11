clear all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Difference between LOS and NLOS
%test 4 LOS with some NLOS data
%test 6 NLOS data (3 people blocking three anchors)
%test 3 LOS data

%fileID= fopen('test5.txt','r');
%formatSpec = '%*s %f %*s %f %*s %f %*s';
%A=fscanf(fileID,formatSpec);
xx=[2.45,3.97];
yy=[0.02,0.02];
plot(xx,yy,'r')
%x=A(1:3:length(A));
%y=A(2:3:length(A));
%a=mean(x);
%b=mean(y);
%aLOS=2.2528;
%bLOS=3.0483;
xx(8)= 3.82
yy(8)= 0.43
xx(7)= 3.71
yy(7)= -0.44
xx(1)= 3.59
yy(1)= 0.54
xx(2)= 3.45
yy(2)= 0.69
xx(3)= 3.20
yy(3)= 0.48
xx(4)= 3.25
yy(4)= -0.62
x(1)=2.45
y(1)=-0.01
x(2)=2.59
y(2)=0.04
x(3)= 3.11
y(3)= 0.02
x(4)= 3.03
y(4)= 0.01
xx(5)= 3.3
yy(5)= 0.35
x(5)= 2.59
y(5)= -0.02
x(6)= 2.94
y(6)= 0.01
x(7)= 2.82
y(7)= 0.07
x(8)= 2.71
y(8)= -0.01
xx(6)= 3.93
yy(6)= -0.27
% x(10)= 4.07
% y(10)= -0.3
% x(11)= 4.19
% y(11)= -0.19
%sigmaa = 0.00055336;
%sigmab = 0.0002645;


hold on
scatter(x,y,'b')
scatter(xx,yy,'m')
legend('True trajectory','Measured data LOS','Measured data NLOS')
%scatter(mvnrnd(a,sigmaa,200),mvnrnd(b,sigmab,200),'g')
%scatter(aLOS,bLOS,'y','filled')
%scatter(a,b,'r','filled')

hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fir the data into a gaussian mixture

X = [x  y];
figure
scatter(X(:,1),X(:,2))
hold on
scatter(a,b,'r','filled')

options = statset('Display','final');
gm = fitgmdist(X,2,'Options',options);
gmPDF = @(x,y)pdf(gm,[x y]);

hold on
h = ezcontour(gmPDF,[2.1 2.4],[3 3.1]);
title('Scatter Plot and PDF Contour')
hold off



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %fit normal distribution to data

pd = fitdist(x,'Normal')

x_values = 0:0.01:5;
f = pdf(pd,x_values);
figure
plot(x_values,f,'LineWidth',2)

pd = fitdist(y,'Normal')

x_values = -1:0.001:1;
f = pdf(pd,x_values);
figure
plot(x_values,f,'LineWidth',2)

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for i=1:1:20%floor(length(xcor)/Nb)
% %     m(i) = mean(xcor((i-1)*Nb+1:i*Nb));
% %     n(i) = mean(ycor((i-1)*Nb+1:i*Nb));
%     m(i)=mvnrnd(X_real(i),0.08,1);
%     n(i)=mvnrnd(Y_real(i),0.08,1);
%     z(1,i)=m(i);
%     z(2,i)=n(i);
%     x_est(1,1)=m(1);
%     x_est(2,1)=n(1);
%     x_est(3,1)=0.0;
%     x_est(4,1)=0.0;
%     x_fin(1)=m(1);
%     y_fin(1)=n(1);
%     formatSpec = 'x: %1.11f \ny: %1.11f \nz: 0.20000000298 \n---\n';
%     fprintf(formatSpec,m(i),n(i))
% end-
% for i=20:1:70%floor(length(xcor)/Nb)
% %     m(i) = mean(xcor((i-1)*Nb+1:i*Nb));
% %     n(i) = mean(ycor((i-1)*Nb+1:i*Nb));
%     m(i)=mvnrnd(X_real(i),0.08,1);
%     n(i)=mvnrnd(Y_real(i),0.08,1);
%     z(1,i)=m(i);
%     z(2,i)=n(i);
%     x_est(1,1)=m(1);
%     x_est(2,1)=n(1);
%     x_est(3,1)=0.0;
%     x_est(4,1)=0.0;
%     x_fin(1)=m(1);
%     y_fin(1)=n(1);
%     formatSpec = 'x: %1.11f \ny: %1.11f \nz: 0.20000000298 \n---\n';
%     fprintf(formatSpec,m(i),n(i))
% end
% for i=70:1:80%floor(length(xcor)/Nb)
% %     m(i) = mean(xcor((i-1)*Nb+1:i*Nb));
% %     n(i) = mean(ycor((i-1)*Nb+1:i*Nb));
%     m(i)=mvnrnd(X_real(i),0.08,1);
%     n(i)=mvnrnd(Y_real(i),0.08,1);
%     z(1,i)=m(i);
%     z(2,i)=n(i);
%     x_est(1,1)=m(1);
%     x_est(2,1)=n(1);
%     x_est(3,1)=0.0;
%     x_est(4,1)=0.0;
%     x_fin(1)=m(1);
%     y_fin(1)=n(1);
%     formatSpec = 'x: %1.11f \ny: %1.11f \nz: 0.20000000298 \n---\n';
%     fprintf(formatSpec,m(i),n(i))
% end
%%
clear all

sigmaa = 0.00055336;
sigmab = 0.0002645;

fileID2= fopen('test5.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);

x=B(1:3:length(B));
y=B(2:3:length(B));
figure
hold on
x_N=0.00;
%x_P = ones(26,20);
%y_P = ones(26,20);

for i=1:1:26
    m(i)= mean(x((i-1)*12+1:i*12));
    n(i) = mean(y((i-1)*12+1:i*12));
    %x_P(i) = m(i) +  sqrt(x_N)*randn;
    %y_P(i) = n(i) + sqrt(x_N)*randn;
    
    
    scatter(mvnrnd(m(i),sigmaa,20),mvnrnd(n(i),sigmab,20),'g')
    scatter(m(i),n(i),'r','filled')
end

%figure
%scatter(x,y,'b')
%hold on
%scatter(m,n,'r','filled')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Particle Filter
clear all

sigmaa = 0.0015;
sigmab = 0.013;

fileID2= fopen('test5.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:3:length(B));
y=B(2:3:length(B));

Nb_particles=300;
x_N=0.0; %covariance matrix

%define the vector of particles
Weight=1/Nb_particles*ones(Nb_particles,1);
Final_Weight= 1/Nb_particles*ones(Nb_particles,1);
Nb=1;

for i=1:1:floor(length(xcor)/Nb)
    m(i) = xcor(i); %mean(xcor((i-1)*Nb+1:i*Nb));
    n(i) = y(i);%2*mean(y((i-1)*Nb+1:i*Nb));
end


%m(1)=3.78;
%n(1)=0;
x_est(1)=m(1);
y_est(1)=n(1);

%x_P= mvnrnd(m(1),sigmaa,Nb_particles);
%y_P= mvnrnd(n(1),sigmab,Nb_particles);

for i=1:1:floor(length(xcor)/Nb)-1
     
        X_1=[x_est(i) ; y_est(i)];
        X_2=[m(i) ; n(i)];
        
    if( my_distance(X_1, X_2, 'L2') < 0.1)
        est(i)=1;
    else
        est(i)=0;
    end
    
    x_P= mvnrnd(x_est(i),sigmaa,Nb_particles);
    y_P= mvnrnd(y_est(i),sigmab,Nb_particles);
    
%      figure
%      hold on
%      scatter(x_P,y_P,'b','filled')
%      scatter(x_est(i),y_est(i),'r','filled')
%      scatter(m,n,'y','filled')
%      hold off
    
    %We get the new estimation
    for j=1:1:Nb_particles
        %x_P(j) = x_P(j) + sqrt(x_N)*randn;
        %y_P(j) = y_P(j) + sqrt(x_N)*randn;
        %%Weight(j) = exp(-(x_P(j)-m(i+1))^2/(2*sigmaa) - (y_P(j)-n(i+1))^2/(2*sigmab));
        dist(j)=(x_P(j)-m(i+1))^2 + (y_P(j)-n(i+1))^2;
    end
    
    [~,index]=min(dist);
    
    %%if(sum(Weight)>0)
        %%Final_Weight=Weight/sum(Weight);
            %x_est(i+1)=x_P'*Final_Weight;
            %y_est(i+1)=y_P'*Final_Weight;
            %else
        
            %x_est(i+1)=2*x_est(i)-x_est(i-1);
            %y_est(i+1)=2*y_est(i)-y_est(i-1);
    %%end
    
    %%x_est(i+1)=x_P'*Final_Weight;
    %%y_est(i+1)=y_P'*Final_Weight;


                        x_est(i+1)= x_P(index);
                        y_est(i+1)= y_P(index);
%      figure
      hold on
      %scatter(x_P,y_P,'g')
      scatter(x_est(i+1),y_est(i+1),'r','filled')
      scatter(m(i+1),n(i+1),'y','filled')
      pause;
%      hold off
end
% figure
% hold on
% scatter(m,n,'b')
% scatter(x_est,y_est,'r','filled')
% %scatter(2,2)
% hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kalmann filter
clear all
close all

fileID2= fopen('humana2.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:3:length(B));
ycor=B(2:3:length(B));

%Define variables.
dt=0.2;
Nb_particles=400;
HexAccel_noise_mag=0.2;
% Ex=[dt^4/4 0 dt^3/2 0;...
%     0 dt^4/4 0 dt^3/2; ...
%     dt^3/2 0 dt^2 0; ...
%     0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2;
Ex=[dt^4/4 0 dt^2/2 0;...
    0 dt^4/4 0 dt^2/2; ...
    dt^2/2 0 1 0; ...
    0 dt^2/2 0 1].*HexAccel_noise_mag^2;
P=Ex;

% x_est=[x,y,Vx,Vy]'
x_est = zeros(4, 1);
p_est = P;

j=1;
k=1;
l=1;
ii=1;
Nb=1;
bool0=1;
bool=0;
bool1=0;
bool2=0;
bool3=0;
Weight_par0=0;
Weight_par=0;
Weight_par1=0;
Weight_par2=0;
Weight_par3=0;
err0(1)=0;
err1(1)=0;
err2(1)=0;
err3(1)=0;
err4(1)=0;
err_raw(1)=0;

for i=1:1:floor(length(xcor)/Nb)
    m(i) = mean(xcor((i-1)*Nb+1:i*Nb));
    n(i) = mean(ycor((i-1)*Nb+1:i*Nb));
    z(1,i)=m(i);
    z(2,i)=n(i);
    x_est(1,1)=m(1);
    x_est(2,1)=n(1);
    x_est(3,1)=0.0;
    x_est(4,1)=0.0;
    x_fin(1)=m(1);
    y_fin(1)=n(1);
end

x_P= mvnrnd(x_fin,1,Nb_particles);
y_P= mvnrnd(y_fin,1,Nb_particles);

x_P1= mvnrnd(x_fin,1,Nb_particles);
y_P1= mvnrnd(y_fin,1,Nb_particles);

x_P2= mvnrnd(x_fin,1,Nb_particles);
y_P2= mvnrnd(y_fin,1,Nb_particles);

x_P3= mvnrnd(x_fin,1,Nb_particles);
y_P3= mvnrnd(y_fin,1,Nb_particles);

x_P4= mvnrnd(x_fin,1,Nb_particles);
y_P4= mvnrnd(y_fin,1,Nb_particles);
%test 9
% xval = 2.15;
% ymin = 1;
% ymax = 5.7;

% test human
xval = 1.95;
ymin = 1.2;
ymax = 4.5;

%test 23
% xval = 2.2;
% ymin = 1.1;
% ymax = 5.5;

%test 30
% xval = 2.25;
% ymin = 1.2;
% ymax = 4.6;

%test 22
% xval = 2.15;
% ymin = 0.6;
% ymax = 5.3;


hold on

%plot the real trajectory
xx=[xval,xval];
yy=[ymin,ymax];
%plot(xx,yy)

y_real=linspace(ymin,ymax,floor(length(xcor)/Nb));

%plot the four anchors
scatter(0,0,'s','filled')
scatter(3.75,0.05,'s','filled')
scatter(3.93,5.81,'s','filled')
scatter(0.05,5.76,'s','filled')

%plot the four people
rectangle('Position',[1.25 4.5 0.4 0.25])
rectangle('Position',[2.6 4.5 0.4 0.25])
rectangle('Position',[1.25 1 0.4 0.25])
rectangle('Position',[2.6 1 0.4 0.25])

for i=1:1:floor(length(xcor)/Nb)-1
    %kalman filter
    [x_est(:,i+1),p_est] = kalm(m(i),n(i),x_est(:,i),p_est);
    
    %real points
    X_real=[xval;y_real(i+1)];
    err_raw(i+1)=my_distance(X_real,z(:,i+1),'L2');

    %Particle filter
    if(bool0==1)
        [x_fin(i+1),y_fin(i+1),x_P,y_P] = particle_fil(x_fin(i),y_fin(i),x_est(:,i+1)...
            ,Nb_particles,m(i+1),n(i+1),x_P,y_P);

        X_F=[x_fin(i+1);y_fin(i+1)];
        
        if(my_distance(X_F,z(:,i+1),'L2')>0.9)
            Weight_par0=Weight_par0+1;
        end
        if(Weight_par0>10)
            bool0=0;
        end
        
        err0(i+1)=my_distance(X_real,X_F,'L2');
    end
    
    %Creates several particle filter
    if( i>=2 && my_distance(z(:,i-1),z(:,i),'L2')>0.25  && i<6)
        if (i==2)
            bool=1;
            x_fin1(1)=m(i);
            y_fin1(1)=n(i);
        end
        if (i==3)
            bool1=1;
            for j=1:1:i-1
                if(my_distance(z(:,j),z(:,i),'L2')<0.25)
                    bool1=0;
                    break;
                end
            end
            if(bool1==1)
                x_fin2(1)=m(i);
                y_fin2(1)=n(i);
            end
        end
        if (i==4)
            bool2=1;
            for j=1:1:i-1
                if(my_distance(z(:,j),z(:,i),'L2')<0.25)
                    bool2=0;
                    break;
                end
            end
            if(bool2==1)
                x_fin3(1)=m(i);
                y_fin3(1)=n(i);
            end
        end
        if (i==5)
            bool3=1;
            for j=1:1:i-1
                if(my_distance(z(:,j),z(:,i),'L2')<0.25)
                    bool3=0;
                    break;
                end
            end
            if(bool3==1)
                x_fin4(1)=m(i);
                y_fin4(1)=n(i);
            end
            
        end
        
    end
    
    
    if(bool==1)
        [x_fin1(j+1),y_fin1(j+1)] = particle_fil(x_fin1(j),y_fin1(j)...
            ,x_est(:,i+1),Nb_particles,m(i+1),n(i+1),x_P1,y_P1);
        X_F1=[x_fin1(j+1) y_fin1(j+1)];
        if(my_distance(X_F1,z(:,i+1),'L2')>0.7)
            Weight_par=Weight_par+1;
        end
        if(Weight_par>10)
            bool=0;
        end
                
        err1(j+1)=  my_distance(X_real,X_F1,'L2');
        j=j+1;
    end
    
    if(bool1==1)
        [x_fin2(k+1),y_fin2(k+1)] = particle_fil(x_fin2(k),y_fin2(k)...
            ,x_est(:,i+1),Nb_particles,m(i+1),n(i+1),x_P2,y_P2);
        X_F2=[x_fin2(k);y_fin2(k)];
        
        if(my_distance(X_F2,z(:,i),'L2')>0.7)
            Weight_par1=Weight_par1+1;
        end
        
        if(Weight_par1>10)
            bool1=0;
        end
        
        err2(k+1)=  my_distance(X_real,X_F2,'L2');
        
        k=k+1;
    end
    
    if(bool2==1)
        [x_fin3(l+1),y_fin3(l+1)] = particle_fil(x_fin3(l),y_fin3(l)...
            ,x_est(:,i+1),Nb_particles,m(i+1),n(i+1),x_P3,y_P3);
        X_F3=[x_fin3(l+1);y_fin3(l+1)];
        
        if(my_distance(X_F3,z(:,i+1),'L2')>0.7)
            Weight_par2=Weight_par2+1;
        end
        
        if(Weight_par2>10)
            bool2=0;
        end
        
        err3(l+1)=  my_distance(X_real,X_F3,'L2');
        l=l+1;

    end
    
    if(bool3==1)
        [x_fin4(ii+1),y_fin4(ii+1)] = particle_fil(x_fin4(ii),y_fin4(ii)...
            ,x_est(:,i+1),Nb_particles,m(i+1),n(i+1),x_P4,y_P4);
        X_F4=[x_fin4(ii+1);y_fin4(ii+1)];
        
        if(my_distance(X_F4,z(:,i+1),'L2')>0.7)
            Weight_par3=Weight_par3+1;
        end
        
        if(Weight_par3>10)
            bool3=0;
        end
        
        
        err4(ii+1)=  my_distance(X_real,X_F4,'L2');
        ii=ii+1;
    end
    
    
    scatter(xval,y_real(i),'g','filled')
    %scatter(x_est(1,i),x_est(2,i),'r','filled')
    scatter(m(i),n(i),'y','filled')
    
    
    if(bool0==1)
        scatter(x_fin(i),y_fin(i),'b','filled')
    end
    if(bool==1)
        %scatter(x_fin1(j-1),y_fin1(j-1),'g','filled')
    end
    if(bool1==1)
        %scatter(x_fin2(k-1),y_fin2(k-1),'c','filled')
    end
    if(bool2==1)
        %scatter(x_fin3(l-1),y_fin3(l-1),'k','filled')
    end
    if(bool3==1)
        %scatter(x_fin4(ii-1),y_fin4(ii-1),'m','filled')
    end
    pause;
    %legend('True trajectory','Measured data','particle filter 1','particle filter 2');

end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all



fileID2= fopen('test10.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:3:length(B));
ycor=B(2:3:length(B));
xcor=mean(xcor);
ycor=mean(ycor);

fileID2= fopen('test11.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

fileID2= fopen('test12.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

fileID2= fopen('test13.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

fileID2= fopen('test14.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

fileID2= fopen('test15.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

fileID2= fopen('test16.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

fileID2= fopen('test17.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

fileID2= fopen('test18.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

fileID2= fopen('test19.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

fileID2= fopen('test20.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor1=B(1:3:length(B));
ycor1=B(2:3:length(B));
xcor1=mean(xcor1);
ycor1=mean(ycor1);
xcor = [xcor;xcor1];
ycor = [ycor;ycor1];

%Define variables.
dt=0.2;
HexAccel_noise_mag= 0.15;
Nb_particles=200;
u=0.00; %acceleration magnitude
% Initialize state transition matrix
A=[ 1 0 1 0;...     % [x  ]
    0 1 0 1;...     % [y]
    0 0 1 0;...     % [Vx  ]
    0 0 0 1 ];     % [Vy]
B=[(dt^2/2); (dt^2/2); dt; dt];
C=[1 0 0 0; 0 1 0 0];
Ez=[1 0; 0 1];
% Ex=[dt^4/4 0 dt^3/2 0;...
%     0 dt^4/4 0 dt^3/2; ...
%     dt^3/2 0 dt^2 0; ...
%     0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2;
Ex=[dt^4/4 0 dt^2/2 0;...
    0 dt^4/4 0 dt^2/2; ...
    dt^2/2 0 1 0; ...
    0 dt^2/2 0 1].*HexAccel_noise_mag^2;
P=Ex;


    x_est = zeros(4, 1);             % x_est=[x,y,Vx,Vy]'
    p_est = P;


Nb=1;

for i=1:1:floor(length(xcor)/Nb)
    m(i) = mean(xcor((i-1)*Nb+1:i*Nb));
    n(i) = mean(ycor((i-1)*Nb+1:i*Nb));
    z(1,i)=m(i);
    z(2,i)=n(i);
    x_est(1,1)=m(1);
    x_est(2,1)=n(1);
    x_est(3,1)=0.1;
    x_est(4,1)=0.1;
    x_fin(1)=m(1);
    y_fin=n(1);
end

%test 23
xval = 2.2;
ymin = 1.15;
ymax = 4.5;

hold on

%plot the real trajectory
xx=[xval,xval];
yy=[ymin,ymax];
plot(xx,yy)

%plot the four anchors
scatter(0,0,'s','filled')
scatter(3.75,0.05,'s','filled')
scatter(3.93,5.81,'s','filled')
scatter(0.05,5.76,'s','filled')

%plot the four people 
rectangle('Position',[1.25 4.5 0.4 0.25])
rectangle('Position',[2.6 4.5 0.4 0.25])
rectangle('Position',[1.25 1 0.4 0.25])
rectangle('Position',[2.6 1 0.4 0.25])


for i=1:1:floor(length(xcor)/Nb)-1
    % Predicted state and covariance
    x_prd(:,i) = A * x_est(:,i) + B*u;
    p_prd = A * p_est * A' + Ex;
    % Kalman gain
    K = p_prd*C'*inv(C*p_prd*C'+Ez);

    % Estimated state and covariance
    x_est(:,i+1) = x_prd(:,i) + K * (z(:,i+1) - C * x_prd(:,i));
    p_est = p_prd - K * C * p_prd;
    % Compute the estimated measurements
    y(:,i) = C * x_est(:,i);
    
    x_P= mvnrnd(x_fin(i),0.1*abs(x_est(3,i)),Nb_particles);
    y_P= mvnrnd(y_fin(i),0.1*abs(x_est(4,i)),Nb_particles);
    
%     x_P= mvnrnd(z(1,i),0.0017,Nb_particles);
%     y_P= mvnrnd(z(2,i),0.013,Nb_particles);
%     
        for j=1:1:Nb_particles

        %Weight(j) = exp(-(x_P(j)-y(1,i))^2/(2*0.0017^2) - (y_P(j)-y(2,i))^2/(2*0.013^2));
        dist(j)=(x_P(j)-m(i+1))^2 + (y_P(j)-n(i+1))^2;
        end
        
        [~,index]=min(dist);
        
        %if(sum(Weight)>0)
        %Final_Weight=Weight/sum(Weight);
        %end
        x_fin(i+1)= x_P(index);
        y_fin(i+1)= y_P(index);

      % scatter(x_fin,y_fin,'g','filled')
      %scatter(x_P,y_P,'g')
      scatter(y(1,i),y(2,i),'r','filled')
      scatter(m(i),n(i),'y','filled')
      scatter(x_fin(i),y_fin(i),'b','filled')
      %pause;
end

% figure
% hold on
% scatter(m,n,'b')
% scatter(y(1,:),y(2,:),'r','filled')
% scatter(x_fin,y_fin,'g','filled')
% hold off

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Kalmann filter
clear all
close all

fileID2= fopen('human12.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:3:length(B));
ycor=B(2:3:length(B));

%Define variables.
dt=0.2;
Nb_particles=50;
HexAccel_noise_mag=0.15;
% Ex=[dt^4/4 0 dt^3/2 0;...
%     0 dt^4/4 0 dt^3/2; ...
%     dt^3/2 0 dt^2 0; ...
%     0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2;
Ex=[dt^4/4 0 dt^2/2 0;...
    0 dt^4/4 0 dt^2/2; ...
    dt^2/2 0 1 0; ...
    0 dt^2/2 0 1].*HexAccel_noise_mag^2;
P=Ex;

% x_est=[x,y,Vx,Vy]'
x_est = zeros(4, 1);

for j=1:Nb_particles
    p_est(:,:,j) = P;
end
  
j=1;
speed = zeros(2,Nb_particles);
Nb=1;
bool0=1;

Weight_par0=0;

err0(1)=0;


for i=1:1:floor(length(xcor)/Nb)
    m(i) = mean(xcor((i-1)*Nb+1:i*Nb));
    n(i) = mean(ycor((i-1)*Nb+1:i*Nb));
    z(1,i)=m(i);
    z(2,i)=n(i);
    x_est(1,1)=m(1);
    x_est(2,1)=n(1);
    x_est(3,1)=0.0;
    x_est(4,1)=0.0;

end

x_P= mvnrnd(m(1),0.1,Nb_particles);
y_P= mvnrnd(n(1),0.1,Nb_particles);
x_fin(1)=mean(x_P);
y_fin(1)=mean(y_P);
%test 9
% xval = 2.15;
% ymin = 1;
% ymax = 5.7;

% test human
xval = 1.95;
ymin = 1.2;
ymax = 4.3;

%test 23
% xval = 2.2;
% ymin = 1.1;
% ymax = 5.5;

%test 30
% xval = 2.25;
% ymin = 1.2;
% ymax = 4.6;

%test 22
% xval = 2.15;
% ymin = 0.6;
% ymax = 5.3;


 figure
 hold on

%plot the real trajectory
xx=[xval,xval];
yy=[ymin,ymax];
%plot(xx,yy)

y_real=linspace(ymin,ymax,floor(length(xcor)/Nb));

%plot the four anchors
scatter(0,0,'s','filled')
scatter(3.75,0.05,'s','filled')
scatter(3.93,5.81,'s','filled')
scatter(0.05,5.76,'s','filled')

%plot the four people
rectangle('Position',[1.25 4.5 0.4 0.25])
rectangle('Position',[2.6 4.5 0.4 0.25])
rectangle('Position',[1.25 1 0.4 0.25])
rectangle('Position',[2.6 1 0.4 0.25])

for i=1:1:floor(length(xcor)/Nb)-1
    %kalman filter
    %[x_est(:,i+1),p_est] = kalm(m(i),n(i),x_est(:,i),p_est);
    
    %real points
    X_real=[xval;y_real(i+1)];
    err_raw(i+1)=my_distance(X_real,z(:,i+1),'L2');
    
    %Particle filter
    [x_fin(i+1),y_fin(i+1),x_P,y_P,p_est,speed] = particle_fil(x_fin(i),y_fin(i),p_est...
        ,Nb_particles,m(i+1),n(i+1),x_P,y_P,speed);
    
    X_F=[x_fin(i+1);y_fin(i+1)];
    
    err0(i+1)=my_distance(X_real,X_F,'L2');
    
    
    %scatter(xval,y_real(i),'g','filled')
    %scatter(x_est(1,i),x_est(2,i),'r','filled')
    scatter(m(i),n(i),'y','filled')
    scatter(x_fin(i),y_fin(i),'b','filled', 'handlevisibility', 'off')
    

    pause;
    %clf;
    %legend('True trajectory','Measured data','particle filter 1','particle filter 2');

end
%%

%Kalmann filter
clear all


fileID2= fopen('Alessio1.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:3:length(B));
ycor=B(2:3:length(B));

%Define variables.
dt=0.2;
Nb_particles=1500;
HexAccel_noise_mag=0.25;
Ex=[dt^4/4 0 dt^2/2 0;...
    0 dt^4/4 0 dt^2/2; ...
    dt^2/2 0 1 0; ...
    0 dt^2/2 0 1].*HexAccel_noise_mag^2;
P=Ex;

x_est = zeros(4, 1);
p_est = P;
err0(1)=0;
err1(1)=0;
figure
hold on
[X_real,Y_real]=draw_map();


for i=1:1:floor(length(xcor))
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

x_P= mvnrnd(x_fin,1,Nb_particles);
y_P= mvnrnd(y_fin,1,Nb_particles);


for i=1:1:79
    %kalman filter
    [x_est(:,i+1),p_est] = kalman(m(i+1),n(i+1),x_est(:,i),p_est);
    
    %real points
    X_=[X_real(i+1);Y_real(i+1)];
    err_raw(i+1)=my_distance(X_,z(:,i+1),'L2');
    err1(i+1)=my_distance(X_,[x_est(1,i+1);x_est(2,i+1)],'L2');
    
    %Particle filter
    [x_fin(i+1),y_fin(i+1),x_P,y_P] = particle_filter(x_fin(i),y_fin(i),x_est(:,i+1)...
        ,Nb_particles,m(i+1),n(i+1),x_P,y_P);
    
    X_F=[x_fin(i+1);y_fin(i+1)];
    
    err0(i+1)=my_distance(X_,X_F,'L2');
    
    
    scatter(X_real(i),Y_real(i),'g','filled', 'handlevisibility', 'on')
    %scatter(x_est(1,i),x_est(2,i),'r','filled')
    scatter(m(i),n(i),'y','filled', 'handlevisibility', 'on')
    scatter(x_fin(i),y_fin(i),'b','filled', 'handlevisibility', 'on')
    plot([x_fin(i) X_real(i)],[y_fin(i) Y_real(i)],'b','LineWidth',1,'handlevisibility','off')

    legend('Particles','real trajectory','measurement point','Bigest Cluster estimation');
    pause
    cla

end



%%
%Kalmann filter
clear all
figure

fileID2= fopen('humana8.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:3:length(B));
ycor=B(2:3:length(B));

%Define variables.
dt=0.2;
Nb_particles=150;
HexAccel_noise_mag=0.2;
% Ex=[dt^4/4 0 dt^3/2 0;...
%     0 dt^4/4 0 dt^3/2; ...
%     dt^3/2 0 dt^2 0; ...
%     0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2;
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

[X_real,Y_real]=draw_map(length(xcor));

for i=1:1:floor(length(xcor))
    m(i)=xcor(i);
    n(i)=ycor(i);
    z(1,i)=m(i);
    z(2,i)=n(i);
    %x_fin(1)=m(1);
    %y_fin(1)=n(1);
end


x_P= mvnrnd(m(1),0.5,Nb_particles);
y_P= mvnrnd(n(1),0.5,Nb_particles);
x_fin(1)=mean(x_P);
y_fin(1)=mean(y_P);

for i=1:1:length(xcor)-1
    
    %real points
    X_=[X_real(i+1);Y_real(i+1)];
    err_raw(i+1)=my_distance(X_,z(:,i+1),'L2');
    
    %Particle filter
    [x_fin(i+1),y_fin(i+1),x_P,y_P,p_est,speed] = particle_fil(p_est...
        ,Nb_particles,m(i+1),n(i+1),x_P,y_P,speed);
    
    X_F=[x_fin(i+1);y_fin( i+1)];
    
    err0(i+1)=my_distance(X_,X_F,'L2');
    
    
    %scatter(X_real(i),Y_real(i),'g','filled', 'handlevisibility', 'on')
    scatter(m(i),n(i),'r','filled', 'handlevisibility', 'off')
    scatter(x_fin(i),y_fin(i),'b','filled', 'handlevisibility', 'off')
    plot([x_fin(i) X_real(i)],[y_fin(i) Y_real(i)],'b','LineWidth',1,'handlevisibility','off')
    %plot([m(i) X_real(i)],[n(i) Y_real(i)],'r','LineWidth',1,'handlevisibility','on')
    drawnow
    pause(0.25)
    cla
    
    %legend('True trajectory','Measured data','particle filter 1','particle filter 2');
    
end
mean(err_raw)
mean(err0)

%%
%Kalmann filter
clear all
figure

class=2;
fileID2= fopen('alessiokio2.txt','r');
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
% Ex=[dt^4/4 0 dt^3/2 0;...
%     0 dt^4/4 0 dt^3/2; ...
%     dt^3/2 0 dt^2 0; ...
%     0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2;
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

[X_real,Y_real]=draw_map();

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
    z1(2,i)=m(2,i);
    %x_fin(1)=m(1);
    %y_fin(1)=n(1);
end


x_P= mvnrnd(0.5*(m(1,1)+m(2,1)),0.5,Nb_particles);
y_P= mvnrnd(0.5*(n(1,1)+n(2,1)),0.5,Nb_particles);
x_fin(1)=mean(x_P);
y_fin(1)=mean(y_P);

for i=1:1:length(xcor)-1
    
    %real points
    %X_=[X_real(i+1);Y_real(i+1)];
    %err_raw(i+1)=my_distance(X_,z(:,i+1),'L2');
    
%     dist_meas=my_distance(z(:,i+1),z1(:,i+1),'L2');
%     if(dist_meas>1.5)
%         dist_meas = my_distance(z(:,i+1),z(:,i),'L2');
%         dist_meas1= my_distance(z1(:,i+1),z(:,i),'L2');
%         if(dist_meas>dist_meas1)
%             m(1,i+1)=m(2,i+1);
%             n(1,i+1)=n(2,i+1);
%         else
%             m(2,i+1)=m(1,i+1);
%             n(2,i+1)=n(1,i+1);            
%         end
%     end
    
    %Particle filter
    [x_fin(i+1),y_fin(i+1),x_P,y_P,p_est,speed] = particle_fil_d(class,p_est...
        ,Nb_particles,m(:,i+1),n(:,i+1),x_P,y_P,speed);
      
    %X_F=[x_fin(i+1);y_fin( i+1)];
    
    %err0(i+1)=my_distance(X_,X_F,'L2');
    
    
    %scatter(X_real(i),Y_real(i),'g','filled', 'handlevisibility', 'on')
    scatter(m(1,i),n(1,i),'r','filled', 'handlevisibility', 'on')
    scatter(m(2,i),n(2,i),'r','filled', 'handlevisibility', 'on')
    quiver(x_fin(i),y_fin(i),cos(ori(i)-pi/2),sin(ori(i)-pi/2),'b','filled','handlevisibility','on')
    scatter(x_fin(i),y_fin(i),'b','filled', 'handlevisibility', 'off')
    %plot([x_fin(i) X_real(i)],[y_fin(i) Y_real(i)],'b','LineWidth',1,'handlevisibility','off')
    %plot([m(i) X_real(i)],[n(i) Y_real(i)],'r','LineWidth',1,'handlevisibility','on')
    drawnow
    pause(0.25)
    cla
    
    %legend('True trajectory','Measured data','particle filter 1','particle filter 2');
    
end
%%
%Kalmann filter
clear all
figure

fileID2= fopen('tagontag6.txt','r');
formatSpec = '%*s %f %*s %f %*s %f %*s';
B=fscanf(fileID2,formatSpec);
xcor=B(1:6:length(B));
ycor=B(2:6:length(B));
xcor1=B(5:6:length(B));
ycor1=B(5:6:length(B));

%Define variables.
dt=0.2;
Nb_particles=150;
HexAccel_noise_mag=0.2;
% Ex=[dt^4/4 0 dt^3/2 0;...
%     0 dt^4/4 0 dt^3/2; ...
%     dt^3/2 0 dt^2 0; ...
%     0 dt^3/2 0 dt^2].*HexAccel_noise_mag^2;
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

[X_real,Y_real]=draw_map();

for i=1:1:floor(length(xcor))
    m(i)=xcor(i);
    n(i)=ycor(i);
    m1(i)=xcor1(i);
    n1(i)=ycor1(i);
    z(1,i)=m(i);
    z(2,i)=n(i);
    %x_fin(1)=m(1);
    %y_fin(1)=n(1);
end


x_P= mvnrnd(m(1),0.5,Nb_particles);
y_P= mvnrnd(n(1),0.5,Nb_particles);
x_fin(1)=mean(x_P);
y_fin(1)=mean(y_P);

x_P1= mvnrnd(m1(1),0.5,Nb_particles);
y_P1= mvnrnd(n1(1),0.5,Nb_particles);
x_fin1(1)=mean(x_P1);
y_fin1(1)=mean(y_P1);

for i=1:1:length(xcor)
    
    %real points
    %X_=[X_real(i+1);Y_real(i+1)];
    %err_raw(i+1)=my_distance(X_,z(:,i+1),'L2');
    
    %Particle filter
    [x_fin(i+1),y_fin(i+1),x_P,y_P,p_est,speed] = particle_fil(p_est...
        ,Nb_particles,m(i+1),n(i+1),x_P,y_P,speed);
    
    [x_fin1(i+1),y_fin1(i+1),x_P1,y_P1,p_est1,speed1] = particle_fil(p_est1...
        ,Nb_particles,m1(i+1),n1(i+1),x_P1,y_P1,speed1);
    
    %X_F=[x_fin(i+1);y_fin( i+1)];
    
    %err0(i+1)=my_distance(X_,X_F,'L2');
    
    
    %scatter(X_real(i),Y_real(i),'g','filled', 'handlevisibility', 'on')
    scatter(m(i),n(i),'r','filled', 'handlevisibility', 'on')
    scatter(x_fin(i),y_fin(i),'b','filled', 'handlevisibility', 'off')
    scatter(x_fin1(i),y_fin1(i),'m','filled', 'handlevisibility', 'off')
    %plot([x_fin(i) X_real(i)],[y_fin(i) Y_real(i)],'b','LineWidth',1,'handlevisibility','off')
    %plot([m(i) X_real(i)],[n(i) Y_real(i)],'r','LineWidth',1,'handlevisibility','on')
    drawnow
    pause(0.25)
    cla
    
    %legend('True trajectory','Measured data','particle filter 1','particle filter 2');
    
end




