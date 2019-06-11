function [x_final,y_final,x_P,y_P] = particle_filter(x_fin,y_fin,x_est,Nb_particles,m,n,x_P,y_P)
%Particle Filter

persistent number;

if (isempty(number)) %It does not already exist
number = 1;
end

number=number + 1;

scatter(x_P,y_P,'k','.');

for j=1:1:Nb_particles
    x_P(j) = x_P(j) + x_est(3,1);
    %x_P(j) =  mvnrnd(x_P(j),0.1*abs(x_est(3,1)),1);
    y_P(j) = y_P(j) + x_est(4,1);
    %y_P(j) =  mvnrnd(y_P(j),0.1*abs(x_est(4,1)),1);
    %weight(j) = exp(-((x_P(j)-m)^2)/0.15 - ((y_P(j)-n)^2)/0.15);
    weight(j) = exp(-((x_P(j)-m)^2)/0.15 - ((y_P(j)-n)^2)/0.15)+exp(-((x_P(j)-x_est(1,1))^2)/0.15 - ((y_P(j)-x_est(2,1))^2)/0.15);
end

weight=weight/sum(weight);

x_final= weight*x_P;
y_final= weight*y_P;

%u(1)=rand/Nb_particles;
%c(1)=weight(1);

for i = 1 : Nb_particles%floor(Nb_particles*(1-rand))
    %c(i)= c(i-1) + weight(i);
    x_P(i) = x_P(find(rand <= cumsum(weight),1)) + 0.001*(rand-0.5);
    y_P(i) = y_P(find(rand <= cumsum(weight),1)) + 0.001*(rand-0.5);
end

if number==100
    x_P= mvnrnd(x_final,0.5,Nb_particles);
    y_P= mvnrnd(y_final,0.5,Nb_particles);
    weight=1/Nb_particles*ones(1,Nb_particles);
    number =0;
end

%fit_gaussian(x_P,y_P);

% i=1;
% for j=1:Nb_particles
%    u(j)=u(1) + (j-1)/Nb_particles;
%    while u(j)>c(i)
%        i=i+1;
%    end
%    x_P(j)=x_P(i);
%    y_P(j)=y_P(i);
%    weight(j)=1/Nb_particles;
% end

end