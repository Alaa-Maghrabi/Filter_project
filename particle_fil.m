function [x_final,y_final,x_P,y_P,p_est,speed] = particle_fil(p_est,Nb_particles,m,n,x_P,y_P,speed)
%Particle Filter

scatter(x_P,y_P,'k','.');

for j=1:1:Nb_particles
    XX=[x_P(j);y_P(j);speed(1,j);speed(2,j)];
    [XX,p_est(:,:,j)] = kalm(m,n,XX,p_est(:,:,j));
    x_P(j)=XX(1,1);
    y_P(j)=XX(2,1);
    speed(1,j)=XX(3,1);
    speed(2,j)=XX(4,1);
    weight(j) = exp(-(x_P(j)-m)^2/0.15 - (y_P(j)-n)^2/0.15);
end

weight=weight/sum(weight);

x_final= weight*x_P;
y_final= weight*y_P;

for i = 1 : Nb_particles
    x_P(i) = x_P(find(rand <= cumsum(weight),1)) + 0.001*rand;
    y_P(i) = y_P(find(rand <= cumsum(weight),1)) + 0.001*rand;
end

end