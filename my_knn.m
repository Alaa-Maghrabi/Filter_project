function [ y_est ] =  my_knn(X_train,  y_train, X_test, k, type)
%MY_KNN Implementation of the k-nearest neighbor algorithm
%   for classification.
%
%   input -----------------------------------------------------------------
%   
%       o X_train  : (N x M_train), a data set with M_test samples each being of dimension N.
%                           each column corresponds to a datapoint
%       o y_train  : (1 x M_train), a vector with labels y \in {0,1} corresponding to X_train.
%       o X_test   : (N x M_test), a data set with M_test samples each being of dimension N.
%                           each column corresponds to a datapoint
%       o k        : number of 'k' nearest neighbors
%       o type   : (string), type of distance {'L1','L2','LInf'}
%
%   output ----------------------------------------------------------------
%
%       o y_est   : (1 x M_test), a vector with estimated labels y \in {0,1} 
%                   corresponding to X_test.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sizes of datasets
 [~, M_train] = size(X_train);
 [N, M_test] = size(X_test);
 dist=zeros(1,M_train); 
 
for i=1:1:M_test
   
    %calculate the distances
   for j=1:1:M_train
       dist(1, j) = my_distance(X_train(:,j), X_test(:,i), type);
   end
   
   %sorting them
       [dis_sor,ind]=sort(dist);
       
   %choosing the closest ones
       y_close = y_train(ind(1:k));
       
       y_est(i) = mode(y_close,2);
       
    
end

end