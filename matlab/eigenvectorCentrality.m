function [V] = eigenvectorCentrality(A,regularization)
% Eigenvector Centrality implemented below

% Let's get the Perron Frobenius eigenvector of A
[V,~] = eigs(A+regularization,1);
