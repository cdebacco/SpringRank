function [c] = katzCentrality(A,regularization)
% Katz Centrality implemented below

% number of nodes
N = max(size(A));

% Get the eigenvalues of A+regularization
% Take their absolute value. 
v = abs(eig(A+regularization));

% Set the damping factor as half the inverse largest eigenvalue
alph = 0.5/max(v);
% Set the scaling parameter b to 1
b = 1;

% x = a Ac + b
% x - a Ac = b1
% Ix - a Ac = b1
% (I-aA)c = b1
% c = (I-aA) \ b1

c = (eye(N)-alph*A) \ (b*ones(N,1)) ;


