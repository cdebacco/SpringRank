% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% c = eigenvectorCentrality(A)
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
%   OUTPUTS:
% c is the Katz centrality with a small regularization for dangling nodes

function [c] = eigenvectorCentrality(A)
% Katz Centrality implemented below
err = 1e-6; %regularization for trees
v = abs(eig(A+err));
a = 0.5/max(v);
b = 1;
N = max(size(A));
c = b* ( (eye(N)-a*A)\ones(N,1) );
c = c/sum(c);