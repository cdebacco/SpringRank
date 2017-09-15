% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% s = davidScore(A)
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
%   OUTPUTS:
% s is the Nx1 vector of Davids Score

function [s] = davidScore(A)

P = A./(A+A'); % Pij = Aij / (Aij + Aji)
P(isnan(P)) = 0;
P(1:size(P,1)+1:end) = 0; % ensure there are no entries on the diagonal
w = sum(P,2);
l = sum(transpose(P),2);
w2 = P*w;
l2 = transpose(P)*l;
s = w+w2-l-l2;
