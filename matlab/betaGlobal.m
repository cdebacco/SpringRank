% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% b = betaGlobal(A,s,varargin)
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
% s is the Nx1 vector of node positions (ranks)
%   OUTPUTS:
% b is the optimal inverse temperature (beta) under the GLOBAL accuracy,
%   which we call \sigma_\ell in the paper.

function [b] = betaGlobal(A,s,varargin)
global M r
M = A;
r = s;
b = fzero(@f,0.1);
end

function [y] = f(b)
global M r
n = length(r);
y = 0;
for i=1:n
    for j=1:n
        d = r(i) - r(j);
        pij = (1+exp(-2*b*d))^(-1);
        y = y + d*(M(i,j) - (M(i,j)+M(j,i))*pij);
    end
end
end