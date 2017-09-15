% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% [A,P] = generativeModel(c,b,s)
%   INPUTS:
% c is the overall sparsity constant
% b is the inverse temperature (called beta in the paper)
% s is the Nx1 vector of planted node positions
%   OUTPUTS:
% A is a directed network adjacency matrix; A(i,j)=1 if i dominates j, e.g.
% P is a full NxN matrix; P(i,j) = expected number of edges from i to j

function [A,P] = generativeModel(c,b,s)

% number of vertices
N = length(s);
% preallocate P
P = zeros(N,N);
for i=1:N
    for j=1:N
        % compute expected number of edges from i to j
        P(i,j) = c*exp(-b/2*(s(i)-s(j)-1)^2);
    end
end
% draw poisson random numbers from P and return
A = poissrnd(P);