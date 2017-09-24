% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% b = betaLocal(A,s)
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
% s is the Nx1 vector of node positions (ranks)
%   OUTPUTS:
% b is the optimal inverse temperature (beta) under the LOCAL accuracy,
%   which we call \sigma_a in the paper.

function [b] = betaLocal(A,s)

global M r
M = A;
r = s;
b = fminbnd(@negacc,1e-6,1000);
end

function [a] = negacc(b)
global M r
m = sum(sum(M));
n = length(r);
y = 0;
for i=1:n
    for j=1:n
        d = r(i) - r(j);
        y = y + abs( M(i,j)- (M(i,j) + M(j,i))*((1+exp(-2*b*d))^(-1)) );
    end
end
a = y/m-1;
end