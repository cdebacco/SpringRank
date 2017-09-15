% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% H = springRankHamiltonian(s,A,mu)
%   INPUTS:
% s is a N-vector of node positions
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
% mu can be a scalar or a NxN matrix
%   OUTPUTS:
% H is the scalar spring energy of the system;
%   NOTE: assumes spring rest length of 1

function [H] = springRankHamiltonian(s,A,mu)
% assuming A is sparse, faster to go through entries of A.
[r,c,v] = find(A);
% preallocate container of hamiltonian values.
h = zeros(size(v));

% NOTE: i = r(n) and j = c(n)

% Probably could be faster
if length(mu)==1 % SCALAR spring constant
    for n = 1:length(v)
        h(n) = v(n) * (s(r(n))-s(c(n))-1)^2;
    end
    h = h*mu;
else % MATRIX of spring constants
    for n = 1:length(v)
        h(n) = v(n) * mu(r(n),c(n)) * (s(r(n))-s(c(n))-1)^2;
    end
end

% sum up, divide by 2, and return
H = sum(h)/2;