% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% r = colleyMatrix(A)
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
%   OUTPUTS:
% r is the Nx1 vector of the Colley Matrix ranks

function [r] = colleyMatrix(A)
%Aij = i beats j
%therefore out-degree = sum over j = wins
%therefore in-degree = sum over i = losses
A(1:max(size(A))+1:end) = 0;
wins = sum(A,2);
losses = sum(A,1)';
total = wins+losses;
matches = A+transpose(A);

n = size(A,1);
C = zeros(n,n);
b = zeros(n,1);
for i=1:n
    b(i) = 1 + (wins(i)-losses(i))/2;
    for j=1:n
        if i==j
            C(i,j) = 2+total(i);
        else
            C(i,j) = -matches(i,j);
        end
    end
end
r = C\b;