% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% g = btl(A,tol)
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
% tol is the accuracy tolerance desired for successive iterations
%   OUTPUTS:
% s is the Nx1 vector of Davids Score
%   NOTE: implementation of a regularized version (for dangling node)
%   version of the algorithm presented in 
% Hunter DR (2004) MM algorithms for generalized Bradley-Terry models. 
% Annals of Statistics pp. 384?406

function [g] = btl(A,tol)

A(1:max(size(A))+1:end) = 0;
N = size(A,1);
g = rand(1,N); % random initial guesss
wins = sum(A,2);
matches = A+transpose(A);
totalMatches = sum(matches);
g_prev = rand(1,N);
eps=1e-6;
while norm(g-g_prev) > tol
    g_prev = g;
    for i=1:N
        if totalMatches(i)>0
            q = matches(i,:)./(g_prev(i)+g_prev);
            q(i) = [];
            g(i) = (wins(i)+eps)/sum(q);
        else
            g(i) = 0;
        end
    end
    g = g/sum(g);
end
g = transpose(g);
end