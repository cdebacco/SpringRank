% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% [p,H0,H] = pvalueNullModel(A,n_repetitions)
%
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
% n_repetitions is the number of randomizations that you would like to use
%   to calculate the p-value. Higher numbers mean a better estimate of
%   probability that the p-value is meant to represent. Lower numbers may 
%   be required for very large networks whose randomizations are expensive
%   or slow.
%
%   OUTPUTS:
% p is the p-value described in the paper for the probability that a
%   network A whose edge directions are randomized would have a lower
%   ground-state energy than the original network A
% H0 is the ground state energy of A
% H is the vector of ground state energies associated with randomizations 
%   of the directions of A.

function [p,H0,H] = pvalueNullModel(A,n_repetitions)

% First determine what kind of matrix we have. Integer or non-integer?
[~,~,v] = find(A);
if sum(mod(v,1)==0)==length(v)
    isInteger = 1;
    Abar = A+A';
else
    isInteger = 0;
end

% Preallocate
H = zeros(n_repetitions,1);
H0 = springRankHamiltonian(springRank(A),A,1);

% Iterate over repetitions
for n = 1:n_repetitions
    % Two different randomization schemes, depending on whether or not the
    % network is integer or scalar.
    if isInteger
        B = randomEdgeDirectionsInt(Abar);
    else
        B = randomEdgeDirectionsScalar(A);
    end
    H(n) = springRankHamiltonian(springRank(B),B,1);
end

p = sum(H<H0)/n_repetitions;

end

function [B] = randomEdgeDirectionsInt(Abar)
N = max(size(Abar));
[r,c,v] = find(triu(Abar,1));
up = binornd(v,0.5);
down = v - up;
B = sparse([r;c],[c;r],[up;down],N,N);
end

function [B] = randomEdgeDirectionsScalar(A)
N = max(size(A));
[r,c,v] = find(A);
for i=1:length(v)
    if rand < 0.5
        temp = r(i);
        r(i) = c(i);
        c(i) = temp;
    end
end
B = sparse(r,c,v,N,N);
end