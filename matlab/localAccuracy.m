% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% a = localAccuracy(A,s,b)
%
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
% s is the Nx1 vector of node positions (ranks)
% b is the inverse temperature (called beta in the paper)
%
%   OUTPUTS:
% a is the local accuracy (called \sigma_a in the paper)

function a = localAccuracy(A,s,b)
% total edges
m = sum(sum(A));

% number of vertices
n = length(s);

% accumulate accuracy of predictions
y = 0;
for i=1:n
    for j=1:n
        d = s(i) - s(j);
        p = (1+exp(-2*b*d))^(-1);
        y = y + abs( A(i,j)-( A(i,j) + A(j,i) )*p );
    end
end

% cleanup
a = 1-0.5*y/m;