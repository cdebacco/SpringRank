% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% y = globalAccuracy(A,s,b)
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
% y is the global accuracy (called \sigma_L in the paper)

function y = globalAccuracy(A,s,b)

% number of nodes
n = length(s);
% accumulate the log likelihood score elementwise
y = 0;
for i=1:n
    for j=1:n
        d = s(i) - s(j);
        p = (1+exp(-2*b*d))^(-1);
        if p==0 || p==1
            % do nothing
        else
            y = y + A(i,j)*log(p)+A(j,i)*log(1-p);
        end
    end
end
end