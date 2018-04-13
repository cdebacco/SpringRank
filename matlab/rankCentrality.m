% Rank Centrality
% Implemented by Dan Larremore, University of Colorado Boulder
% April 8, 2018
%
% Based on the manuscript
% Rank Centrality: Ranking from Pairwise Comparisons
% Sahand Negahban, Sewoong Oh, Devavrat Shah
% 2017
%
function [rc] = rankCentrality(A)
% In their text, a_ij = # of times j is preferred over i.
% In the SpringRank paper, we usually assume the opposite. 
% Here, we'll use the authors' direction, but note that whenever we call
% this code, we'll have to pass the transpose of A. 

% Note that there are no self-loops in this model, so we will check, 
% discard, and warn 
N = size(A,1);
if sum(diag(A)) > 0
%     fprintf('Warning: self-loops detected (and ignored)\n')
    A(1:N+1:end) = 0;
end


% see Eq 5 of https://arxiv.org/pdf/1209.1688.pdf
% We're going to regularize.
% They suggest epsilon = 1. 
% This seems extreme?

% Not listed in the paper, but this is important. We have to regularize the
% matrix A before we compute dmax. 
regularization = 1;
A = A+regularization;

% Find dmax
dout = sum(A,2);
dmax = max(dout);

% Eq 5

P = (1/dmax) * A./(A+transpose(A));

% But we need to make sure that the matrix remains stochastic by making the
% rows sum to 1. Without regularization, Eq 1 says P(i,i) = 1 - dout(i)/dmax;
% Instead, we're going to just do this "manually"

P(1:N+1:end) = 0;
for i=1:N
    P(i,i) = 1-sum(P(i,:));
end

[V,~] = eigs(transpose(P),1);

rc = V / sum(V);

end