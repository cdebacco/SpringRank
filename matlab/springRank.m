% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% s = springRank(A)
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j. 
%   A(i,j) = # of times that j endorsed i.
%   OUTPUTS:
% s is a N-vector of node positions according to SpringRank

function [s] = springRank(A)

% Input check for sparsity
if issparse(A)
    % Great. We like sparse.
else
    A = sparse(A);
    % Still works, but slower
    warning('Input matrix not sparse; much faster if A is a sparse matrix.')
end

% Number of vertices
N = size(A,2);
% out-degree and in-degree VECTORS
dout = sum(A,2);
din = sum(A,1);
% out-degree and in-degree DIAGONAL MATRICES
Dout = diag(dout(1:N-1));
Din = diag(din(1:N-1));
% out-degree and in-degree for vertex N
dNout = dout(N);
dNin = din(N);

B = Dout + Din - A(1:N-1,1:N-1) - transpose(A(1:N-1,1:N-1))...
    - repmat(A(N,1:N-1),N-1,1) - repmat(transpose(A(1:N-1,N)),N-1,1);
b = diag(Dout)-diag(Din)+dNout-dNin;

% Sparse solve. Use [t,~] to suppress warnings. 
[t,~] = bicgstab(B,b,1e-12,200);

% ranks
s = [t;0];

% adjust mean of each component to be 0
[nComponents,~,members] = networkComponents(A);
for n = 1:nComponents
    s(members{n}) = s(members{n})-mean(s(members{n}));
end
