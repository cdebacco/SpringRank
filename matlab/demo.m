% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% This code will:
% 1. Generate a network via generativeModel.m
% 2. Use methods to find rankings from the resulting network
%       springRank.m
%       mvr.m
%       colleyMatrix.m
%       davidScore.m
%       btl.m
%       eigenvectorCentrality.m
% 

%% 1. Create a network

% choose preferred degree
deg = 5;
% choose inverse temperature beta
bet = 0.2;
% choose number of nodes N
N = 102;

% create planted ranks
s0 = [normrnd(4,sqrt(1),N/3,1);...
    normrnd(0,sqrt(0.5),N/3,1);...
    normrnd(-4,sqrt(2),N/3,1)];

% precompute spring energy matrix to keep degree fixed below
spr = zeros(N,N);
for i=1:N
    for j=1:N
        spr(i,j) = (s0(i) - s0(j) -1)^2;
    end
end
spr(1:N+1:end) = 0;
% expected degrees at inverse temperature B(i)
M = summ(exp(-bet/2*spr));
% adjust c to give expected degrees = k*n;
c = deg*N/M;
% create network using generative model
A = generativeModel(c,bet,s0);

% find number of components and analyze only the largest component
% THIS CODE CALLS networkComponents.m by Dan Larremore
% Included in the SpringRank github repo
% or mathworks.com/matlabcentral/fileexchange/42040-find-network-components
[nComponents,sizes,members] = networkComponents(A);
if nComponents > 1
    % extract largest component
    A = A(members{1},members{1});
    s0 = s0(members{1});
end

%remove self loops
A(1:N+1:end) = 0; 

fprintf('Network created.\n');
%% 2. Analyze the network
fprintf('Running SpringRank...\t\t')
tic
s_spr = springRank(A);
fprintf('%f seconds.\n',toc)

fprintf('Running BTL...\t\t\t')
tic
s_btl = btl(A,1e-3);
fprintf('%f seconds.\n',toc)

fprintf('Running Colley Matrix...\t')
tic
s_col = colleyMatrix(A);
fprintf('%f seconds.\n',toc)

fprintf('Running David Score...\t\t')
tic
s_dav = davidScore(A);
fprintf('%f seconds.\n',toc)

fprintf('Running Katz centrality...\t')
tic
s_eig = eigenvectorCentrality(A);
fprintf('%f seconds.\n',toc)

fprintf('Running MVR...\t\t\t')
tic
s_mvr = mvr(A,5);
fprintf('%f seconds.\n',toc)

fprintf('Running PageRank...\t\t')
tic
s_pag = pageRank(A,0.8,1e-12);
fprintf('%f seconds.\n',toc)