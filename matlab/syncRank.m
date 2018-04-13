function [sy] = syncRank(A)

N = size(A,2);

% 1. Form C
% Whenever Aij > Aji we set Cij = 1
% Whenever Aij < Aji we set Cij = -1
% Else, Cij = 0;
% This means that Cij = sign(Aij - Aji)
C = sign(A-transpose(A));

% 2. Form Theta
T = pi*C/(N-1);

% 3. Form H
H = spalloc(N,N,nnz(T));
H(T~=0) = exp(1i*T(T~=0));

% 4. Form Dinv
Dinv = diag(1./sum(abs(H)));

% 5. Form fancyH
fancyH = Dinv*H;

% 6. Leading eigenvector of fancyH
[V,~] = eigs(fancyH,1);

% 7. Get angles in complex plane.
angles = angle(V);

% 8. Get order from angles
[~,idx] = sort(angles,'ascend');
sy(idx) = 1:N;

% 10. Choose the rank permutation that minimizes violations.
viols = zeros(N,1);
for ii=1:N
    sy_perm = mod(sy + ii - 2,N) + 1;
    idx_perm(sy_perm) = 1:N;
    viols(ii) = sum(sum(triu(A(idx_perm,idx_perm))));
end
best = find(viols==min(viols));

sy = (mod(sy + best(1) -2,N)+1)';
