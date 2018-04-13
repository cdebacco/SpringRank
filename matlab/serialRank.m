function [serr] = serialRank(A)

% In serialRank, C(i,j) = 1 if j was preferred over i. 
% The way that we construct the adj matrices in SpringRank, it's the opp.
% Therefore, pass in A transpose and implement serialRank faithfullly from
% the formulas.


% According to the SerialRank tutorial, we have
% If Aij = Aji = 0, then Qij = Qji = 1/2
% Else, Qij = 1/mij * sum (cij + 1)/2 where the sum is over the mij outcomes
% Noting that mij from SerialRank is Aij + Aji,
% And noting that the sum is (Aji - Aji + Aji + Aij)
% Else, Qij = ( 1 / 2(Aij+Aji) ) * (A_ij - A_ji + A_ij + A_ji)
% i.e.  Qij = ( 1 / 2(Aij+Aji) ) * 2A_ij
% i.e.  Qij = ( 1 / (Aij+Aji) ) * A_ij
% i.e.  Qij =  A_ij / (Aij+Aji)

N = size(A,1);
Q = A;
M = A+A';
Q(A~=0) = A(A~=0)./M(A~=0);
Q(M==0) = 1/2;

% Sij = sum over k, where B(i,k) and B(j,k) are nonzero, the quantity
% (1 - abs(Q(i,k)-Q(j,k)/2)
% and (sum) the quantity 1/2 otherwise.
% Hmmm. I'm not sure how to do this without a for loop... Probably there's
% a way, but I don't feel like optimizing someone else's alg. :/
% Also, S is dense and symmetric, I think.

S = zeros(N);
for i=1:N
    for j=i+1:N
        % Form the indicator Bik * Bjk == 0
        indicator = M(i,:).*M(j,:);
        S(i,j) = S(i,j) + sum(indicator==0)/2;
        nonzero = find(indicator>0);
        if ~isempty(nonzero)
            S(i,j) = S(i,j) + length(nonzero) - sum(abs(Q(i,nonzero) - Q(j,nonzero)))/2;
        end
    end
end
S = S+S';
L = diag(sum(S)) - S;
[V,~] = eigs(L,2,'smallestabs');
serr = V(:,2);

% [~,idx] = sort(W,'descend');
% serr(idx) = 1:N;
% serr = serr';