% Parameter M adjacency matrix where M_i,j represents the link from 'j' to 'i', such that for all 'j'
%     sum(i, M_i,j) = 1
% Parameter d damping factor
% Parameter v_quadratic_error quadratic error for v
% Return v, a vector of ranks such that v_i is the i-th rank from [0, 1]

function v = pageRank(A, d, v_quadratic_error)
N = max(size(A)); % N is equal to either dimension of M and the number of documents
M = zeros(N,N);
for j=1:N
    if sum(A(:,j)) > 0
        M(:,j)=A(:,j)/sum(A(:,j));
    else
        M(:,j) = 1/N;
    end
end
% v = rand(N, 1);
v = ones(N,1);
v = v ./ norm(v, 1);   % This is now L1, not L2
last_v = ones(N, 1) * inf;
M_hat = (d .* M) + (((1 - d) / N) .* ones(N, N));

while(norm(v - last_v, 2) > v_quadratic_error)
	last_v = v;
	v = M_hat * v;
    v = v/norm(v,1);
    % removed the L2 norm of the iterated PR
end