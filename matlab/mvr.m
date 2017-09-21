% SpringRank
% CODE  ->  https://github.com/cdebacco/SpringRank
% PAPER ->  http://danlarremore.com/pdf/SpringRank_2017_PrePrint.pdf
% Code by Daniel Larremore
% University of Colorado at Boulder
% BioFrontiers Institute & Dept of Computer Science
% daniel.larremore@colorado.edu
% http://danlarremore.com
%
% [order,violations,A] = mvr(A)
%   INPUTS:
% A is a NxN matrix representing a directed network
%   A can be weighted (integer or non-integer)
%   A(i,j) = # of dominance interactions by i toward j.
%   A(i,j) = # of times that j endorsed i.
% n_samples is an integer number of independent replicates of the MVR MCMC
% search procedure.
%   OUTPUTS:
% best_ranks is a vector of ranks. ONE IS BEST. N IS WORST
% best_violations is the number of violations
% best_A is the reordered matrix whose lower triangle contains min. viols.

function [best_ranks,best_violations,best_A] = mvr(A,n_samples)

best_violations = size(A,1)^2;

for n = 1:n_samples
    [ranks,violations,A] = mvr_single(A);
    if violations < best_violations
        best_violations = violations;
        best_ranks = ranks;
        best_A = A;
    end
end

end

function [ranks,violations,A] = mvr_single(A)
violations = compute_violations(A);

N = size(A,1);
%order = shuffle(1:N);
order =1:N;
A(:,:) = A(order,:);
A(:,:) = A(:,order);

step = 1;
fails = 0;
% fprintf('Initial\t\t\t\tviolations\t%i\n',violations);
hist_viols(step) = violations;
hist_viols_backup(step) = violations;
hist_fails(step) = fails;

% RANDOM STEPS - Randomly swap till N^2 failures in a row.
while 1
    i = randi(N); % choose random node
    j = randi(N); % choose second random node.
    while j==i % make sure different
        i = randi(N);
        j = randi(N);
    end
    dx = compute_violations_change(A,i,j);
    if dx < 0
        order([i,j]) = order([j,i]);
        A([i,j],:) = A([j,i],:);
        A(:,[i,j]) = A(:,[j,i]);
        step = step+1;
        hist_swaps(step,:) = [i,j];
        hist_fails(step) = fails;
        hist_viols(step) = hist_viols(step-1)+dx;
        violations = compute_violations(A);
        hist_viols_backup(step) = violations;
        %         fprintf('swap %i ~ %i \t --> %i\tviolations\t%i\t%i\n',i,j,dx,violations,fails)
        fails = 0;
    else
        fails = fails+1;
    end
    if fails == N^2
        A(1,:);
        %         fprintf('----- Too much fails -----\n');
        break
    end
end

% DETERMINISTIC STEPS - Find any local steps deterministically by search.
while 1
    dxbest = 0;
    for i=1:N-1
        for j=i+1:N
            dx = compute_violations_change(A,i,j);
            if dx < dxbest
                bestSwap = [i,j];
                dxbest = dx;
            end
        end
    end
    if dxbest==0
        %         fprintf('---- no improvement, exiting ----\n');
        [~,ranks] = sort(order);
        return;
    end
    i = bestSwap(1);
    j = bestSwap(2);
    order([i,j]) = order([j,i]);
    %     before = compute_violations(A);
    A([i,j],:) = A([j,i],:);
    A(:,[i,j]) = A(:,[j,i]);
    %     after = compute_violations(A);
    step = step+1;
    hist_swaps(step,:) = [i,j];
    hist_viols(step) = hist_viols(step-1)+dxbest;
    violations = compute_violations(A);
    hist_viols_backup(step) = violations;
    %     fprintf('swap %i ~ %i \t --> %i\tviolations\t%i\n',i,j,dxbest,violations)
end

end

function dx = compute_violations_change(A,ii,jj)
% Let's arbitrarily choose i to fall (larger rank number) and j to rise
% (smaller rank number).
i = min(ii,jj);
j = max(ii,jj);
dx= full(-sum(A(j,i:j-1)) ...
    +sum(A(i,i+1:j)) ...
    -sum(A(i+1:j-1,i)) ...
    +sum(A(i+1:j-1,j)));
end

function x = compute_violations(B)
x = full(sum(sum(tril(B,-1))));
end
