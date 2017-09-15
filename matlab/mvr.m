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
%   OUTPUTS:
% order is a N-vector; MVR ranking of node indices; order(1) is best;
% violations is the number of violations
% A is the reordered matrix whose lower triangle contains min. viols.

function [order,violations,A] = mvr(A)

violations = compute_violations(A,reps);
    
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
        order;
%         fprintf('----- Too much fails -----\n');
        break
    end
end

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
