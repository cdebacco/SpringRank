% Code by Daniel Larremore
% Santa Fe Institute
% larremore@santafe.edu
% http://danlarremore.com

% evaluate the local accuracy of edge direction prediction

function a = localAccuracy_BTL(A,g)
m = sum(sum(A));
n = length(g);
y = 0;
for i=1:n
    for j=1:n
        p = g(i)/(g(i)+g(j)); % BTL probability
        if isnan(p)
            % do nothing
        else
            y = y + abs( A(i,j)-( A(i,j) + A(j,i) )*p );
        end
    end
end
a = 1-0.5*y/m;
