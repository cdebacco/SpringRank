% Code by Daniel Larremore
% Santa Fe Institute
% larremore@santafe.edu
% http://danlarremore.com

% evaluate the local accuracy of edge direction prediction

function y = globalAccuracy_BTL(A,g)
n = length(g);
y = 0;
for i=1:n
    for j=1:n
        p = g(i)/(g(i)+g(j)); % BTL probability
        if p==0 || p==1 || isnan(p)
            % do nothing
        else
            y = y + A(i,j)*log(p)+A(j,i)*log(1-p);
        end
    end
end

