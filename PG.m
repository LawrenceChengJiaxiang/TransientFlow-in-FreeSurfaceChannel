function [X] = PG(A,B)
m = length(B);
for k = 1:m-1
    for i = k+1:m 
        g = A(i,k)/A(k,k);
        B(i,1) = B(i,1) - g*B(k,1);
        for j = k:m
            A(i,j) = A(i,j) - g*A(k,j);
        end
    end
end
X(1,m) = B(m,1)/A(m,m);
for i = m-1:-1:1
    sum = 0;
    for j = i+1:m
        sum = sum + A(i,j)*X(j);
    end
    X(1,i) = (B(i,1)-sum)/A(i,i);
end