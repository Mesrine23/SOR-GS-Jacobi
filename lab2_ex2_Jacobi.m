clear all
disp('Jacobi');
n=10;
a_i=1;
b_i=2;
A = full(gallery('tridiag',n,-a_i,4,-b_i));
b=sum(A,2);


tol=0.000001/2;

[M,N] = size(A);

for m=1:M
    row = abs(A(m,:));
    d = sum(row) - row(m);
    if row(m) <= d
        error('[A] is not diagonally dominant, while it should be.');
    end
end

x = b;

D = diag(diag(A));

itcount = 0;
max_it = 50;
err = inf;
flag=0;

while err>tol
    %update x
    if(itcount>max_it)
        disp(['Oxi sigklisi meta apo ' num2str(max_it)]);
        flag=1;
        break;
    end
    dx = D\(b-A*x);
    x = x + dx;
    
    itcount = itcount + 1;
    %compute error
    err = max(abs(dx./x));
end

%display answer
if(flag==0)
    disp(['After ' num2str(itcount) ' iterations...']);
    %disp(x);
end

