clear all
disp('Gauss-Seidel');
n=100;
a_i=1;
b_i=2;
A = full(gallery('tridiag',n,-a_i,4,-b_i));
b=sum(A,2);

x = b;
err = inf;
tol=0.000001/2;
max_it = 50;
flag=0;
itcount=0;

while(err>tol)
    if(itcount>max_it)
        disp(['Oxi sigklisi meta apo ' num2str(max_it)]);
        flag=1;
        break;
    end
    xold = x;
    for i=1:size(A,1)
        sum=0;
        for j=1:i-1
            sum = sum + A(i,j)*x(j);
        end
        for j=i+1:size(A,1)
            sum = sum + A(i,j)*xold(j);
        end
        x(i) = (1/A(i,i)) * (b(i)-sum);
    end
    itcount = itcount+1;
    err = abs(xold-x(i));
end

if(flag==0)
    disp(['Found solution after ' num2str(itcount) ' iterations']);
    %disp(x);
end
