clear all

disp('EX 3.1-3.2');
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    n=10;
    a_i=1;
    b_i=2;
    A = full(gallery('tridiag',n,-a_i,4,-b_i));
    b=sum(A,2);

    tol=0.000001/2;

    CL=-tril(A, -1);
    %disp('CL');disp(CL);

    CU=-triu(A, 1);
    %disp('CU');disp(CU);


    I=eye(n);
    %disp(I);
    D=diag(diag(A));
    %disp('D');disp(D);

    D1=inv(D);
    %disp('D1'); disp(D1);

    L=D1*CL;
    %disp('L');disp(L);

    U=D1*CU;
    %disp('U');disp(U);


    B=L+U;
    %disp('B');  disp(B);
    %idiotimes tou pinaka B (Jacobi)
        disp('idiotimes');
        x=eig(B);     %disp(x);
        rB=max(abs(x));
    omega=2.0/(1.0+sqrt(1-rB*rB));
    
    disp('rB - fasmatikh aktina'); disp(rB);
    
    loop = 1;
    itcount = [];

for omega = 0.1:0.1:1.9
    
    x0=b;
    x1=x0;
    disp('omega'); disp(omega);
    itcount(loop)=0;
    maxits=50;

    while itcount(loop)<=maxits
       x0=x1;
       %x1=inv(I-omega*L)*((1-omega)*I+omega*U)*x0+omega*inv(I-omega*L)*D1*b;
       x1=inv(I-omega*U)*((1-omega)*I+omega*L)*x0+omega*inv(I-omega*U)*D1*b;
       nm=norm(x1-x0, Inf); 
       if nm<tol  
           iter=itcount(loop);
           disp('siglisi se'); disp(iter); disp('epanalipseis'); 
          break;   
       end
       itcount(loop)=itcount(loop)+1;
    end
    
    if nm>tol
        disp('oxi siglisi meta apo'); disp(maxits); disp('epanalipseis');
        itcount(loop) = itcount(loop) - 1;
    end
    %disp('x1'); disp(x1);
    disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~');
    loop = loop + 1;
    %clear all
end

omega = 0.1:0.1:1.9;
plot(omega,itcount)


