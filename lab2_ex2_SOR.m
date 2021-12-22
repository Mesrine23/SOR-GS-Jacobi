clear all
disp('SOR');
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~');
for n=[10 100 1000]    
        disp('n= '); disp(n);
        a=1;
        b=2;
        A = full(gallery('tridiag',n,-a,4,-b));
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

        disp('rB - fasmatikh aktina'); disp(rB);
        omega=2.0/(1.0+sqrt(1-rB*rB));


    %for omega = 0.1:0.1:1.9

        x0=b;
        x1=x0;
        disp('omega'); disp(omega);
        itcount=0;
        maxits=50;

        while itcount<=maxits
           x0=x1;
           %x1 = inv(I-omega*L)*((1-omega)*I+omega*U)*x0 + omega*inv(I-omega*L)*D1*b;
           x1=inv(I-omega*U)*((1-omega)*I+omega*L)*x0+omega*inv(I-omega*U)*D1*b;
           nm=norm(x1-x0, Inf); 
           if nm<tol  
               iter=itcount;
               disp('siglisi se'); disp(iter); disp('epanalipseis'); 
               %disp('x1'); disp(x1);
              break;   
           end
           itcount=itcount+1;
        end
        if nm>tol
            disp('oxi siglisi meta apo'); disp(maxits); disp('epanalipseis');
        end
        %disp('x1'); disp(x1);
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~');

        %clear all
    %end
    clear all
end