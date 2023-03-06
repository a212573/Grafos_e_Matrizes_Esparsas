a=(1/2)^(1/2); %alpha
A=[0   -1    0    0    0    1    0    0    0    0    0    0    0
   0    0    1    0    0    0    0    0    0    0    0    0    0
  -a    0    0    1    a    0    0    0    0    0    0    0    0
  -a    0   -1    0   -a    0    0    0    0    0    0    0    0  
   0    0    0   -1    0    0    0    1    0    0    0    0    0
   0    0    0    0    0    0   -1    0    0    0    0    0    0 
   0    0    0    0   -a   -1    0    0    a    1    0    0    0
   0    0    0    0    a    0    1    0    a    0    0    0    0
   0    0    0    0    0    0    0    0    0   -1    0    0    1
   0    0    0    0    0    0    0    0    0    0    1    0    0
   0    0    0    0    0    0    0   -1   -a    0    0    a    0
   0    0    0    0    0    0    0    0   -a    0   -1   -a    0
   0    0    0    0    0    0    0    0    0    0    0   -a    1];

b=[ 0
   10
    0
    0
    0
    0
    0
    15
    0
    20
    0
    0
    0];
P=[0    0    1    0    0    0    0    0    0    0    0    0    0
   1    0    0    0    0    0    0    0    0    0    0    0    0
   0    1    0    0    0    0    0    0    0    0    0    0    0
   0    0    0    1    0    0    0    0    0    0    0    0    0
   0    0    0    0    1    0    0    0    0    0    0    0    0
   0    0    0    0    0    0    1    0    0    0    0    0    0
   0    0    0    0    0    1    0    0    0    0    0    0    0
   0    0    0    0    0    0    0    0    0    0    1    0    0
   0    0    0    0    0    0    0    1    0    0    0    0    0
   0    0    0    0    0    0    0    0    1    0    0    0    0
   0    0    0    0    0    0    0    0    0    1    0    0    0
   0    0    0    0    0    0    0    0    0    0    0    1    0
   0    0    0    0    0    0    0    0    0    0    0    0    1];


n=size(A);
n=n(1);
PA=P*A;
if det(A)~=0
    U=PA;
    i=1; 
    for r=randn(0,n-1)
        i=i+1;
        k=i-1; 
        for r=randn(0,n-i+1)
            k=k+1;
            j=0;
            k3=U(k,i-1)/U(i-1,i-1);
            if U(i-1,i-1)==0 
                k3=0;
            end
            for r=randn(0,n)
                j=j+1;
                U(k,j)=U(k,j)-U(i-1,j)*k3;

            end
        end    
    end
    L=PA/U;
end


x_sup=1;%descrobrindo tamanho de ENV_sup e ENVlin_sup
j=1;
for r=randn(0,n-1)
    j=j+1;
    i=0;
    for r=randn(0,j-1)
        i=i+1;
        if U(i,j)~=0
            x_sup=x_sup+j-i;
            break
        end
    end                
end
ENV_sup=zeros(1,x_sup);%gerando ENV_sup e ENVlin_sup
ENVlin_sup=zeros(1,x_sup);
ie=0;%indice de ENV
il=0;%indice de ENVlin
j=1;
for r=randn(0,n-1)
    j=j+1;
    i=0;
    for r=randn(0,j-1)
        i=i+1;
        if U(i,j)~=0
            for r=randn(0,j-i)
                ie=ie+1;
                ENV_sup(1,ie)=U(i,j);
                il=il+1;
                ENVlin_sup(1,il)=i;
                i=i+1;
            end
            break
        end
    end                
end
ENVcol_sup=zeros(1,n+1);%gerando ENVcol_sup
ENVcol_sup(1,1)=1;
ENVcol_sup(1,2)=1;
ic=2;%indice de ENVcol_sup
j=1;
for r=randn(0,n-1)
    ic=ic+1;
    j=j+1;
    i=0;
    for r=randn(0,j-1)
        i=i+1;
        if U(i,j)~=0
            ENVcol_sup(1,ic)=ENVcol_sup(1,ic-1)+j-i;
            break
        end
        if i==(j-1) && U(i,j)==0
            ENVcol_sup(1,ic)=ENVcol_sup(1,ic-1);
        end 
    end                
end

DIAG_U=zeros(1,n);%gerando DIAG_U
i=0;
for r=randn(0,n)
    i=i+1;
    DIAG_U(i)=U(i,i);
end
DIAG_U

ENV_sup
ENVcol_sup
ENVlin_sup


Lt=L';%encontrando a transposta de L
x_inf=1;%descrobrindo tamanho de ENV_inf e ENVcol_inf
j=1;
for r=randn(0,n-1)
    j=j+1;
    i=0;
    for r=randn(0,j-1)
        i=i+1;
        if Lt(i,j)~=0
            x_inf=x_inf+j-i;
            break
        end
    end                
end
ENV_inf=zeros(1,x_inf);%gerando ENV_inf e ENVcol_inf
ENVcol_inf=zeros(1,x_inf);
ie=0;%indice de ENV
il=0;%indice de ENVcol
j=1;
for r=randn(0,n-1)
    j=j+1;
    i=0;
    for r=randn(0,j-1)
        i=i+1;
        if Lt(i,j)~=0
            for r=randn(0,j-i)
                ie=ie+1;
                ENV_inf(1,ie)=Lt(i,j);
                il=il+1;
                ENVcol_inf(1,il)=i;
                i=i+1;
            end
            break
        end
    end                
end
ENVlin_inf=zeros(1,n+1);%gerando ENVlin_inf
ENVlin_inf(1,1)=1;
ENVlin_inf(1,2)=1;
ic=2;%indice de ENVlin_inf
j=1;
for r=randn(0,n-1)
    ic=ic+1;
    j=j+1;
    i=0;
    for r=randn(0,j-1)
        i=i+1;
        if Lt(i,j)~=0
            ENVlin_inf(1,ic)=ENVlin_inf(1,ic-1)+j-i;
            break
        end
        if i==(j-1) && Lt(i,j)==0
            ENVlin_inf(1,ic)=ENVlin_inf(1,ic-1);
        end 
    end                
end

DIAG_L=zeros(1,n);%gerando DIAG_L
i=0;
for r=randn(0,n)
    i=i+1;
    DIAG_L(i)=L(i,i);
end
DIAG_L

ENV_inf
ENVlin_inf
ENVcol_inf

f= inv(L*U)*P*b;

spy(L);

figure(1)
spy(A)
title("Matriz A")
figure(2)
spy(P*A)
title("Matriz P*A")
figure(3)
spy(L)
title("Matriz L")
figure(4)
spy(U)
title("Matriz U")


