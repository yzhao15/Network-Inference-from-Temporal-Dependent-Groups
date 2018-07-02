function [theta,A,u,rho,alpha,beta,gamma,S]=EM(G)


T=size(G,1);
n=size(G,2);

alpha=0;
beta=0;
gamma=0;

rho=sum(G);
rho=rho/sum(rho);
u=log(rho)+2;

A=zeros(n);

for i=1:(n-1)
    for j=(i+1):n
        A(i,j)=(sum(G(:,i).*G(:,j))+sum(G(:,i).*G(:,j)))/(sum(G(:,i))+sum(G(:,j)));
    end
end

A=A+A';
theta=zeros(n);

for i=1:n
    A(i,i)=1;
    theta(i,i)=inf;
end

for i=1:n
    for j=1:n
        if ((A(i,j)~=0)&&(A(i,j)~=1))
            theta(i,j)=log(A(i,j)/(1-A(i,j)));
        end
        if (A(i,j)==0)
            theta(i,j)=-inf;
        end
        if (A(i,j)==1)
            theta(i,j)=inf;
        end
        
    end
end

diff=1;

while (diff>1e-3)
    [S,B]=E_step(G,theta,u,alpha,beta,gamma);
    [theta1,A1,u1,rho1,alpha1,beta1,gamma1]=M_step(G,S,B);
    

    diff_vec=[max(max(abs(A1-A))) abs(rho-rho1) abs(alpha-alpha1) abs(beta-beta1) abs(gamma-gamma1)];
    diff=max(diff_vec);
    theta=theta1;
    u=u1;
    alpha=alpha1;
    beta=beta1;
    gamma=gamma1;
    A=A1;
    rho=rho1;
    
end


