function [A,rho,S] = HM_EM(G)

T=size(G,1);
n=size(G,2);

rho=sum(G);
rho=rho/sum(rho);

A=zeros(n);

for i=1:(n-1)
    for j=(i+1):n
        A(i,j)=(sum(G(:,i).*G(:,j))+sum(G(:,i).*G(:,j)))/(sum(G(:,i))+sum(G(:,j)));
    end
end

A=A+A';

for i=1:n
    A(i,i)=1;
end

diff=1;

while (diff>1e-6)
    
    S=zeros(T,n);
    for t=1:T
        for i=1:n
            if G(t,i)==1
                S(t,i)=rho(i)*prod((A(i,:).^G(t,:)).*(1-A(i,:)).^(1-G(t,:)));
            end
        end
        S(t,:)=S(t,:)/sum(S(t,:));
    end
    rho1=sum(S)/T;
    A1=zeros(n);
    for i=1:(n-1)
        for j=(i+1):n
            if (sum(S(:,i))+sum(S(:,j)))>0
                A1(i,j)=(sum(S(:,i).*G(:,j))+sum(S(:,j).*G(:,i)))/(sum(S(:,i))+sum(S(:,j)));
            end
        end
    end
    A1=A1+A1';
    for i=1:n
        A1(i,i)=1;
    end
    
    diff_vec=[max(max(abs(A1-A))) abs(rho-rho1)];
    diff=max(diff_vec);
    A=A1;
    rho=rho1;
    
end

