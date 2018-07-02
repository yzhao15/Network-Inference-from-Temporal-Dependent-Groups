function [theta_hat,A_hat,u_hat,rho_hat,alpha_hat,beta_hat,gamma_hat]=M_step(G,S,B)

T=size(G,1);
n=size(G,2);
D1=zeros(n);
D2=zeros(n);
D3=zeros(n);
D4=zeros(n);
D5=zeros(n);
D6=zeros(n);


for i=1:n
    for j=1:n
        
        if (i~=j)
            D1(i,j)=S(1,i)*G(1,j)+sum(S(2:T,i).*(1-G(1:T-1,i)).*G(2:T,j));
            D2(i,j)=S(1,i)*(1-G(1,j))+sum(S(2:T,i).*(1-G(1:T-1,i)).*(1-G(2:T,j)));
            D3(i,j)=sum(S(2:T,i).*G(1:T-1,i).*G(1:T-1,j).*G(2:T,j));
            D4(i,j)=sum(S(2:T,i).*G(1:T-1,i).*G(1:T-1,j).*(1-G(2:T,j)));
            D5(i,j)=sum(S(2:T,i).*G(1:T-1,i).*(1-G(1:T-1,j)).*G(2:T,j));
            D6(i,j)=sum(S(2:T,i).*G(1:T-1,i).*(1-G(1:T-1,j)).*(1-G(2:T,j)));
        end
    end
end

freq=zeros(n);
for i=1:n
    for j=1:n
        freq(i,j)=sum(G(:,j).*S(:,i))+sum(G(:,i).*S(:,j));
    end
end

A_hat=zeros(n);
ind1=zeros(n);
ind0=zeros(n);
for i=1:n
    ind1(i,i)=1;
end

ind0(freq==0)=1;

for i=1:n
    for j=1:n
        if ((i~=j)&&(freq(i,j)~=0))
            A_hat(i,j)=(D1(i,j)+D3(i,j)+D5(i,j)+D1(j,i)+D3(j,i)+D5(j,i))/(D1(i,j)+D2(i,j)+D3(i,j)+D4(i,j)+D5(i,j)+D6(i,j)+D1(j,i)+D2(j,i)+D3(j,i)+D4(j,i)+D5(j,i)+D6(j,i));
            
        end
    end
end

for i=1:n
    A_hat(i,i)=1;
end
theta_hat=zeros(n);

for i=1:n
    for j=1:n
        if ((A_hat(i,j)>=1e-6)&&(A_hat(i,j)<=1-(1e-6)))
            theta_hat(i,j)=log(A_hat(i,j)/(1-A_hat(i,j)));
        end
        if (A_hat(i,j)<1e-6)
            theta_hat(i,j)=-inf;
            ind0(i,j)=1;
        end
        if (A_hat(i,j)>1-(1e-6))
            theta_hat(i,j)=+inf;
            ind1(i,j)=1;
        end
        
    end
end

beta_hat=0;
gamma_hat=0;


diffOut=1;

while (diffOut>1e-6)
    beta_hat1=beta_hat;
    gamma_hat1=gamma_hat;
    norm_derv=1;
    while (norm_derv>1e-6)
        derv1=[0; 0];
        for i=1:n
            for j=1:n
                if ((i~=j)&&(ind1(i,j)~=1)&&(ind0(i,j)~=1))
                    derv1=derv1+[D3(i,j)-(D3(i,j)+D4(i,j))*exp(theta_hat(i,j)+beta_hat1)/(1+exp(theta_hat(i,j)+beta_hat1)) ; D5(i,j)-(D5(i,j)+D6(i,j))*exp(theta_hat(i,j)+gamma_hat1)/(1+exp(theta_hat(i,j)+gamma_hat1)) ];
                end
%                 if isnan(derv1)
%                     i
%                     j
%                 end
            end
        end
        derv2=zeros(2);
        for i=1:n
            for j=1:n
                if ((i~=j)&&(ind1(i,j)~=1)&&(ind0(i,j)~=1))
                    derv2(1,1)=derv2(1,1)-(D3(i,j)+D4(i,j))*exp(theta_hat(i,j)+beta_hat1)/(1+exp(theta_hat(i,j)+beta_hat1))^2;
                    derv2(2,2)=derv2(2,2)-(D5(i,j)+D6(i,j))*exp(theta_hat(i,j)+gamma_hat1)/(1+exp(theta_hat(i,j)+gamma_hat1))^2;
                end
            end
        end


        for k=1:2
            if abs(derv2(k,k))<1e-6
                derv1(k)=0;
                derv2(k,k)=1;
            end
        end
       
        adjustNew= [beta_hat1; gamma_hat1]-derv2\derv1;
        beta_hat1New=adjustNew(1);
        gamma_hat1New=adjustNew(2);
        norm_derv=max(abs(derv1));
        beta_hat1=beta_hat1New;
        gamma_hat1=gamma_hat1New;
        
    end
    %%should check diff
    diffOut=max(abs(beta_hat-beta_hat1),abs(gamma_hat-gamma_hat1));
    beta_hat=beta_hat1;
    gamma_hat=gamma_hat1;
    
    
    
    theta_hat1=theta_hat;
    for i=1:(n-1)
        for j=(i+1):n
            if ((i~=j)&&(ind1(i,j)~=1)&&(ind0(i,j)~=1))
                norm_derv=1;
                
                while (norm_derv>1e-6)
                    derv1=D1(i,j)-(D1(i,j)+D2(i,j))*exp(theta_hat1(i,j))/(1+exp(theta_hat1(i,j)))+D3(i,j)-(D3(i,j)+D4(i,j))*exp(theta_hat1(i,j)+beta_hat)/(1+exp(theta_hat1(i,j)+beta_hat))+D5(i,j)-(D5(i,j)+D6(i,j))*exp(theta_hat1(i,j)+gamma_hat)/(1+exp(theta_hat1(i,j)+gamma_hat))+D1(j,i)-(D1(j,i)+D2(j,i))*exp(theta_hat1(j,i))/(1+exp(theta_hat1(j,i)))+D3(j,i)-(D3(j,i)+D4(j,i))*exp(theta_hat1(j,i)+beta_hat)/(1+exp(theta_hat1(j,i)+beta_hat))+D5(j,i)-(D5(j,i)+D6(j,i))*exp(theta_hat1(j,i)+gamma_hat)/(1+exp(theta_hat1(j,i)+gamma_hat));
                    derv2=  -(D1(i,j)+D2(i,j))*exp(theta_hat1(i,j))/(1+exp(theta_hat1(i,j)))^2-(D3(i,j)+D4(i,j))*exp(theta_hat1(i,j)+beta_hat)/(1+exp(theta_hat1(i,j)+beta_hat))^2-(D5(i,j)+D6(i,j))*exp(theta_hat1(i,j)+gamma_hat)/(1+exp(theta_hat1(i,j)+gamma_hat))^2-(D1(j,i)+D2(j,i))*exp(theta_hat1(j,i))/(1+exp(theta_hat1(j,i)))^2-(D3(j,i)+D4(j,i))*exp(theta_hat1(j,i)+beta_hat)/(1+exp(theta_hat1(j,i)+beta_hat))^2-(D5(j,i)+D6(j,i))*exp(theta_hat1(j,i)+gamma_hat)/(1+exp(theta_hat1(j,i)+gamma_hat))^2;
                    if abs(derv1)<1e-6
                        derv1=0;
                        derv2=1;
                    end                  
                    theta_hat1New=theta_hat1(i,j)-derv1/derv2;
    
                    %diff=abs(exp(theta_hat1New)/(1+exp(theta_hat1New))-exp(theta_hat1(i,j))/(1+exp(theta_hat1(i,j))));
                    norm_derv=abs(derv1);
                    theta_hat1(i,j)=theta_hat1New;
                    theta_hat1(j,i)=theta_hat1New;
                end
            end
        end
    end
    

    A_hat=exp(theta_hat)./(1+exp(theta_hat));
    A_hat1=exp(theta_hat1)./(1+exp(theta_hat1));
    diffOut=max(diffOut,max(max(abs(A_hat-A_hat1))));
    theta_hat=theta_hat1;
    %theta_hat

end


for i=1:n
    for j=1:n
        if (ind1(i,j)==1)
            
            A_hat(i,j)=1;
            
            
        end
        if (ind0(i,j)==1)
            A_hat(i,j)=0;
            
        end
    end
end
for i=1:n
    A_hat(i,i)=1;
end

diffOut=1;
alpha_hat=2;
u_hat=zeros(1,n)+1;
freq=sum(S);
u_hat(freq==0)=-inf;
%u_hat(n)=0;
while (diffOut>1e-6)
    diffOut=0;
    for r=1:n
        if u_hat(r)~=-inf
            u_hat1=u_hat;
            alpha_hat1=alpha_hat;
            norm_derv=1;
            while (norm_derv>1e-6)
                derv1=0;
                derv2=0;
                rho=exp(u_hat1);
                derv1=derv1+S(1,r)-sum(S(1,:))*(rho(r)/sum(rho));
                derv2=derv2-sum(S(1,:))*(rho(r)*sum(rho)-rho(r)*rho(r))/(sum(rho))^2;
                
                for j=1:n
                    v=u_hat1;
                    v(j)=v(j)+alpha_hat1;
                    rho=exp(v);
                    derv1=derv1+B(r,j)-sum(B(:,j))*(rho(r)/sum(rho));
                    derv2=derv2-sum(B(:,j))*(rho(r)*sum(rho)-rho(r)*rho(r))/(sum(rho))^2;
                    
                end
                u_hat1(r)=u_hat1(r)-derv1/derv2;
                norm_derv=abs(derv1);
                
            end
            diffOut=max(abs(exp(u_hat1)./sum(exp(u_hat1))-exp(u_hat)./sum(exp(u_hat))));
            u_hat=u_hat1;
        end
    end
    %diffOut
    u_hat1=u_hat;
    alpha_hat1=alpha_hat;
    norm_derv=1;
    while (norm_derv>1e-6)
        derv1=0;
        derv2=0;
        derv1=derv1+sum(diag(B));
        
        for j=1:n
            v=u_hat1;
            v(j)=v(j)+alpha_hat1;
            rho=exp(v);
            derv1=derv1-sum(B(:,j))*(rho(j)/sum(rho));
            
        end
        
        
        for j=1:n
            v=u_hat1;
            v(j)=v(j)+alpha_hat1;
            rho=exp(v);
            derv2=derv2-sum(B(:,j))*(rho(j)*sum(rho)-rho(j)*rho(j))/(sum(rho))^2;
            
        end
        alpha_hat1=alpha_hat1-derv1/derv2;
        norm_derv=abs(derv1);
    end
    diffOut=max(diffOut,abs(alpha_hat1-alpha_hat));
    alpha_hat=alpha_hat1;
end

rho_hat=exp(u_hat)/sum(exp(u_hat));
A_hat=exp(theta_hat)./(1+exp(theta_hat));
for i=1:n
    A_hat(i,i)=1;
end

for i=1:n
    for j=1:n
        if theta_hat(i,j)==inf
            A_hat(i,j)=1;
        end
    end
end


