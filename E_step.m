function [S,B]=E_step(G,theta,u,alpha,beta,gamma)

T=size(G,1);
n=size(G,2);

A=exp(theta)./(1+exp(theta));
for i=1:n
    A(i,i)=1;
end

for i=1:n
    for j=1:n
        if theta(i,j)==inf
            A(i,j)=1;
        end
    end
end


aa=zeros(T,n);
t=1;
rho=exp(u)/sum(exp(u));
Phi=zeros(n,n);
for j=1:n  %%%i: row index, time t j: column index, time t-1
    v=u;
    v(j)=v(j)+alpha;
    Phi(:,j)=exp(v)/sum(exp(v));
end

for i=1:n
    if G(t,i)==1
        pr_z=rho(i);
        pr_G_cond=1;
        for j=1:n
            if (j~=i)
                
                pr_G_cond=pr_G_cond*(G(t,j)*A(i,j)+(1-G(t,j))*(1-A(i,j)));
            end
        end
        aa(t,i)=pr_z*pr_G_cond;
    end
end

for t=2:T
    for i=1:n
        if G(t,i)==1
            pr_G_cond=1;
            if G(t-1,i)==0
                for j=1:n
                    if (j~=i)
                        
                        pr_G_cond=pr_G_cond*(G(t,j)*A(i,j)+(1-G(t,j))*(1-A(i,j)));
                    end
                end
            else
                for j=1:n
                    if (j~=i)
                        if G(t-1,j)==1
                            pr=exp(theta(i,j)+beta)./(1+exp(theta(i,j)+beta));
                            if isnan(pr)>0
                                pr=1;
                            end
                            pr_G_cond=pr_G_cond*(G(t,j)*pr+(1-G(t,j))*(1-pr));
                        else
                            pr=exp(theta(i,j)+gamma)./(1+exp(theta(i,j)+gamma));
                            if isnan(pr)>0
                                pr=1;
                            end
                            pr_G_cond=pr_G_cond*(G(t,j)*pr+(1-G(t,j))*(1-pr));
                        end
                        
                    end
                end
            end
            for i_prev=1:n
                pr_z=Phi(i,i_prev);
                aa(t,i)=aa(t,i)+aa(t-1,i_prev)*pr_z*pr_G_cond;

%                 aa(t,i)
%                 pr_z
%                 pr_G_cond
            end
            
        end
    end
    aa(t,:)=aa(t,:)/sum(aa(t,:));

end

bb=zeros(T,n);

bb(T,:)=1;

for t=(T):(-1):2
    for i_prev=1:n
        
        for i=1:n
            if G(t,i)==1
                
                pr_z=Phi(i,i_prev);
                
                pr_G_cond=1;
                if G(t-1,i)==0
                    for j=1:n
                        if (j~=i)
                            
                            pr_G_cond=pr_G_cond*(G(t,j)*A(i,j)+(1-G(t,j))*(1-A(i,j)));
                        end
                    end
                else
                    for j=1:n
                        if (j~=i)
                            if G(t-1,j)==1
                                pr=exp(theta(i,j)+beta)./(1+exp(theta(i,j)+beta));
                                if isnan(pr)>0
                                    pr=1;
                                end
                                pr_G_cond=pr_G_cond*(G(t,j)*pr+(1-G(t,j))*(1-pr));

                            else
                                pr=exp(theta(i,j)+gamma)./(1+exp(theta(i,j)+gamma));
                                if isnan(pr)>0
                                    pr=1;
                                end
                                pr_G_cond=pr_G_cond*(G(t,j)*pr+(1-G(t,j))*(1-pr));

                            end
                            
                        end
                    end
                end
                bb(t-1,i_prev)=bb(t-1,i_prev)+bb(t,i)*pr_z*pr_G_cond;
            end
        end
        bb(t,:)=bb(t,:)/sum(bb(t,:));
    end
end

pr_single=aa.*bb;
sum_single= repmat(sum(pr_single,2),1,n);
pr_single=pr_single./sum_single;

cc=zeros(T,n);
cc(1,:)=nan;

t=T;

for i=1:n
    if G(t,i)==1
        pr_G_cond=1;
        if G(t-1,i)==0
            for j=1:n
                if (j~=i)
                    
                    pr_G_cond=pr_G_cond*(G(t,j)*A(i,j)+(1-G(t,j))*(1-A(i,j)));
                end
            end
        else
            for j=1:n
                if (j~=i)
                    if G(t-1,j)==1
                        pr=exp(theta(i,j)+beta)./(1+exp(theta(i,j)+beta));
                        if isnan(pr)>0
                            pr=1;
                        end
                        pr_G_cond=pr_G_cond*(G(t,j)*pr+(1-G(t,j))*(1-pr));
                    else
                        pr=exp(theta(i,j)+gamma)./(1+exp(theta(i,j)+gamma));
                        if isnan(pr)>0
                            pr=1;
                        end
                        pr_G_cond=pr_G_cond*(G(t,j)*pr+(1-G(t,j))*(1-pr));
                    end
                    
                end
            end
        end
        cc(t,i)=pr_G_cond;
        
    end
end

for t=(T-1):(-1):2
    for i=1:n
        if G(t,i)==1
            pr_G_cond=1;
            if G(t-1,i)==0
                for j=1:n
                    if (j~=i)
                        
                        pr_G_cond=pr_G_cond*(G(t,j)*A(i,j)+(1-G(t,j))*(1-A(i,j)));
                    end
                end
            else
                for j=1:n
                    if (j~=i)
                        if G(t-1,j)==1
                            pr=exp(theta(i,j)+beta)./(1+exp(theta(i,j)+beta));
                            if isnan(pr)>0
                                pr=1;
                            end
                            pr_G_cond=pr_G_cond*(G(t,j)*pr+(1-G(t,j))*(1-pr));
                        else
                            pr=exp(theta(i,j)+gamma)./(1+exp(theta(i,j)+gamma));
                            if isnan(pr)>0
                                pr=1;
                            end
                            pr_G_cond=pr_G_cond*(G(t,j)*pr+(1-G(t,j))*(1-pr));
                        end
                        
                    end
                end
            end
            
            for i_after=1:n
                pr_z=Phi(i_after,i);
                
                cc(t,i)=cc(t,i)+cc(t+1,i_after)*pr_z*pr_G_cond;
            end
            
        end
        
    end
    cc(t,:)=cc(t,:)/sum(cc(t,:));
end

pr_double=zeros(n,n,T);
pr_double(:,:,1)=nan;

for t=2:T
    for i=1:n
        for j=1:n
            
            pr_double(i,j,t)=aa(t-1,j)*Phi(i,j)*cc(t,i);
        end
    end
end
sum_double=sum(sum(pr_double,1),2);
for t=1:T
    pr_double(:,:,t)=pr_double(:,:,t)/sum_double(t);
end

S=pr_single;
B=sum(pr_double(:,:,2:T),3);


