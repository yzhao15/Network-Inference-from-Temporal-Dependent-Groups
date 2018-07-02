clear


fileName='chimp.xlsx';

G1=xlsread(fileName);
n=size(G1,2)-1;
T=size(G1,1);

%%The last six chimps don't belong to any cluster
clusterLabel=[zeros(1,13)+1 zeros(1,3)+2 zeros(1,2)+3 zeros(1,5)+4 zeros(1,3)+5 zeros(1,2)+6 zeros(1,3)+7 zeros(1,2)+8 zeros(1,4)+9 zeros(1,2)+10 zeros(1,6)];



day=G1(:,(n+1));
G1=G1(:,1:n);


daySeq=find(day);

K=size(daySeq,1);

G=zeros(K,n);


daySeq=[daySeq; T+1]; 

for k=1:K
    t1=daySeq(k);
    t2=daySeq(k+1)-1;
    if (t1==t2)
        G(k,:)=G1(t1,:);
    end
    
    if (t2>t1)
        jaccard=zeros(1,t2-t1+1);
        for t=t1:t2
            jaccard(t-t1+1)=cal_jaccard(G(k-1,:),G1(t,:));
        end
        
        [a,b]=max(jaccard);
        G(k,:)=G1(b+t1-1,:);
    end
    
end

ind=find(max(G));

G=G(:,ind);
clusterLabel=clusterLabel(ind);

%%Data cleaning until this line

[A_indp, rho_indp,S_indp]=HM_EM(G);
[theta_EM,A_EM,u_EM,rho_EM,alpha_EM,beta_EM,gamma_EM,S_EM]=EM(G);

n=size(G,2);


figure
PlotGreyScale_try(A_indp,clusterLabel)
figure
PlotGreyScale_try_colorbar(A_EM,clusterLabel)

fileName=['estimate_para'];

theta=theta_EM;
u=u_EM;
alpha=alpha_EM;
beta=beta_EM;
gamma=gamma_EM;
save(fileName,'theta','u','alpha','beta','gamma');


T=size(S_EM,1);
n=size(S_EM,2);

z_SM=zeros(1,T);

for t=1:T;
    [~,b]=max(S_EM(t,:));
    z_SM(t)=b;
    
    
end

R=zeros(1,T);

for t=2:T
    if G(t-1,z_SM(t))==0
        R(t)=1;
    end
end
figure
PlotGreyScale_nonSquare_try(G,z_SM,R);

