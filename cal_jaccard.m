function jaccard=cal_jaccard(G1,G2)

    jaccard=sum((G1==1)&(G2==1))/sum((G1==1)|(G2==1));
