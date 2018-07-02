function PlotGreyScale_nonSquare_try(A,z,R)

nodes=size(A,2);
T=size(A,1);

if nargin==1
    nodeNames=num2str(1:nodes);
end %if

for j=1:T
    if R(j)==1
        
        line([1+0.2 nodes+1-0.2], [j j], 'color','blue','LineWidth',2);
    end
    for i=1:nodes
        %a=(1-A(j,i));  % 1 is white and 0 is black
        if A(j,i)==0
            a=1;
        else
            a=0.5;
        end
        rectangle('position', [i, j, 1, 1],'facecolor',[a a a],'edgecolor','white')
        if i==z(j)
            rectangle('position', [i, j, 1, 1],'facecolor','red','edgecolor','white')
        end
    end %j
end %i
axis([1 nodes+1 1 T+1 ])
axis ij
%
set(gca,'XTick',(1:nodes)+0.5);
set(gca,'XTickLabel','');
set(gca, 'XAxisLocation', 'top')
set(gca,'YTick',(1:T)+0.5);
set(gca,'YTickLabel','');

xlabel('Chimpanzees')
ylabel('Groups over time')