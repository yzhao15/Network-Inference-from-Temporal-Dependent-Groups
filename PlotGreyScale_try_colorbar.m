% 25APR15
% Author: Charles Weko

% 
% This is a function to generate the grey scale visualization that we have
% been using.  I am putting it into MATLAB to reduce the amount of time
% that I have to spend screwing around with it.
%

function PlotGreyScale_try_colorbar(A,clusterLabel)

colormap(flipud(gray(256)));
colorbar;


nodes=size(A,2);
T=size(A,1);

if nargin==1
    nodeNames=num2str(1:nodes);
end %if

for i=1:T
    for j=1:nodes     
        a=1-A(i,j);  % 1 is white and 0 is black
        rectangle('position', [i, j, 1, 1],'facecolor',[a a a],'edgecolor','white')
    end %j
end %i

K=max(clusterLabel);
for k=1:K
    ind=find(clusterLabel==k);
    x=linspace(min(ind),max(ind)+1);
    y=linspace(max(ind)+1,max(ind)+1);
    line(x,y,'Color','red','LineWidth',2);
    y=linspace(min(ind),max(ind)+1);
    x=linspace(max(ind)+1,max(ind)+1);
    line(x,y,'Color','red','LineWidth',2);
    x=linspace(min(ind),max(ind)+1);
    y=linspace(min(ind),min(ind));
    line(x,y,'Color','red','LineWidth',2);
    y=linspace(min(ind),max(ind)+1);
    x=linspace(min(ind),min(ind));
    line(x,y,'Color','red','LineWidth',2);

end

axis([1 T+1 1 nodes+1])
 axis ij
 
 set(gca,'XTick',(1:nodes)+0.5);
 set(gca,'XTickLabel','');
set(gca, 'XAxisLocation', 'top')
 set(gca,'YTick',(1:nodes)+0.5);
 set(gca,'YTickLabel','');
 xlabel('Chimpanzees');
 ylabel('Chimpanzees');
