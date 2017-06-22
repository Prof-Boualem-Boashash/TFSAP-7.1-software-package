function [featuresOrdered,rank]=filterFisher(train, group)

%This code allows to compute the Fisher criterion for ranking features.
%Note that in this implementation the class should start from 1 and indexed
% till nc, where nc represents te number of classes.
nc = length(unique(group));

for k=1:size(train,2)
    for c=1:nc
        posC=find(group==(c));
        tempNum(c)=nc*((mean(train(posC,k))-mean(train(:,k)))^2);
        tempDen(c)=nc*var(train(posC,k));
    end
    rank(k)=sum(tempNum)/sum(tempDen);
end
[a,featuresOrdered]=sort(rank,'descend');