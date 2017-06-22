% function [Xk,xk,Jbest] = SBS_test(Xkp1,dataTrain, LabelTrain, dataTest, LabelTest, gridWidth)
function [Xk,xk,Jbest] = SBS_test(Xkp1,D, kernel)


% D=size(dataTrain,2);
index = 1:D;
Jbest = 0;
Dk = length(Xkp1);

for nB = 1:Dk
    X = zeros(1,D);
    Xtmp = Xkp1;
    Xtmp(nB) = [];
    X(Xtmp) = 1;
%     bestcv = J_cout(dataTrain(:,index(X==1)), LabelTrain, dataTest(:,index(X==1)), LabelTest, gridWidth);
    bestcv = J_cost(index(X==1),kernel);

    J_fctn_k = bestcv;
    
    
    if J_fctn_k > Jbest
        Jbest = J_fctn_k;
        nB_hat = Xkp1(nB);
    end
end


if ~isempty(Xkp1)
    xk = nB_hat;
    Xk = setdiff(Xkp1,xk);
else
    xk = [];
    Xk = [];
end

