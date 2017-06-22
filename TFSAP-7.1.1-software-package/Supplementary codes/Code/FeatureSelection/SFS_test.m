% function [Xkp1,xkp1,Jbest] = SFS_test(Xk,dataTrain, LabelTrain, dataTest, LabelTest)
function [Xkp1,xkp1,Jbest] = SFS_test(Xk,D,kernel)


% D=size(dataTrain,2);
index = 1:D;
YDmk = setdiff(1:D,Xk);
Jbest = 0;
Dmk = length(YDmk);
for nF = 1:Dmk
    X = zeros(1,D);
    Xtmp = union(Xk,YDmk(nF));
    X(Xtmp)=1;
%     bestcv = J_cout(dataTrain(:,index(X==1)), LabelTrain, dataTest(:,index(X==1)), LabelTest, gridWidth);
    bestcv = J_cost(index(X==1),kernel);

    J_fctn_k = bestcv;
    
    
    if J_fctn_k > Jbest
        Jbest = J_fctn_k;
        nF_hat = YDmk(nF);
    end
end

xkp1 = nF_hat;
Xkp1 = union(Xk,xkp1);

