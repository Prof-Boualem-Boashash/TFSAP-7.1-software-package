% function [Xbest, J_hat_best]=SFFS_test(dataTrain, LabelTrain, dataTest, LabelTest,gridWidth)
function [Xbest, J_hat_best]=SFFS_test(D,kernel)

% D=size(dataTrain,2);
% J_hat = zeros(1,D);
J_hat_best = 0;
Xk = [];
Xbest = [];
k = 0;


% SFS(k=2)
while k<2
    [Xkp1,xkp1,Jbest] = SFS_test(Xk,D,kernel);
    Xk = Xkp1;
    J_hat(k+1) = Jbest;
    k = k+1;
end


% SFFS (2<=k<=d<=D)
while k<D
    %Inclusion
[Xkp1,xkp1,JbestF] = SFS_test(Xk,D,kernel);

    %Test
[Xr,xr,JbestBr] = SBS_test(Xkp1,D,kernel);
    if xr==xkp1
        Xk = Xkp1;
        k = k+1;
        J_hat(k) = JbestF;
    else

        if JbestBr <= J_hat(k)
            Xk = Xkp1;
            k = k+1;
            J_hat(k) = JbestF;
        else

            if k==2
                Xk = Xr;
                J_hat(k) = JbestBr;
            else

                %Exclusion
                Xs = Xr;
                while k>2
                    [Xsm1,xsm1,JbestBs] = SBS_test(Xs,D,kernel);
                    if JbestBs <= J_hat(k-1)
                        break;
                    else
                        Xs = Xsm1;
                        k = k-1;
                    end
                end
                %si k=2 ou JbestBs <= J_hat(k-1)
                Xk = Xs;
                J_hat(k) = JbestBs;
            end
        end
    end

    if J_hat(k) > J_hat_best
        J_hat_best = J_hat(k)
        Xbest = Xk
    end

end
