function ACC = LOO(length_features_spec, signal_features_spec, signal_class_spec, index)


Y=[];
Output=[];


Np = length(length_features_spec);



for tt=1:Np
    mask=[];
    TEST=[];
    normal_feature_vector=[];
    artifact_feature_vector=[];
    seizure_feature_vector=[];
    k1=0;
    k2=0;
    k3=0;
    
    %Data for learning
    for ii=1:Np
        if ii~=tt % one patient out
            
            for jj=1:length_features_spec(ii)
                if  signal_class_spec(ii,jj)==3
                    k1=k1+1;
                    z(1:length(index))=signal_features_spec(ii,jj,index);
                    normal_feature_vector(k1,:)=z;
                elseif signal_class_spec(ii,jj)==4
                    k2=k2+1;
                    z(1:length(index))=signal_features_spec(ii,jj,index);
                    artifact_feature_vector(k2,:)=z;
                else
                    k3=k3+1;
                    z(1:length(index))=signal_features_spec(ii,jj,index);
                    seizure_feature_vector(k3,:)=z;
                end  
            end
           
        end
    end
    
    Training=[seizure_feature_vector( :,:);normal_feature_vector(:,:);artifact_feature_vector(:,:)];
    Group     = [ones(k3,1); zeros(k1,1);-ones(k2,1)];    
        
    for jj=1:(length_features_spec(tt))
        z(1:length(index))=signal_features_spec(tt,jj,index);
        TEST(jj,:)=z;
        if or(signal_class_spec(tt,jj)==3,signal_class_spec(tt,jj)==4) %non seizure and artefact
            mask(jj)=0;
        else 
            mask(jj)=1;
        end
    end
    
    Y=[Y mask];
    
    B = TreeBagger(50,Training,Group, 'Method', 'classification');
%         Output1 =  multisvm(Training,Group,TEST,1, 1, 1);
%         Output1 = Output1-2;
        Output1 =  B.predict(TEST);
    Output1 = str2double(Output1);

%     SVMStruct1 =  svmtrain(Training,Group,'boxconstraint',1,'kernel_function','rbf','rbf_sigma',1,'method','SMO','SMO_opts',options);
%     Output1 = svmclassify(SVMStruct1,TEST);%1 vs 2

    indNan = isnan(Output1);
    Output1(indNan) = 0;
    

    
    toto(tt) = sum(mask(:)==Output1(:))/length(mask);
    
    
    
    Output=[Output Output1'];   
end

%% Results
ACC = sum(Y==Output);
% index
