%% Parameters
Np = 36;
Nf = 89;
N_maxEpoch = 986;

signal_features_norm = zeros(Np,N_maxEpoch,Nf);

num_patient = Data(:,end);
data_norm_new = zeros(1,size(data_norm,1),size(data_norm,2)-1);

for k=1:size(data_norm,1)
    for l=1:size(data_norm,2)-1
        data_norm_new(1,k,l) = data_norm(k,l);
    end
end

for i_patient = 1:Np
    signal_features_norm(i_patient,1:length_features(i_patient),:) = data_norm_new(1,find(num_patient==i_patient),:);
end