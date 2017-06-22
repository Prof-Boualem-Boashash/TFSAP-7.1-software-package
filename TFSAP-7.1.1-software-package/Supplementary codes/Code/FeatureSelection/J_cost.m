% function Perf = J_cout(dataTrain, LabelTrain, dataTest, LabelTest, gridWidth)
function Perf = J_cost(indexR,kernel)
cur_dir=pwd;
load([cur_dir '\result_' kernel])
rng(1);
[~, ~, ~, Perf] = performanceEstimationRandomForests(length_features, signal_features_reduced, signal_class, indexR);
% Perf = dataMix_Perf_RandomForest(X_Artefact,X_Normal,X_suppression,X_burst,X_seizure,20,indexR);
% Perf = LOO(length_features, signal_features_reduced, signal_class, indexR);
% Perf = classificationEstimation(data(:,indexR),groups,2);
% Perf = sum(predictedOutput(:)==Group(:))*100/length(Group)