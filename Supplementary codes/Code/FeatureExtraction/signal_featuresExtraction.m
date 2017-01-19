function [ TFC ] = signal_featuresExtraction(tfd)
%*************************************************
%
% Copyright 2016 Boualem Boashash
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Authors: Boualem Boashash         (boualem.boashash@gmail.com)
%          Samir Ouelha  			(samir_ouelha@hotmail.fr)
%
% The following references should be cited whenever this script is used:
% [1] B. Boashash, Samir Ouelha, Designing high-resolution time-frequency
% and time-scale distributions for the analysis and  classification of non 
% stationary signals: a tutorial review with features performance comparison
%, Digital Signal Processing, In Press.
% [2] B. Boashash, Samir Ouelha, Efficient Software Matlab Package to compute
% Time-Frequency Distributions and related Time-Scale methods with extraction
% of signal characteristics, SoftwareX, In Press.
%
% This study was funded by grants from the ARC and QNRF NPRP 6-885-2-364
% and NPRP 4-1303-2-517 
%
%
% Description:
% This function extracts different features from a TFD.


%% time-frequency entropy flatness
[N,M]=size(tfd);
TF_f1 = (prod(abs(tfd(:)).^(1/(N*M))))/mean(abs(tfd(:))); 
tfd_temp=tfd;

%% time-frequency Renyi entropy
tfdn=tfd/sum(tfd(:));
TF_f2=0.5*log(sum((tfdn(:).^3)));  % Renyi entropy

%% Time-frequency flux
LL=1;
QQ=1;

TF_11=abs(tfd(1:end-LL,1:end-QQ)-tfd(LL+1:end,QQ+1:end)); 
TF_f3 = sum(TF_11(:))/sum(tfd(:));

LL=1;

TF_10=abs(tfd(1:end-LL,:)-tfd(LL+1:end,:));
TF_f4 = sum(TF_10(:))/sum(tfd(:));

QQ=1;

TF_01=abs(tfd(:,1:end-QQ)-tfd(:,QQ+1:end));
TF_f5=sum(TF_01(:))/sum(tfd(:)); 


%% Energy concentration
P=2;
TF_f6=10*log10((sum(abs(tfdn(:)).^(1/P)))^(P));


%% Shannon entropy
% tfd(tfd<=0)=eps;
% tfd=tfd/sum(tfd(:));
tfdG = TransformToimgray(tfd,256);
TF_f7 = entropy(tfdG);


%% Statistical feature
tfd=tfd_temp;
TF_t1=mean(tfd(:));
TF_t2=std(tfd(:));
TF_t3=TF_t2/TF_t1;
TF_t4=skewness(tfd(:));
TF_t5=kurtosis(tfd(:));


TFC=[TF_f1 TF_f2 TF_f3 TF_f4 TF_f5 TF_f6 TF_f7 TF_t1 TF_t2 TF_t3 TF_t4 TF_t5  ];


%% Image features
tfdBin = binary_segmented_TF_image(tfd,20);
TF_IM=TF_image_features(tfdBin,1);


%% Concatenation
TFC=[TFC  TF_IM]; % TF translated, TF Image, TF IF related

end