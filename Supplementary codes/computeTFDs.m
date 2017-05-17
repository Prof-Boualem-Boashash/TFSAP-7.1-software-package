%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Authors: Prof. Boualem Boashash   (boualem.boashash@gmail.com)
%          Samir Ouelha             (samir_ouelha@hotmail.fr)
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
% Last Modification: 25-12-2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
%
% Computes TFDs and makes signals ready for further processing
% Kernel can take the following values:
% kernel = {'wvd' 'spec' 'embd' 'ckd' 'mdd' 'awvd' 'scal' 'spawvd'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tfd = computeTFDs(kernel, varargin)
tfd = [];
TFDs = {'WVD','SPEC','EMBD','CKD','MDD','AWVD','SCALO','SPAWVD'};

%% Main Inputs Checkup
if(nargin == 0), error_msg(1); return;
elseif(nargin >= 1)
    if(isempty(kernel)), error_msg(2); return;
    elseif(~sum(strcmpi(kernel,TFDs)))
        error_msg(3); disp(TFDs'); return;
    end
end

%% Auxiliary Inputs Checkup
switch lower(kernel)
    case 'wvd'
        if(nargin == 1), tres = 1;
        elseif(nargin == 2)
            if(isempty(varargin{1})), tres = 1;
            elseif(length(varargin{1})>1), error_msg(4); return;
            elseif(~isnumeric(varargin{1})), error_msg(4); return;
            elseif(mod(varargin{1},1)), error_msg(4); return;
            elseif(varargin{1}<=0), error_msg(4); return;
            else tres = varargin{1};
            end
        else error_msg(5); return;
        end
    case 'spec'
        if(nargin == 1), tres = 1; window_length = 65;
        elseif(nargin == 2)
            if(isempty(varargin{1})), tres = 1;
            elseif(length(varargin{1})>1), error_msg(4); return;
            elseif(~isnumeric(varargin{1})), error_msg(4); return;
            elseif(mod(varargin{1},1)), error_msg(4); return;
            elseif(varargin{1}<=0), error_msg(4); return;
            else tres = varargin{1};
            end
            window_length = 65;
        elseif(nargin == 3)
            if(isempty(varargin{1})), tres = 1;
            elseif(length(varargin{1})>1), error_msg(4); return;
            elseif(~isnumeric(varargin{1})), error_msg(4); return;
            elseif(mod(varargin{1},1)), error_msg(4); return;
            elseif(varargin{1}<=0), error_msg(4); return;
            else tres = varargin{1};
            end
            if(isempty(varargin{2})), window_length = 65;
            elseif(length(varargin{2})>1), error_msg(6); return;
            elseif(~isnumeric(varargin{2})), error_msg(6); return;
            elseif(mod(varargin{2},1)), error_msg(6); return;
            elseif(~mod(varargin{2},2)), warning_msg(1,varargin{2}-1); window_length = varargin{2}-1;
            else window_length = varargin{2};
            end
        else error_msg(5); return;
        end
    case 'embd'
        if(nargin == 1), tres = 1; alpha = 0.05; beta = 0.05;
        elseif(nargin == 2)
            if(isempty(varargin{1})), tres = 1;
            elseif(length(varargin{1})>1), error_msg(4); return;
            elseif(~isnumeric(varargin{1})), error_msg(4); return;
            elseif(mod(varargin{1},1)), error_msg(4); return;
            elseif(varargin{1}<=0), error_msg(4); return;
            else tres = varargin{1};
            end
            alpha = 0.05; beta = 0.05;
        elseif(nargin == 3)
            if(isempty(varargin{1})), tres = 1;
            elseif(length(varargin{1})>1), error_msg(4); return;
            elseif(~isnumeric(varargin{1})), error_msg(4); return;
            elseif(mod(varargin{1},1)), error_msg(4); return;
            elseif(varargin{1}<=0), error_msg(4); return;
            else tres = varargin{1};
            end
            if(isempty(varargin{2})), alpha = 0.05;
            elseif(length(varargin{2})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{2})), error_msg(7); return;
            elseif(varargin{2}<0), error_msg(7); return;
            else alpha = varargin{2};
            end
            beta = 0.05;
        elseif(nargin == 4)
            if(isempty(varargin{1})), tres = 1;
            elseif(length(varargin{1})>1), error_msg(4); return;
            elseif(~isnumeric(varargin{1})), error_msg(4); return;
            elseif(mod(varargin{1},1)), error_msg(4); return;
            elseif(varargin{1}<=0), error_msg(4); return;
            else tres = varargin{1};
            end
            if(isempty(varargin{2})), alpha = 0.05;
            elseif(length(varargin{2})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{2})), error_msg(7); return;
            elseif(varargin{2}<0), error_msg(7); return;
            else alpha = varargin{2};
            end
            if(isempty(varargin{3})), beta = 0.05;
            elseif(length(varargin{3})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{3})), error_msg(7); return;
            elseif(varargin{3}<0), error_msg(7); return;
            else beta = varargin{3};
            end
        else error_msg(5); return;
        end
    case 'ckd'
        if(nargin == 1), C = 1; D = 0.04; E = 0.04;
        elseif(nargin == 2)
            if(isempty(varargin{1})), C = 1;
            elseif(length(varargin{1})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{1})), error_msg(7); return;
            elseif(varargin{1}<0), error_msg(7); return;
            else C = varargin{1};
            end
            D = 0.04; E = 0.04;
        elseif(nargin == 3)
            if(isempty(varargin{1})), C = 1;
            elseif(length(varargin{1})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{1})), error_msg(7); return;
            elseif(varargin{1}<0), error_msg(7); return;
            else C = varargin{1};
            end
            if(isempty(varargin{2})), D = 0.04;
            elseif(length(varargin{2})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{2})), error_msg(7); return;
            elseif(varargin{2}<0), error_msg(7); return;
            else D = varargin{2};
            end
            E = 0.04;
        elseif(nargin == 4)
            if(isempty(varargin{1})), C = 1;
            elseif(length(varargin{1})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{1})), error_msg(7); return;
            elseif(varargin{1}<0), error_msg(7); return;
            else C = varargin{1};
            end
            if(isempty(varargin{2})), D = 0.04;
            elseif(length(varargin{2})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{2})), error_msg(7); return;
            elseif(varargin{2}<0), error_msg(7); return;
            else D = varargin{2};
            end
            if(isempty(varargin{3})), E = 0.04;
            elseif(length(varargin{3})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{3})), error_msg(7); return;
            elseif(varargin{3}<0), error_msg(7); return;
            else E = varargin{3};
            end
        else error_msg(5); return;
        end
    case 'mdd'
        if(nargin == 1), C = 1; thr = 0.2;
        elseif(nargin == 2)
            if(isempty(varargin{1})), C = 1;
            elseif(length(varargin{1})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{1})), error_msg(7); return;
            elseif(varargin{1}<0), error_msg(7); return;
            else C = varargin{1};
            end
            thr = 0.2;
        elseif(nargin == 3)
            if(isempty(varargin{1})), C = 1;
            elseif(length(varargin{1})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{1})), error_msg(7); return;
            elseif(varargin{1}<0), error_msg(7); return;
            else C = varargin{1};
            end
            if(isempty(varargin{2})), thr = 0.2;
            elseif(length(varargin{2})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{2})), error_msg(7); return;
            elseif(varargin{2}<0), error_msg(7); return;
            else thr = varargin{2};
            end
        else error_msg(5); return;
        end
    case 'awvd'
        if(nargin > 1), error_msg(5); return; end
    case 'scalo'
        if(nargin > 1), error_msg(5); return; end
    case 'spawvd'
        if(nargin == 1), nh = 16; ng = 6;
        elseif(nargin == 2)
            if(isempty(varargin{1})), nh = 16;
            elseif(length(varargin{1})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{1})), error_msg(7); return;
            elseif(varargin{1}<0), error_msg(7); return;
            else nh = varargin{1};
            end
            ng = 6;
        elseif(nargin == 3)
            if(isempty(varargin{1})), nh = 16;
            elseif(length(varargin{1})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{1})), error_msg(7); return;
            elseif(varargin{1}<0), error_msg(7); return;
            else nh = varargin{1};
            end
            if(isempty(varargin{2})), ng = 6;
            elseif(length(varargin{2})>1), error_msg(7); return;
            elseif(~isnumeric(varargin{2})), error_msg(7); return;
            elseif(varargin{2}<0), error_msg(7); return;
            else ng = varargin{2};
            end
        else error_msg(5); return;
        end
end

%% Adding all Required Paths
addpath(genpath([pwd '\EEGDataSet\Filtered']))
addpath(genpath([pwd '\code\TS analysis']))
%replace "path_to_TFSAP_ToolBox" with the actual path to the toolbox on
%your PC
%addpath('path_to_TFSAP_ToolBox');

%% Constant Parameters
fs = 32;   % sampling frequency for preprocessed data
T  = 8;    % epoch time in seconds
N  = T*fs; % epoch samples
Np = 36;    % number of patients

%% Main
for i_patient = 1:Np
    % Load filtered data
    load(['preprocessed_data_' num2str(i_patient)],'sig_channel_averaged','eegStateN');
    % number of epochs of sample N each in the averaged channel signal
    number_epochs = floor(length(sig_channel_averaged)/N);
    if floor(length(sig_channel_averaged)/N)==length(sig_channel_averaged)/N
        number_epochs = number_epochs-1;
    end
    disp([num2str(i_patient) ': the number of epochs for this patient is ' num2str(number_epochs)])
    for i_epoch = 1:number_epochs
        sig_start = N*(i_epoch-1) + 1;
        sig_end   = N*(i_epoch-1) + N;
        sig_cur_epoch = sig_channel_averaged(sig_start:sig_end);
        switch lower(kernel)
            case 'wvd'
                tfd = quadtfd(sig_cur_epoch, N-1, tres, 'wvd',N);
            case 'spec'
                tfd = quadtfd(sig_cur_epoch, N-1, tres, 'specx',window_length,'hamm',N);
            case 'embd'
                tfd = quadtfd(sig_cur_epoch, N-1, tres, 'emb',alpha, beta, N);
            case 'ckd'
                tfd = cmpt(sig_cur_epoch,'ckd', C,D,E);
            case 'mdd'
                tfd = myMDD(sig_cur_epoch(:)',C,thr,0,1);
            case 'awvd'   % Affine Wigner-Ville
                tfd = awvd(sig_cur_epoch, N);
            case 'scalo'  % Scalogram
                tfd = scalogram(sig_cur_epoch, 1, N);
            case 'spawvd' % Pseudo Smoothed AWVD
                tfd = spawd(sig_cur_epoch, nh, ng, N);
        end
        eval([kernel '_epoch_' num2str(i_epoch) '=tfd;']);
    end
    pathPatient = [pwd '\EEGDataSet\TFDs\' kernel '\' kernel 'preprocessed_TFD_' num2str(i_patient)];
    save(pathPatient);
    clear([kernel '_epoch_*']);
end
end
%% Supplementary Functions
function error_msg(n)
switch n,
    case 1, fprintf(2,'ERROR (No input kernel)\n');
    case 2, fprintf(2,'ERROR (Input kernel is empty)\n');
    case 3, fprintf(2,'ERROR (Unknown TFD kernel)\n'); fprintf(2,'Kindly choose one of the following TFDs\n');
    case 4, fprintf(2,'ERROR (Time resolution must be a 1x1 positive integer)\n');
    case 5, fprintf(2,'ERROR (Extra inputs, please check the function help)\n');
    case 6, fprintf(2,'ERROR (Window length must be a 1x1 odd integer)\n');
    case 7, fprintf(2,'ERROR (Smoothing parameter must be a 1x1 positive number)\n');
end
end
function warning_msg(n, x)
switch n,
    case 1, fprintf('WARNING (Window length must be odd, thus it is truncated to %d)\n',x);
end
end
