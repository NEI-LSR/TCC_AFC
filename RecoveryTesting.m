% Comparing model outputs

clear, clc, close all

realOrSimData = 'sim';

if strcmp(realOrSimData,'sim')

    data = generateSimulatedData;

elseif strcmp(realOrSimData,'real')

    addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))
    DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';
    if ~exist('data','var')
        data = combineData_mat(DataDir); % TODO Switch to csv
    end

    data.trialdata.nBig = 64;
    data.trialdata.nSmall = 4;
    data.trialdata.nTrials = 98104;
    warning('Assuming nBig, nSmall, and nTrials')

end

%% 

rn = 1;

%%

ParameterEstimator_caller(rn,data,'single-ssnu');

%%

ParameterEstimator_caller(rn,data,'single-sg');

%%

ParameterEstimator_caller(rn,data,'simultaneous');

%%

ParameterEstimator_caller(rn,data,'iterative');

