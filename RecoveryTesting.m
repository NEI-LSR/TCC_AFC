function [nll,aic] = RecoveryTesting(realOrSimData,rn)

% Comparing model outputs

if ~exist('realOrSimData','var')
    realOrSimData = 'sim';
    warning('No value for `realOrSimData` provided. Defaulting to simulated data')
end

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

if ~exist('rn','var')
    rn = 0;
    warning('Using default rn (0)')
end

%%

[~,~,~,nll.dPrimeAndGaussianWidth,~] = ParameterEstimator_caller(rn,data,'dPrimeAndGaussianWidth');
aic.dPrimeAndGaussianWidth = aicbic(-nll.dPrimeAndGaussianWidth,2);

%%

[~,~,~,nll.single_ssnu,~] = ParameterEstimator_caller(rn,data,'single-ssnu');
aic.single_ssnu = aicbic(-nll.single_ssnu,64);

%%

[~,~,~,nll.single_sg,~] = ParameterEstimator_caller(rn,data,'single-sg');
aic.single_sg = aicbic(-nll.single_sg,64);

%%

[~,~,~,nll.simultaneous,~] = ParameterEstimator_caller(rn,data,'simultaneous');
aic.simultaneous = aicbic(-nll.simultaneous,130);

%%

% [~,~,~,nll.iterative,~] = ParameterEstimator_caller(rn,data,'iterative');
% aic.iterative = aicbic(-nll.iterative,130);

%%

end