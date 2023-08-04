function ParameterEstimator_caller(rn)

%clear, clc, close all

% Currently some manual handling is required
% If running the model to fit dprime and gaussianwidth:
% comment out `'dPrime',       data.trialdata.dPrime,...    'gaussianWidth',data.trialdata.gaussianWidth,...` 
% in ParameterEstimator.m 

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))

% data = generateSimulatedData;

DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';

%rn = 0;

%%
data = combineData_mat(DataDir); % TODO Switch to csv

data.trialdata.nBig = 64;
data.trialdata.nTrials = 98104;
data.trialdata.nSmall = 4;

%% Add in stimulus remapping

load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\ssnu\4_230803-161235.mat',...
    'x')
data.trialdata.stimulusRemappingPol = x;

%%

params.dPrime               = true;
params.gaussianWidth        = true;

% [x,aic,bic] = ParameterEstimator(data,params,rn);

% disp(x)

%%

% TODO See if we can define f at this level

% data.trialdata.dPrime           = x(1);
% data.trialdata.gaussianWidth    = x(2);

data.trialdata.dPrime           = 1.4899;
data.trialdata.gaussianWidth    = 39.0070;

params = rmfield(params,'dPrime'); % TODO It would be nice if we could just set these to false rather than have to delete them
params = rmfield(params,'gaussianWidth');

% params.stimulusRemapping = true;
params.skewedGaussians = true;

[x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);

% disp(x)
disp(aic)
disp(nll_x)

% TODO [#C] automate the switch in naming so that they go into the right folders
save(['C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\',...
    'sg_ssnu\',num2str(rn),'_',datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


end