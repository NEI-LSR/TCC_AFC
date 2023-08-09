function ParameterEstimator_caller(rn,data)

clear, clc, close all

if ~exist('rn','var')
    rn = 0;
    warning('Using default value for `rn` (0)')
end

addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))

%% Load or generate data

% data = generateSimulatedData;

DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';
if ~exist('data','var')
    data = combineData_mat(DataDir); % TODO Switch to csv
end

data.trialdata.nBig = 64;
data.trialdata.nSmall = 4;
data.trialdata.nTrials = 98104;
warning('Assuming nBig, nSmall, and nTrials')


%% Iterative fitting

currentValues.dPrime            = 1;
currentValues.gaussianWidth     = 60;
currentValues.stimulusRemapping = NaN;
currentValues.skewedGaussians   = NaN;

for i = 1:5

    params.dPrime               = true;
    params.gaussianWidth        = false;
    params.stimulusRemapping    = false;
    params.skewedGaussians      = false;

    [currentValues.dPrime,...
        aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,...
        'dPrime',               currentValues.dPrime,...
        'gaussianWidth',        currentValues.gaussianWidth,...
        'stimulusRemappingPol', currentValues.stimulusRemapping,...
        'skewedGaussians',      currentValues.skewedGaussians)

    save(['C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\',...
        num2str(rn),'_',datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


    params.dPrime               = false;
    params.gaussianWidth        = true;
    params.stimulusRemapping    = false;
    params.skewedGaussians      = false;

    [currentValues.gaussianWidth,...
        aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,...
        'dPrime',               currentValues.dPrime,...
        'gaussianWidth',        currentValues.gaussianWidth,...
        'stimulusRemappingPol', currentValues.stimulusRemapping,...
        'skewedGaussians',      currentValues.skewedGaussians)

    save(['C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\',...
        num2str(rn),'_',datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


    params.dPrime               = false;
    params.gaussianWidth        = false;
    params.stimulusRemapping    = true;
    params.skewedGaussians      = false;

    [currentValues.stimulusRemapping,...
        aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,...
        'dPrime',               currentValues.dPrime,...
        'gaussianWidth',        currentValues.gaussianWidth,...
        'stimulusRemappingPol', currentValues.stimulusRemapping,...
        'skewedGaussians',      currentValues.skewedGaussians)

    save(['C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\',...
        num2str(rn),'_',datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


    params.dPrime               = false;
    params.gaussianWidth        = false;
    params.stimulusRemapping    = false;
    params.skewedGaussians      = true;

    [currentValues.skewedGaussians,...
        aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,...
        'dPrime',               currentValues.dPrime,...
        'gaussianWidth',        currentValues.gaussianWidth,...
        'stimulusRemappingPol', currentValues.stimulusRemapping,...
        'skewedGaussians',      currentValues.skewedGaussians)

    save(['C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\',...
        num2str(rn),'_',datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

end

%% Add in stimulus remapping

% load('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\ssnu\4_230803-161235.mat',...
%     'x')
%
% params.dPrime           = false;
% params.gaussianWidth    = false;
%
% [x,aic,bic] = ParameterEstimator(data,params,rn,'stimulusRemappingPol', x);

%%

% data.trialdata.dPrime           = x(1);
% data.trialdata.gaussianWidth    = x(2);

data.trialdata.dPrime           = 1.4899;
data.trialdata.gaussianWidth    = 39.0070;

params.dPrime           = false;
params.gaussianWidth    = false;

%

[x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);

% disp(x)
disp(aic)
disp(nll_x)

%%

% TODO [#C] automate the switch in naming so that they go into the right folders
save(['C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses\TCCModels\',...
    'ssnu_THEN_sg\',...
    num2str(rn),'_',datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


end