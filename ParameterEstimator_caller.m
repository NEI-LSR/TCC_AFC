function [x,aic,bic,nll_x,x0] = ParameterEstimator_caller(rn,data,fittingType)

SaveDir = '.'; % Location for models to be saved out 

if ~exist('rn','var')
    rn = 0;
    warning('Using default value for `rn` (0)')
end

% addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories'))

%% Load or generate data

% data = generateSimulatedData;

% DataDir = 'C:\Users\cege-user\Documents\MacaqueColorCategories\Data\';
% if ~exist('data','var')
%     data = combineData_mat(DataDir); % TODO Switch to csv
% end
%
% data.trialdata.nBig = 64;
% data.trialdata.nSmall = 4;
% data.trialdata.nTrials = 98104;
% warning('Assuming nBig, nSmall, and nTrials')

%% Overwrite `data` with simulated responses, for testing purposes

% remove `response` to get new data,
% given fixed parameters
% aka, just use the cues/choices etc
%
% rng(rn)
%
% stimulusRemapping   = rand(1,data.trialdata.nBig);
% skewedGaussians     = rand(1,data.trialdata.nBig);
%
% choiceInds  = cell2mat(data.trialdata.choices);
% cueInd      = cell2mat(data.trialdata.cues);
%
% [nll_data,data] = GenerativeModel([],... % anonymous function so that we can pass additional parameters
%     'choiceInds',       choiceInds',...
%     'cueInd',           cueInd,...
%     'nTrials',          data.trialdata.nTrials,...
%     'nBig',             data.trialdata.nBig, ...
%     'nSmall',           data.trialdata.nSmall,...
%     'dPrime',               1,...
%     'gaussianWidth',        60,...
%     'stimulusRemappingPol', stimulusRemapping,...
%     'skewedGaussians',      skewedGaussians);
%
% plotSimilarityMatrix(data.trialdata.similarityMatrix)
%
% data.trialdata.nBig = 64;
% data.trialdata.nSmall = 4;
% data.trialdata.nTrials = 98104;
% warning('Assuming nBig, nSmall, and nTrials')

%% Key:

% 1 - Free Similarity Matrix
% 2 - dPrime
% 3 - stimulus remapping (cartesian)
% 4 - stimulus remapping (polar, angle in degrees)
% 5 - gaussianWidth
% 6 - skewedGaussians
% 7 - offsetGaussian

%% dPrimeAndGaussianWidth

if strcmp(fittingType,'dPrimeAndGaussianWidth')

    params = [0,1,0,0,1,0,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);

    save([SaveDir, filesep,...
        fittingType,num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

end

%% Single fitting - ssnu

if strcmp(fittingType,'single-ssnu')

    params = [0,1,0,0,1,0,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);

    % save([SaveDir, filesep,...
    %     fittingType,num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    %     '-regexp', '^(?!(data)$).')
    % % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

    disp(x(1))
    disp(x(2))

    params = [0,0,0,1,0,0,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
        'dPrime',           x(1),...
        'gaussianWidth',    x(2));

    save([SaveDir, filesep,...
        fittingType,num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

end

%% Single fitting - sg

if strcmp(fittingType,'single-sg')

    params = [0,1,0,0,1,0,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);

    % save([SaveDir, filesep,...
    %     fittingType,num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    %     '-regexp', '^(?!(data)$).')
    % % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

    params = [0,0,0,0,0,1,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
        'dPrime',           x(1),...
        'gaussianWidth',    x(2));

    save([SaveDir, filesep,...
        fittingType,num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

end

%% Single fitting - og

if strcmp(fittingType,'single-og')

    params = [0,1,0,0,1,0,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);

    % save([SaveDir, filesep,...
    %     fittingType,num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
    %     '-regexp', '^(?!(data)$).')
    % % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

    disp(x(1))
    disp(x(2))
    
    params = [0,0,0,0,0,0,1];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
        'dPrime',           x(1),...
        'gaussianWidth',    x(2));

    save([SaveDir, filesep,...
        fittingType,num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

end

%% GaussianWidthsAndSkewedGaussians

if strcmp(fittingType,'GaussianWidthsAndSkewedGaussians')

    params = [0,0,0,0,1,1,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,[],[],...
        'dPrime',2);

    save([SaveDir, filesep,...
        fittingType,num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

end

%% Simultaneous fitting

if strcmp(fittingType,'simultaneous')

    params = [0,1,0,1,1,1,0];

    [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn);

    save([SaveDir, filesep,...
        fittingType,num2str(rn),datestr(now,'yymmdd-HHMMSS'),'.mat'],...
        '-regexp', '^(?!(data)$).')
    % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

end

%% Iterative fitting

if strcmp(fittingType,'iterative')

    rng(rn);
    % Here we want the result to be reproducible, but we want each loop to
    % be different, and so we call the rng _here_ instead of passing it to
    % `ParamaterEstimator`

    % Initialise starting values
    currentValues.dPrime            = 1;
    currentValues.gaussianWidth     = 60;
    currentValues.stimulusRemapping = NaN; % NaN defaults to no stimulus remapping
    currentValues.skewedGaussians   = NaN; % NaN defaults to no skew in the gaussians

    for i = 1:3

        params = [0,1,0,0,0,0,0]; % dPrime

        [currentValues.dPrime,...
            aic,bic,nll_x,x0] = ParameterEstimator(data,params,randi(100,1),[],[],...
            'dPrime',               currentValues.dPrime,...
            'gaussianWidth',        currentValues.gaussianWidth,...
            'stimulusRemappingPol', currentValues.stimulusRemapping,...
            'skewedGaussians',      currentValues.skewedGaussians);

        save([SaveDir, filesep,...
            fittingType, datestr(now,'yymmdd-HHMMSS'),'.mat'],...
            '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


        params = [0,0,0,0,1,0,0]; % gaussianWidth

        [currentValues.gaussianWidth,...
            aic,bic,nll_x,x0] = ParameterEstimator(data,params,randi(10,1),[],[],...
            'dPrime',               currentValues.dPrime,...
            'gaussianWidth',        currentValues.gaussianWidth,...
            'stimulusRemappingPol', currentValues.stimulusRemapping,...
            'skewedGaussians',      currentValues.skewedGaussians);

        save([SaveDir, filesep,...
            fittingType, datestr(now,'yymmdd-HHMMSS'),'.mat'],...
            '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


        params = [0,0,0,1,0,0,0]; % gaussianWidth

        [currentValues.stimulusRemapping,...
            aic,bic,nll_x,x0] = ParameterEstimator(data,params,randi(10,1),[],[],...
            'dPrime',               currentValues.dPrime,...
            'gaussianWidth',        currentValues.gaussianWidth,...
            'stimulusRemappingPol', currentValues.stimulusRemapping,...
            'skewedGaussians',      currentValues.skewedGaussians);

        save([SaveDir, filesep,...
            fittingType, datestr(now,'yymmdd-HHMMSS'),'.mat'],...
            '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)


        params = [0,0,0,0,0,1,0]; % gaussianWidth

        [currentValues.skewedGaussians,...
            aic,bic,nll_x,x0] = ParameterEstimator(data,params,randi(10,1),[],[],...
            'dPrime',               currentValues.dPrime,...
            'gaussianWidth',        currentValues.gaussianWidth,...
            'stimulusRemappingPol', currentValues.stimulusRemapping,...
            'skewedGaussians',      currentValues.skewedGaussians);

        save([SaveDir, filesep,...
            fittingType, datestr(now,'yymmdd-HHMMSS'),'.mat'],...
            '-regexp', '^(?!(data)$).') % save everything except data (https://www.mathworks.com/matlabcentral/answers/101287-how-do-i-save-all-of-the-workspace-variables-except-for-a-certain-specified-variable-name-in-matlab#answer_110635)

    end

    load('iterative230817-153048.mat')
    x = [currentValues.dPrime,...
        currentValues.stimulusRemapping,...
        currentValues.gaussianWidth,...
        currentValues.skewedGaussians'];
end

end