function [x,aic,bic,nll_x,x0] = ParameterEstimator(data,params,rn,dim1,dim2,varargin)

if ~exist('bads.m','file')
    warning('BADS is not installed. Check that you have run `submodule init`.')
end

% Key:
% 1 - Free Similarity Matrix
% 2 - dPrime
% 3 - stimulus remapping (cartesian)
% 4 - stimulus remapping (polar, angle in degrees)
% 5 - gaussianWidth
% 6 - skewedGaussians
% 7 - offsetGaussians

if ~isfield(data.trialdata,'nBig')
    data.trialdata.nBig = 64;
    warning('Assuming nBig == 64')
end
if ~isfield(data.trialdata,'nSmall')
    data.trialdata.nSmall = 4;
    warning('Assuming nSmall == 4')
end
if ~isfield(data.trialdata,'nTrials')
    data.trialdata.nTrials = 98104;
    warning('Assuming nTrials == 98104')
end

nBig = data.trialdata.nBig;

%% Set bounds (lb/ub) and starting values (x0)

rng(rn);

if sum(params) == 1

    if params(1) % Free Similarity Matrix
        lb = zeros(1,nBig^2);
        ub = ones(1,nBig^2);
        x0 = lb + (ub-lb).*rand(1,numel(lb));
    end

    if params(2) % dPrime
        lb = 0.1;
        ub = 100;
        x0 = lb + (ub-lb).*rand(1,numel(lb));
    end

    if params(3) % stimulus remapping (cartesian)
        lb = (ones(1,nBig*2) * -35) + randn(1,nBig*2);
        ub = (ones(1,nBig*2) *  35) + randn(1,nBig*2);
        x0 = zeros(1,nBig*2) + randn(1,nBig*2);
    end

    if params(4) % stimulus remapping (polar, angle in degrees)
        lb = zeros(1,nBig) + 0.01;
        ub = (ones(1,nBig) *  100) + rand(1,nBig);
        x0 = ones(1,nBig) + rand(1,nBig);
    end

    if exist('dim1','var') && ~isempty(dim1) % TODO (hacky, temporary)
        if params(4) % stimulus remapping (polar, angle in degrees)
            lb = zeros(1,dim1) + 0.01;
            ub = (ones(1,dim1) *  100) + rand(1,dim1);
            x0 = ones(1,dim1) + rand(1,dim1);
        end
    end

    if params(5) % gaussianWidth
        lb = 1;
        ub = 100;
        x0 = lb + (ub-lb).*rand(1,numel(lb));
    end

    if params(6) % skewedGaussians
        lb = zeros(1,nBig);
        ub = ones(1,nBig);
        x0 = rand(1,nBig);
    end
    
    if params(7) % offsetGaussians
        lb = ones(1,nBig)*-180;
        ub = ones(1,nBig)*180;
        x0 = (rand(1,nBig)*30)-15; % somewhat arbitrary
    end

    if exist('dim2','var') && ~isempty(dim2) % TODO (hacky, temporary)
        if params(6) % skewedGaussians
            lb = zeros(1,dim2);
            ub = ones(1,dim2);
            x0 = rand(1,dim2);
        end
        
        if params(7) % skewedGaussians
            lb = ones(1,dim2)*-180;
            ub = ones(1,dim2)*180;
            x0 = (rand(1,dim2)*30)-15; % somewhat arbitrary
        end
    end

elseif sum(params) == 2

    % % multi % %

    if params(2) && params(5) % dPrime && gaussianWidth
        lb = [0.1, 1];
        ub = [10, 100];
        x0 = [2, 50];
    end

    % if params(5) && params(6) % gaussianWidth && skewedGaussians
    %     lb = [20, zeros(1,nBig)];
    %     ub = [100, ones(1,nBig)];
    %     % x0 = ones(1,nBig) * 0.5;
    %     x0 = [50,rand(1,nBig)];
    % end

    if params(5) && params(6) % gaussianWidth && skewedGaussians
        lb = [ones(1,nBig), zeros(1,nBig)];
        ub = [ones(1,nBig)*100, ones(1,nBig)];
        % x0 = ones(1,nBig) * 0.5;
        x0 = [(ones(1,nBig)*37) + rand(1,nBig)*6, rand(1,nBig)];
    end

    if params(4) && params(6) % stimulus remapping (polar, angle in degrees) && skewedGaussians
        lb = zeros(1,nBig*2) + 0.001;
        ub = ones(1,nBig*2) + rand(1,nBig*2)/100;
        x0 = ones(1,nBig*2)*0.5 + (rand(1,nBig*2)-0.5);
    end
    
    if params(4) && params(7) % stimulus remapping (polar, angle in degrees) && offsetGaussians (assumes stimulus space reduction)
        lb = [zeros(1,dim1) + 0.01, ones(1,dim2)*-180];
        ub = [(ones(1,dim1) *  100) + rand(1,dim1), ones(1,dim2)*180];
        x0 = [ones(1,dim1) + rand(1,dim1),(rand(1,dim2)*30)-15];
    end

elseif params(2) && params(4) && params(5) && params(6) % "ALL"

    lb = [0.1,  zeros(1,nBig) + 0.001,                     1,      zeros(1,nBig) + 0.001];
    ub = [100,   ones(1,nBig) + rand(1,nBig)/100,          100,    ones(1,nBig) + rand(1,nBig)/100];
    x0 = [2,    ones(1,nBig)*0.5 + (rand(1,nBig)-0.5)/5,   50,     ones(1,nBig)*0.5 + (rand(1,nBig)-0.5)/5];

end


%% What do the x values represent in the fitting?

optimisationMeta = double(zeros([7,2]));

% Which parameters?
optimisationMeta(:,1) = params;

% How many values does each paramter have?
% optimisationMeta(:,2) = [...
%     nBig*nBig;...   % Free Similarity Matrix
%     1;...           % dPrime
%     nBig*2;...      % stimulus remapping (cartesian)
%     nBig;...        % stimulus remapping (polar, angle in degrees)
%     1;...           % gaussianWidth
%     nBig,...        % skewedGaussians
%     ];

% % How many values does each paramter have?
% optimisationMeta(:,2) = [...
%     nBig*nBig;...   % Free Similarity Matrix
%     1;...           % dPrime
%     nBig*2;...      % stimulus remapping (cartesian)
%     nBig;...        % stimulus remapping (polar, angle in degrees)
%     nBig;...        % gaussianWidth
%     nBig;...        % skewedGaussians
%     nBig];

% How many values does each paramter have?
optimisationMeta(:,2) = [...
    nBig*nBig;...   % Free Similarity Matrix
    1;...           % dPrime
    nBig*2;...      % stimulus remapping (cartesian)
    nBig;...        % stimulus remapping (polar, angle in degrees)
    1;...           % gaussianWidth
    nBig;...        % skewedGaussians
    nBig];

if exist('dim1','var') && ~isempty(dim1) % TODO (hacky, temporary)
    optimisationMeta(:,2) = [...
    nBig*nBig;...   % Free Similarity Matrix
    1;...           % dPrime
    nBig*2;...      % stimulus remapping (cartesian)
    dim1;...        % stimulus remapping (polar, angle in degrees)
    1;...           % gaussianWidth
    nBig;...        % skewedGaussians
    nBig];
end

if exist('dim2','var') && ~isempty(dim2) % TODO (hacky, temporary)
    optimisationMeta(:,2) = [...
    nBig*nBig;...   % Free Similarity Matrix
    1;...           % dPrime
    nBig*2;...      % stimulus remapping (cartesian)
    nBig;...        % stimulus remapping (polar, angle in degrees)
    1;...           % gaussianWidth
    nBig;...        % skewedGaussians
    dim2];
end

if exist('dim1','var') && ~isempty(dim1) && exist('dim2','var') && ~isempty(dim2) % TODO (hacky, temporary)
    optimisationMeta(:,2) = [...
    nBig*nBig;...   % Free Similarity Matrix
    1;...           % dPrime
    nBig*2;...      % stimulus remapping (cartesian)
    dim1;...        % stimulus remapping (polar, angle in degrees)
    1;...           % gaussianWidth
    nBig;...        % skewedGaussians
    dim2];
end

%% Hyperparameters

% general options
options = optimoptions('lsqnonlin',...
        'Display','iter-detailed',...
        'FunctionTolerance', 1e-50,...
        'StepTolerance', 1e-50,...
        'MaxIterations', 1500,...
        'PlotFcn',@optimplotx);

% specific options
if params(1)
    options = optimoptions('lsqnonlin',...
        'Display','iter-detailed',...
        'PlotFcn',{@optimplotx_SimilarityMatrix},...
        'FunctionTolerance', 10e-10,...
        'MaxFunctionEvaluations',1e6,...
        'UseParallel',true);
end

if params(4) || params(7)
    options = optimoptions('lsqnonlin',...
        'Display','iter-detailed',...
        'FunctionTolerance', 1e-50,...
        'StepTolerance', 1e-50,...
        'PlotFcn',@optimplotx,...
        'DiffMinChange',0.1,... % TODO Write note about why this is neccessary
        'UseParallel',true); 
end

%% Define and test function

if size(data.trialdata.chosen,1) < size(data.trialdata.chosen,2) % TODO This is a temp hack fix for a difference between real and simulated data
    data.trialdata.chosen = data.trialdata.chosen';
    warning('transposing choices');
end

% data cleaning (remove trials with NaN choices) 
% (I have no idea how such trials occured)

if any(any(isnan(cell2mat(data.trialdata.choices)))) % TODO Move this into the initial data import stage (e.g. combineData.m)
    warning('Detected trials with NaN choices. Removing.')
    [a,~] = ind2sub([data.trialdata.nTrials,data.trialdata.nSmall],...
        find(isnan(cell2mat(data.trialdata.choices))));
    try
        data.trialdata.allchoices(unique(a)) = [];
    catch
    end
    try
        data.trialdata.dirname(unique(a)) = [];
    catch
    end
    data.trialdata.paradigm(unique(a)) = [];
    data.trialdata.choices(unique(a)) = [];
    data.trialdata.chosen_idx(unique(a)) = [];
    data.trialdata.cues(unique(a)) = [];
    try
        data.trialdata.stimCols_raw(unique(a)) = [];
    catch
    end
    data.trialdata.chosen(unique(a)) = [];
    data.trialdata.stimCols(unique(a)) = [];
    data.trialdata.nTrials = size(data.trialdata.cues,1);
end

choiceInds = cell2mat(data.trialdata.choices);
cueInd = cell2mat(data.trialdata.cues);
response = cell2mat(data.trialdata.chosen);

f = @(x)GenerativeModel(x,... % anonymous function so that we can pass additional parameters
    'choiceInds',       choiceInds',...
    'cueInd',           cueInd,...
    'response',         response,...
    'nTrials',          data.trialdata.nTrials,...
    'nBig',             data.trialdata.nBig, ...
    'nSmall',           data.trialdata.nSmall,...
    'optimisationMeta', optimisationMeta,...
    varargin{:});

tic % timing test 1
nll0 = f(x0);
toc

tic % timing test 2 (because the first is always slower)
nll0 = f(x0);
toc

[~,data_x0] = f(x0);

if exist('plotSimilarityMatrix','file') % Part of the MCC repo
    plotSimilarityMatrix(data_x0.trialdata.similarityMatrix)
    title('x0 Similarity Matrix')
else
    figure,
    imagesc(data_x0.trialdata.similarityMatrix') % transposed so cue on x-axis
    axis equal tight
    colormap('gray')
    colorbar
    caxis([0 1])
    xlabel('Choice')
    ylabel('Cue')
    title('x0 Similarity Matrix')
end

%% Recovery - Run optimizer

if length(x0) < 10 && ~exist('dim2','var') % use bads for low dimensional models, TODO hacky way to get it to use lsqnonlin for lowdim sg models
    x = bads(f,x0,lb,ub);
else
    x = lsqnonlin(f,x0,lb,ub,options);

    % TODO Work out how to compute confidence intervals
    % [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqnonlin(f,x0,lb,ub,options);
    % CI = nlparci(x,FVAL,'jacobian',JACOB);
end

% sOpt = optimset('fmin_adam');
% sOpt.GradObj = 'off';
% sOpt.MaxFunEvals = 1e4;
% sOpt.Display = 'iter-detailed';
%
% x = fmin_adam(f,x0,[], [], [], [], [], sOpt);

[nll_x,data_x] = f(x);

if exist('plotSimilarityMatrix','file') % Part of the MCC repo
    plotSimilarityMatrix(data_x.trialdata.similarityMatrix)
    title('x Similarity Matrix')
else
    figure,
    imagesc(data_x.trialdata.similarityMatrix') % transposed so cue on x-axis
    axis equal tight
    colormap('gray')
    colorbar
    caxis([0 1])
    xlabel('Choice')
    ylabel('Cue')
    title('x Similarity Matrix')
end

%% Compute AIC/BIC

% 

% numParam = input(['What number of parameters should be used to compute AIC/BIC?',newline]);
%
% [aic,bic] = aicbic(-nll_x,numParam,data.trialdata.nTrials);
aic = NaN;
bic = NaN;

end
