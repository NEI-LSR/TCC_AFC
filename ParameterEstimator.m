function [x,aic,bic] = ParameterEstimator(data,params,rn)

if ~exist('bads.m','file')
    warning('BADS is not installed. Check that you have run `submodule init`.')
end

nBig = data.trialdata.nBig;

%% Set bounds (lb/ub) and starting values (x0)

if isfield(params,'dPrime')
    lb = 0;
    ub = 10;
    x0 = lb + (ub-lb).*rand(1,numel(lb));
end

if isfield(params,'skewedGaussians')
    lb = zeros(1,nBig);
    ub = ones(1,nBig);
    x0 = rand(nBig,1);
end

if isfield(params,'stimulusRemapping') % Polar 
    lb = zeros(1,nBig) + 0.01;
    ub = (ones(1,nBig) *  100) + rand(1,nBig);
    x0 = ones(1,nBig) + rand(1,nBig);
end

if isfield(params,'stimulusRemappingCart') % Cartesian
    lb = (ones(1,nBig*2) * -35) + randn(1,nBig*2);
    ub = (ones(1,nBig*2) *  35) + randn(1,nBig*2);
    x0 = zeros(1,nBig*2) + randn(1,nBig*2);
end

if isfield(params,'freeSimilarityMatrix')
    lb = zeros(1,nBig^2);
    ub = ones(1,nBig^2);
    x0 = lb + (ub-lb).*rand(1,numel(lb));
end

% % multi % %

if isfield(params,'gaussianWidth') && isfield(params,'skewedGaussians')
    lb = [20, zeros(1,nBig)];
    ub = [100, ones(1,nBig)];
    % x0 = ones(1,nBig) * 0.5;
    x0 = [50,rand(1,nBig)];
end

if isfield(params,'dPrime') && isfield(params,'gaussianWidth')
    lb = [0.1, 1];
    ub = [50, 200];
    x0 = [2, 50];
end

if isfield(params,'stimulusRemapping') && isfield(params,'skewedGaussians')
    lb = zeros(1,nBig*2) + 0.001;
    ub = ones(1,nBig*2) + rand(1,nBig*2)/100;
    x0 = ones(1,nBig*2)*0.5 + (rand(1,nBig*2)-0.5);
end

%% What do the x values represent in the fitting?

% Which parameters?
optimisationMeta(:,1) = [...
    isfield(params,'freeSimilarityMatrix');...
    isfield(params,'dPrime');...
    isfield(params,'stimulusRemappingCart');...
    isfield(params,'stimulusRemapping');...
    isfield(params,'gaussianWidth');...
    isfield(params,'skewedGaussians'),...
    ];
optimisationMeta(:,1) = logical(optimisationMeta(:,1));

disp(optimisationMeta)

% How many values does each paramter have?
optimisationMeta(:,2) = [...
    nBig*nBig;...   % free similarityMatrix
    1;...           % dPrime
    nBig*2;...      % stimulus remapping (cartesian)
    nBig;...        % stimulus remapping (polar, angle in degrees)
    1;...           % gaussianWidth
    nBig,...        % skewedGaussians
    ];

%% Hpyerparameters

if isfield(params,'stimulusRemapping')
    options = optimoptions('lsqnonlin',...
        'Display','iter-detailed',...
        'FunctionTolerance', 1e-50,...
        'StepTolerance', 1e-50,...
        'PlotFcn',@optimplotx,...
        'DiffMinChange',0.1);
end

if isfield(params,'freeSimilarityMatrix')
    options = optimoptions('lsqnonlin',...
        'Display','iter-detailed',...
        'PlotFcn',{@optimplotx_SimilarityMatrix},...
        'FunctionTolerance', 10e-10,...
        'MaxFunctionEvaluations',1e6,...
        'UseParallel',true);
end

f = @(x)GenerativeModel(x,... % anonymous function so that we can pass additional parameters
    'choiceInds',   cell2mat(data.trialdata.choices)',...
    'cueInd',       cell2mat(data.trialdata.cues),...
    'response',     cell2mat(data.trialdata.chosen)',...
    'nTrials',      data.trialdata.nTrials,...
    'nBig',         data.trialdata.nBig, ...
    'nSmall',       data.trialdata.nSmall,...    % 'dPrime',       data.trialdata.dPrime,...    % 'gaussianWidth',data.trialdata.gaussianWidth,...
    'optimisationMeta',optimisationMeta); % (add stimCols to speed up slightly)

tic % timing test
nll0 = f(x0);
disp(nll0)
toc

tic % second timing test, because the first is always slow
nll0 = f(x0);
disp(nll0)
toc

if isfield(params,'freeSimilarityMatrix')
    [~,data_x0] = f(x0);

    if exist('plotSimilarityMatrix','file')
        plotSimilarityMatrix(data_x0.trialdata.similarityMatrix)
    else
        figure,
        imagesc(data_x0.trialdata.similarityMatrix)
        axis equal tight
        colormap('gray')
        colorbar
        caxis([0 1])
        xlabel('Choice')
        ylabel('Cue')
    end
    title('x0 Similarity Matrix')
end


%% Recovery - Run optimizer

rng(rn);

if length(x0) < 3 % use bads for low dimensional models
    x = bads(f,x0,lb,ub);
else
    x = lsqnonlin(f,x0,lb,ub,options);
    % [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqnonlin(f,x0,lb,ub,options);
    % CI = nlparci(x,FVAL,'jacobian',JACOB);
end

[nll_x,data_x] = f(x);

%% Plotting results

% figure, hold on
% scatter(skewedGaussians,x,'k*');
% axis equal square
% plot([0,1],[0,1],'k')

% [nll_x,data_x] = f(x);
% 
% figure,
% imagesc(data_x.trialdata.similarityMatrix)
% axis equal tight
% colormap('gray')
% colorbar
% caxis([0 1])
% xlabel('Choice')
% ylabel('Cue')
% title('Recovered matrix')



% disp(gaussianWidth);
% disp(x(1));

% % angles
% x0_norm = x0(1:nBig)*(360/sum(x0(1:nBig)));
% x_norm  =  x(1:nBig)*(360 /sum(x(1:nBig)));
% 
% figure, hold on, axis tight
% plot(pol(1,:),'DisplayName','Ground truth')
% plot(cumsum([0,x0_norm(1:end-1)]),'DisplayName','x0')
% plot(cumsum([0,x_norm(1:end-1)]),'DisplayName','x')
% legend('Location','best')
% xticks([min(xlim),max(xlim)])
% yticks([min(ylim),max(ylim)])


% % attractors
% figure, hold on, 
% legend('Location','best')
% axis tight
% plot(1:nBig, x(1:nBig),'DisplayName','x')
% % plot...
% 
% figure, hold on, 
% legend('Location','best')
% axis tight
% plot((1:nAttactors)+nBig, attractorWeights,'DisplayName','attractorWeights')
% plot((1:nAttactors)+nBig, x((1:nAttactors)+nBig),'DisplayName','x')
% plot((1:nAttactors)+nBig, x0((1:nAttactors)+nBig),':','DisplayName','x0')


% stimuli scatter
% [stimCols,pol] = generateStimCols('nBig',nBig);
% [stimCols_x0(:,1),stimCols_x0(:,2)] = pol2cart(deg2rad(cumsum([0,x0_norm(1:end-1)])), pol(2,:));
% [stimCols_x(:,1),stimCols_x(:,2)]   = pol2cart(deg2rad(cumsum([0, x_norm(1:end-1)])), pol(2,:));
% 
% figure,
% plot3(stimCols_x0(:,1),stimCols_x0(:,2),1:nBig,'ko-','DisplayName','x_0')
% axis equal, legend('Location','best')
% view(2)
% 
% figure,
% plot3(stimCols_x(:,1),stimCols_x(:,2),1:nBig,'ko-','DisplayName','x_{lsq}')
% axis equal, legend('Location','best')
% view(2)


% % similarity matrix
% [nll_x,data2] = f(x);
% figure,
% imagesc(data2.trialdata.similarityMatrix)
% axis equal tight
% colormap('gray')
% colorbar
% caxis([0 1])
% xlabel('Choice')
% ylabel('Cue')

%% Compute AIC/BIC

numParam = length(x);
warning('Reminder: Do you need to add any other parameters?')

[aic,bic] = aicbic(-nll_x,numParam,data.trialdata.nTrials);

end
