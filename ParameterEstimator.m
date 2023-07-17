clear, 
clc, 
close all

addpath(genpath('C:\toolbox\bads'))
addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses'))

%%
simdata_or_realdata = 'real'; % 'sim' to simulate data, 'real' to load real data

if strcmp(simdata_or_realdata,'sim') % simuldated data
    
    nTrials = 50000;
    nBig = 64;
    nSmall = 4;
   
    % % Create offsets (simulating an inhomogenous colorspace)
    % (copied from `justColSpace.m`)
    
    % [stimCols,pol] = generateStimCols('nBig',nBig);
    % interval = 360/nBig;
    % magnitude = 50;
    % offset = sin(deg2rad(0:interval:360-interval)) * magnitude;
    % 
    % pol(1,:) = pol(1,:) + offset;
    % [cart(:,1),cart(:,2)] = pol2cart(deg2rad(pol(1,:)),pol(2,:));
    % 
    % stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols(1,:);stimCols(2,:)]);
    % colvals = im2double(stimCols_sRGB);
    % 
    % figure,
    % scatter(cart(:,1),cart(:,2),75,colvals,'filled')
    % axis equal
    % axis off
    % saveas(gcf,fullfile(['stimulus-space-non-uniformity_', datestr(now,'yymmdd'), '.svg']))
    % 
    % figure,
    % plot([cart(:,1);cart(:,2)]' - [stimCols(1,:),stimCols(2,:)])

    % - % Attractor Points

%     stimCols_ap = generateStimCols('nBig',nBig);

%     attractorPoints =  [stimCols_ap(1,:),stimCols_ap(2,:)]; %[0,-38,38,0];
%     attractorWeights = zeros(1,nBig);%[20,30];
%     attractorWeights(16) = 20;
%     attractorWeights(40) = 30;

    % attractorPoints =  [0,-38,38,0];
    % attractorWeights = [20,30];

    % - %
    
    
    % % Generate simulated data

    rng(0);
      
%     [nll, data] = GenerativeModel([],...
%         'attractorPoints',attractorPoints,...
%         'attractorWeights',attractorWeights,...
%         'stimulusRemappingCart',[cart(:,1);cart(:,2)]' - [stimCols(1,:),stimCols(2,:)],...
%         'nTrials',nTrials,...
%         'nBig',nBig,...
%         'nSmall',nSmall);

    % r = [0.1, 0.9];
    % skewedGaussians = ((sin(deg2rad(linspace(0,360,nBig)))+1)/2);
    % skewedGaussians = skewedGaussians/(1/diff(r))+r(1);
    % 
    % skewedGaussians = linspace(0.1,0.9,nBig); % linear skew
    % skewedGaussians = ones(1,nBig)*0.5; % no skew
    
    % figure, plot(skewedGaussians); axis tight % arbitrary, just for testing

    dprime = 1 + 2/3;
    SimFunc_sd = 100*(2/3);

    [nll, data] = GenerativeModel([],...
        'dprime',dprime,...
        'SimFunc_sd', SimFunc_sd,...
        'nTrials',nTrials,...
        'nBig',nBig,...
        'nSmall',nSmall,...
        'pltSimFigs',true);

            % 
    
    choiceInds =    cell2mat(data.trialdata.choices)';
    cueInd =        cell2mat(data.trialdata.cues);
    response =      cell2mat(data.trialdata.chosen);
    stimCols =      data.trialdata.stimCols;
    
elseif strcmp(simdata_or_realdata,'real') % real data

    %     filepath = 'Y:\PROJECTS\CausalGlobs\data\';
    filepath = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\data\'; %!!!!!!!!!!!!!!!
    %filename = '210422--211012_Pollux_data';
    %filename = '210517--211108_Castor_data.mat';
    %filename = '220322--220823_Morty_data.mat';
    %filename = '210428--210609_Buster_data';
    filename = 'combinedData.mat';
    
    load([filepath,filename])
    try % for individual data
        data.trialdata = trialdata;
    catch % for combined data
        data.trialdata = cleandata.trialdata;
    end

    choiceInds =    cell2mat(data.trialdata.choices)';
    cueInd =        cell2mat(data.trialdata.cues);
    response =      cell2mat(data.trialdata.chosen);

    aborts = isnan(cell2mat(data.trialdata.chosen));
    
    [~,col] = ind2sub(size(choiceInds),find(isnan(choiceInds(:))));
    trailsWithNaNChoices = false(size(aborts));
    trailsWithNaNChoices(unique(col)) = 1;

    choiceInds =    choiceInds(:,~aborts & ~trailsWithNaNChoices);
    cueInd =        cueInd(~aborts & ~trailsWithNaNChoices);
    response =      response(~aborts & ~trailsWithNaNChoices);

    %nBig = size(data.trialdata.stimCols{1,1},1);
    %nSmall = size(data.trialdata.choices{1,1},2);
    nBig = length(unique(cell2mat(data.trialdata.cues)));
    %nSmall = max(unique(cell2mat(data.trialdata.chosen_idx)));
    nSmall = length(data.trialdata.choices{1,1});
    nTrials = sum(~aborts & ~trailsWithNaNChoices);
    
    % add stimCols to speed up slightly

%     stimCols_ap = generateStimCols('nBig',nBig,'showFig',true);
%     attractorPoints_inds = [3,38];
%     attractorPoints =  [stimCols_ap(1,attractorPoints_inds),stimCols_ap(2,attractorPoints_inds)];

end

%% Recovery - Prepare optimizer

rng('shuffle')

nAttactors = 0;

% % --- for stimulus remapping --- % %
% lb = (ones(1,nBig*2) * -35) + randn(1,nBig*2); % lower bound
% ub = (ones(1,nBig*2) *  35) + randn(1,nBig*2); % upper bound 
% x0 = zeros(1,nBig*2) + randn(1,nBig*2);
% % --- % %

% % --- for free similarityMatrix --- % %
% lb = zeros(1,nBig^2); % lower bound
% ub = ones(1,nBig^2); % upper bound                      
% x0 = lb + (ub-lb).*rand(1,numel(lb)); % random points between lb and ub
% % --- % %

% % % --- for attractor dynamics --- % %
% % nAttactors = 2;
% % lb = [-50,-50,-50,-50,0,0]; % lower bound
% % ub = [50,50,50,50,80,80]; % upper bound                      
% % x0 = lb + (ub-lb).*rand(1,numel(lb)); % random points between lb and ub
% nAttactors = nBig;
% lb = [zeros(1,nBig)]; % lower bound
% ub = [ones(1,nBig)*30]; % upper bound                      
% x0 = lb + (ub-lb).*rand(1,numel(lb))/2; % random points between lb and ub
% % % --- % %

% % ----- for d prime ----- % %
% lb = 0; % lower bound
% ub = 10; % upper bound 
% x0 = lb + (ub-lb).*rand(1,numel(lb)); % random points between lb and ub
% % --- % %

% % % stimulus remapping (cartesian) AND attractor dynamics
% nAttactors = nBig;
% lb = [(ones(1,nBig*2) * -35) + randn(1,nBig*2), zeros(1,nBig)]; % lower bound
% ub = [(ones(1,nBig*2) *  35) + randn(1,nBig*2), ones(1,nBig)*30]; % upper bound 
% x0 = [zeros(1,nBig*2) + randn(1,nBig*2), lb((nBig*2)+1:end) + (ub((nBig*2)+1:end)-lb((nBig*2)+1:end)).*rand(1,numel(lb((nBig*2)+1:end)))/2];
% % %

% stimulus remapping (angle) 
lb = zeros(1,nBig) + 0.01; % lower bound
ub = (ones(1,nBig) *  100) + rand(1,nBig); % upper bound 
x0 = ones(1,nBig) + rand(1,nBig);
%

% % % stimulus remapping (angle) AND attractor dynamics
% nAttactors = nBig;
% lb = [zeros(1,nBig) + 0.01, zeros(1,nBig)]; % lower bound
% ub = [(ones(1,nBig) *  100) + rand(1,nBig), ones(1,nBig)*30]; % upper bound 
% x0 = [ones(1,nBig) + rand(1,nBig), lb((nBig)+1:end) + (ub((nBig)+1:end)-lb((nBig)+1:end)).*rand(1,numel(lb((nBig)+1:end)))/2];
% % %

% % % ----- for d prime AND lambda/sigma----- % %
% lb = [0.1,-1,0.1]; % lower bound
% ub = [7,0,10]; % upper bound 
% x0 = lb + (ub-lb).*rand(1,numel(lb)); % random points between lb and ub
% % % --- % %

% % % ----- for d prime AND lambda/sigma AND attractors ----- % %
% lb = [0.1,0,0,-1,0.1]; % lower bound
% ub = [7,50,50,0,10]; % upper bound 
% x0 = lb + (ub-lb).*rand(1,numel(lb)); % random points between lb and ub
% % % --- % %

% % % ----- for d prime AND lambda/sigma AND stimulus remapping ----- % %
% lb = [0.1,zeros(1,nBig) + 0.01,-1,0.1]; % lower bound
% ub = [7,ones(1,nBig) *  100,0,10]; % upper bound 
% x0 = lb + (ub-lb).*rand(1,numel(lb)); % random points between lb and ub
% 
% x0(1) = 2.443912865562119;
% x0(66) = -0.015041279678144;
% x0(67) = 3.913215487407419;
% x0(2:nBig+1) = ones(1,nBig) + rand(1,nBig);
% 
% % % --- % %

% % % --- skewed gaussians ---
% lb = zeros(1,nBig); % lower bound
% ub = ones(1,nBig); % upper bound 
% x0 = rand(nBig,1);
% %

% % % --- SimFunc_sd and skewed gaussians ---
% lb = [20, zeros(1,nBig)]; % lower bound
% ub = [100, ones(1,nBig)]; % upper bound 
% % x0 = ones(1,nBig) * 0.5;
% x0 = [50,rand(1,nBig)];
% % %

% % % d-prime and SimFunc_sd % %
% lb = [0.1, 1]; % lower bound
% ub = [50, 200]; % upper bound 
% x0 = [2, 50];

optimisationMeta = [...      % what do the x values represent in the fitting? (1st column - is it passed, 2nd column - how many values would it be if it were)
    false,   nBig*nBig;...   % free similarityMatrix
    false,   1;...           % dprime
    false,   nBig*2;...      % stimulus remapping (cartesian)
    true,   nBig;...        % stimulus remapping (polar, angle in degrees)
    false,   nAttactors*2;... % attractor points
    false,   nAttactors;...   % attractor weights
    false,    1;...           % SimFunc_sd
    false,    nBig,...        % skewedGaussians              
    ];

% options = optimoptions('lsqnonlin',...
%     'Display','iter-detailed',...
%      'UseParallel',true,...
%      'FunctionTolerance', 1e-15,...
%      'StepTolerance', 1e-15,...
%      'DiffMinChange',1,...
%      'PlotFcn',@optimplotx_SimilarityMatrix);
%  %    'StepTolerance',1e-1000,...
%  %     'FunctionTolerance', 1e-10
%  %     ...

% options = optimoptions('lsqnonlin',...
%     'Display','iter-detailed',...
%      'FunctionTolerance', 1e-50,...
%      'StepTolerance', 1e-50,...
%      'DiffMinChange',1);%,...
   %  'PlotFcn',@optimplotx_SimilarityMatrix);
 %    'StepTolerance',1e-1000,...
 %     'FunctionTolerance', 1e-10
 %     ...

options = optimoptions('lsqnonlin',...
    'Display','iter-detailed',...
    'FunctionTolerance', 1e-50,...
     'StepTolerance', 1e-50,...
     'PlotFcn',@optimplotx,...
     'DiffMinChange',1);

% options = optimoptions('lsqnonlin',...
%     'Display','iter-detailed',...
%     'PlotFcn',@optimplotx);

if optimisationMeta(1,1)
    options = optimoptions('lsqnonlin',...
        'Display','iter-detailed',...
        'PlotFcn',{@optimplotx_SimilarityMatrix},...
        'FunctionTolerance', 10e-10,...
        'MaxFunctionEvaluations',1e6,...
        'UseParallel',true);
end

 
% f = @(x)GenerativeModel(x,... % anonymous function so that we can pass additional parameters
%     'choiceInds',choiceInds,...
%     'cueInd',cueInd,...
%     'response',response,...
%     'nTrials',nTrials,...
%     'nBig',nBig, ...
%     'nSmall',nSmall,...
%     'attractorPoints',attractorPoints,...
%     'optimisationMeta',optimisationMeta); % add stimCols to speed up slightly

f = @(x)GenerativeModel(x,... % anonymous function so that we can pass additional parameters
    'choiceInds',choiceInds,...
    'cueInd',cueInd,...
    'response',response,...
    'nTrials',nTrials,...
    'nBig',nBig, ...
    'nSmall',nSmall,...
    'dprime',1.4977,...
    'SimFunc_sd',39.1110,...
    'optimisationMeta',optimisationMeta); % add stimCols to speed up slightly

tic % timing test
nll0 = f(x0);
toc

tic % second timing test, because the first is always slow
nll0 = f(x0);
toc

[~,data_x0] = f(x0);

figure,
imagesc(data_x0.trialdata.similarityMatrix)
axis equal tight
colormap('gray')
colorbar
caxis([0 1])
xlabel('Choice')
ylabel('Cue')
title('x0 matrix')


%% Recovery - Run optimizer

x = lsqnonlin(f,x0,lb,ub,options);
% [x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqnonlin(f,x0,lb,ub,options);

% CI = nlparci(x,FVAL,'jacobian',JACOB);
% 
% x = bads(f,x0,lb,ub);


%% Plotting results

% figure, hold on
% scatter(skewedGaussians,x,'k*');
% axis equal square
% plot([0,1],[0,1],'k')

[nll_x,data_x] = f(x);

figure,
imagesc(data_x.trialdata.similarityMatrix)
axis equal tight
colormap('gray')
colorbar
caxis([0 1])
xlabel('Choice')
ylabel('Cue')
title('Recovered matrix')





% disp(SimFunc_sd);
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

%%




%% similarity matrix
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

% numParam = length(x);
numParam = length(x)+2;
warning('Reminder: Do you need to add any other parameters?')

[aic,bic] = aicbic(-nll_x,numParam,nTrials);
% [aic,bic] = aicbic(-nll_x,numParam,nTrials,Normalize=true);
