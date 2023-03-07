% clear, 
% clc, 
% close all

%%
simdata_or_realdata = 'real'; % 'sim' to simulate data, 'real' to load real data

if strcmp(simdata_or_realdata,'sim') % simuldated data
    
    nTrials = 10000;
    nBig = 64;
    nSmall = 4;
    
    
    % % Create offsets (simulating an inhomogenous colorspace)
    % (copied from `justColSpace.m`)
    
%     [stimCols,pol] = generateStimCols('nBig',nBig);
%     interval = 360/nBig;
%     magnitude = 50;
%     offset = sin(deg2rad(0:interval:360-interval)) * magnitude;
%     
%     pol(1,:) = pol(1,:) + offset;
%     [cart(:,1),cart(:,2)] = pol2cart(deg2rad(pol(1,:)),pol(2,:));
%     
%     figure,
%     scatter(cart(:,1),cart(:,2))
%     axis square
%     
%     figure,
%     plot([cart(:,1);cart(:,2)]' - [stimCols(1,:),stimCols(2,:)])
    
    
    % % Generate simulated data
    
    [nll, data] = GenerativeModel([],...
        'nTrials',nTrials, ...
        'nBig',nBig, ...
        'nSmall',nSmall);
    
%     [nll, data] = GenerativeModel([],...
%         'nTrials',nTrials,...
%         'nBig',nBig, ...
%         'nSmall',nSmall,...
%         'dprime',1,...
%         'stimulusRemapping',[cart(:,1);cart(:,2)]' - [stimCols(1,:),stimCols(2,:)]);
    
    % attractorPoints =  [37,-37,0,0];
    % attractorWeights = [0.5,0.6];
    %
    % [nll, data] = generateSimulatedData_TCC_parameterized([],...
    %     'nTrials',nTrials,...
    %     'nBig',nBig, ...
    %     'nSmall',nSmall,...
    %     'dprime',1,...
    %     'attractorPoints',attractorPoints,...
    %     'attractorWeights',attractorWeights,...
    %     'stimulusRemapping',[cart(:,1);cart(:,2)]' - [stimCols(1,:),stimCols(2,:)]);
    
    choiceInds =    cell2mat(data.trialdata.choices)';
    cueInd =        cell2mat(data.trialdata.cues);
    response =      cell2mat(data.trialdata.chosen);
    stimCols =      data.trialdata.stimCols;
    
%     figure,
%     imagesc(data.trialdata.similarityMatrix)
%     axis equal tight
%     colormap('gray')
%     colorbar
%     xlabel('Choice')
%     ylabel('Cue')

elseif strcmp(simdata_or_realdata,'real') % real data

    %     filepath = 'Y:\PROJECTS\CausalGlobs\data\';
    filepath = 'C:\Users\cege-user\Dropbox\Documents\MATLAB\CausalGlobs\data\';
    %filename = '210422--211012_Pollux_data';
    filename = '210517--211108_Castor_data.mat';
    %filename = '220322--220823_Morty_data.mat';
    %filename = '210428--210609_Buster_data';
    
    load([filepath,filename])
    data.trialdata = trialdata;

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

end

%%

%heatmap_analysis(data,'savefigs',false);

% figure, 
% er = response-cueInd; % error
% er(er<-nBig/2) = er(er<-nBig/2)+nBig;
% er(er> nBig/2) = er(er> nBig/2)-nBig;
% error_prob = histcounts(er)/nTrials;
% histogram('BinCounts', error_prob,'BinEdges', linspace(-180,180,42), ...
%     'FaceColor',[0.7,0.7,0.7],'EdgeColor','white')
% xlim([-180,180])
% xticks([-180,0,180])
% %ylim([0,inf])
% yticks([0,max(ylim)])
% set(gca, 'FontSize', 16)
% xlabel('Error (color distance)')
% ylabel('Probability')
% box off

%% Recovery - Prepare optimizer

% % --- for stimulus remapping --- % %
% lb = (ones(1,nBig*2) * -35) + randn(1,nBig*2); % lower bound
% ub = (ones(1,nBig*2) *  35) + randn(1,nBig*2); % upper bound 
% x0 = zeros(1,nBig*2) + randn(1,nBig*2);
% % --- % %

% % --- for free similarityMatrix --- % %
lb = zeros(1,nBig^2); % lower bound
ub = ones(1,nBig^2); % upper bound                      
x0 = lb + (ub-lb).*rand(1,numel(lb)); % random points between lb and ub
% % --- % %

% % --- for attractor dynamics --- % %
% nAttactors = 2;
% lb = [0,0,0,0]; % lower bound
% ub = [nBig,nBig,10,10]; % upper bound                      
% x0 = lb + (ub-lb).*rand(1,numel(lb)); % random points between lb and ub
% % --- % %

% % ----- for d prime ----- % %
% lb = 0; % lower bound
% ub = 10; % upper bound 
% x0 = lb + (ub-lb).*rand(1,numel(lb)); % random points between lb and ub
% % --- % %


optimisationMeta = [...     % what do the x values represent in the fitting? (1st column - is it passed, 2nd column - how many values would it be if it were)
    true,   nBig*nBig;...   % free similarityMatrix
    false,   1;...           % dprime
    false,   nBig*2;...      % stimulus remapping
    false,   4;...           % attractor points
    false,   2;...           % attractor weights
    false,   1;...           % lambda
    false,   1;...           % sigma
    ];

options = optimoptions('lsqnonlin',...
    'Display','iter-detailed',...
     'UseParallel',true,...
     'FunctionTolerance', 1e-15,...
     'StepTolerance', 1e-15,...
     'DiffMinChange',1,...
     'PlotFcn',@optimplotx);
 %    'StepTolerance',1e-1000,...
 %     'FunctionTolerance', 1e-10
 %     ...

if optimisationMeta(1,1)
    options = optimoptions('lsqnonlin',...
        'Display','iter-detailed',...
        'PlotFcn',{@optimplotx_SimilarityMatrix},...
        'FunctionTolerance', 10e-10,...
        'MaxFunctionEvaluations',1e6,...
        'UseParallel',true);
end

 
f = @(x)GenerativeModel(x,... % anonymous function so that we can pass additional parameters
    'choiceInds',choiceInds,...
    'cueInd',cueInd,...
    'response',response,...
    'nTrials',nTrials,...
    'nBig',nBig, ...
    'nSmall',nSmall,...
    'optimisationMeta',optimisationMeta); % add stimCols to speed up slightly

tic % timing test
nll0 = f(x0);
toc

tic % second timing test, because the first is always super slow
nll0 = f(x0);
toc


%% Recovery - Run optimizer

[x,Resnorm,FVAL,EXITFLAG,OUTPUT,LAMBDA,JACOB] = lsqnonlin(f,x0,lb,ub,options);

%x = bads(f,x0,lb,ub,[],PUB);


%%

% disp(x0)
% disp(x)

%%

% figure,
% scatter(cart(:,1) - stimCols(1,:)',x(2:nBig+1))
% 
% figure,
% scatter(cart(:,2) - stimCols(2,:)',x(nBig+2:nBig*2+1))

% figure, hold on
% axis equal tight
% legend('Location','best')
% scatter(cart(:,1),cart(:,2),'DisplayName','true')
% scatter(stimCols(1,:)+x0(1:nBig),stimCols(2,:)+x0(nBig+1:end),'filled','DisplayName','x0')
% scatter(stimCols(1,:)+x(1:nBig),stimCols(2,:)+x(nBig+1:end),'filled','DisplayName','x')
% 
% figure, hold on
% plot(x0,'DisplayName','x0')
% plot(x,'DisplayName','x')
% legend('Location','best')

% figure, hold on
% scatter(data.trialdata.similarityMatrix(:),x,'MarkerFaceColor','k','MarkerFaceAlpha',0.05,'MarkerEdgeAlpha',0)
% xlim([0,1])
% ylim([0,1])
% axis square
% plot([0,1],[0,1])

%%

% [nll, simdata] = GenerativeModel(x,... 
%     'choiceInds',choiceInds,...
%     'cueInd',cueInd,...
%     'response',response,...
%     'nTrials',nTrials,...
%     'nBig',nBig, ...
%     'nSmall',nSmall,...
%     'optimisationMeta',optimisationMeta,...
% 	'forceNumerical',false);
% 
% figure,
% imagesc(simdata.trialdata.similarityMatrix)
% axis equal tight
% colormap('gray')
% colorbar
% caxis([0 1])
% xlabel('Choice')
% ylabel('Cue')

%%

% sm = simdata.trialdata.similarityMatrix;
% sm_cs = zeros(size(sm));
% for i = 1:nBig
%     sm_cs(i,:) = circshift(sm(i,:),(nBig-i)+(nBig/2));
% end
% 
% sm_cs = [sm_cs(:,end),sm_cs]; %wrap around
% 
% figure,
% imagesc(sm_cs)
% axis equal tight
% colormap('gray')
% colorbar
% caxis([0 1])
% xlabel('Choice (relative to cue)')
% ylabel('Cue')
% ticks = [0 nBig*.25 nBig*.5 nBig*.75 nBig]+1;
% ticklabels = {'-180', '-90', '0', '90', '180'};
% xticks(ticks)
% xticklabels(ticklabels) 


%%

% [~,simdata_recoveredValues] = f(x);
% 
% figure,
% imagesc(simdata_recoveredValues.trialdata.similarityMatrix)
% axis equal tight
% colormap('gray')
% colorbar
% caxis([0 1])
% xlabel('Choice')
% ylabel('Cue')
% 
% 
% [~,simdata_realValues] = f([1.827346,2.72365,20]);
% 
% figure,
% imagesc(simdata_realValues.trialdata.similarityMatrix)
% axis equal tight
% colormap('gray')
% colorbar
% caxis([0 1])
% xlabel('Choice')
% ylabel('Cue')

% %% Lagrange multipliers
% 
% figure,
% imagesc(reshape(x,[nBig, nBig]))
% axis equal tight
% colormap('gray')
% colorbar
% 
% figure,
% imagesc(reshape(LAMBDA.upper,[nBig, nBig]))
% axis equal tight
% colormap('gray')
% colorbar
% 
% figure,
% imagesc(reshape(LAMBDA.lower,[nBig, nBig]))
% axis equal tight
% colormap('gray')
% colorbar
