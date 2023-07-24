clear, 
clc, 
close all

addpath(genpath('C:\toolbox\bads'))
addpath(genpath('C:\Users\cege-user\Documents\MacaqueColorCategories\Analyses'))

%%

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

% % stimulus remapping (angle) 
% lb = zeros(1,nBig) + 0.01; % lower bound
% ub = (ones(1,nBig) *  100) + rand(1,nBig); % upper bound 
% x0 = ones(1,nBig) + rand(1,nBig);
% %

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

% stimulus remapping AND skewed gaussians
lb = zeros(1,nBig*2) + 0.001; % lower bound
ub = ones(1,nBig*2) + rand(1,nBig*2)/100; % upper bound 
x0 = ones(1,nBig*2)*0.5 + (rand(1,nBig*2)-0.5);
%

optimisationMeta = [...      % what do the x values represent in the fitting? (1st column - is it passed, 2nd column - how many values would it be if it were)
    false,   nBig*nBig;...   % free similarityMatrix
    false,   1;...           % dprime
    false,   nBig*2;...      % stimulus remapping (cartesian)
    true,   nBig;...        % stimulus remapping (polar, angle in degrees)
    false,   nAttactors*2;... % attractor points
    false,   nAttactors;...   % attractor weights
    false,    1;...           % SimFunc_sd
    true,    nBig,...        % skewedGaussians              
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
     'DiffMinChange',0.1);

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
