function data = generateSimulatedData

nTrials = 50000;
nBig = 64;
nSmall = 4;

% % Create offsets (simulating an inhomogenous colorspace)

[stimCols,pol] = generateStimCols('nBig',nBig);
interval = 360/nBig;
magnitude = 20;
offset = cos(deg2rad(0:interval:360-interval)) * magnitude;

pol(1,:) = pol(1,:) + offset;
[cart(:,1),cart(:,2)] = pol2cart(deg2rad(pol(1,:)),pol(2,:));

stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); stimCols(1,:);stimCols(2,:)]);
colvals = im2double(stimCols_sRGB);

figure,
scatter(cart(:,1),cart(:,2),75,colvals,'filled')
axis equal
axis off
% saveas(gcf,fullfile(['stimulus-space-non-uniformity_', datestr(now,'yymmdd'), '.svg']))

figure, 
plot(diff(pol(1,:)),'k')
axis tight
xlabel('Cue')
ylabel('Interval (degrees)')

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
skewedGaussians = (sin(deg2rad(linspace(0,360,nBig)))/10)+0.5;
% skewedGaussians = skewedGaussians/(1/diff(r))+r(1);
%
% skewedGaussians = linspace(0.1,0.9,nBig); % linear skew
% skewedGaussians = ones(1,nBig)*0.5; % no skew

figure, plot(skewedGaussians); axis tight % arbitrary, just for testing

dPrime = 1;
gaussianWidth = 40;

disp(['dPrime = ',num2str(dPrime)])
disp(['gaussianWidth = ',num2str(gaussianWidth)])

[~, data] = GenerativeModel([],...
    'dPrime',dPrime,...
    'gaussianWidth', gaussianWidth,...
    'nTrials',nTrials,...
    'nBig',nBig,...
    'nSmall',nSmall,...
    'stimulusRemappingCart',[cart(:,1);cart(:,2)]' - [stimCols(1,:),stimCols(2,:)],... % TODO Switch to using 'stimulusRemappingPol', since that's what we're primarily using in the recovery/fitting. The effect is analagous here, but it would be easier to read
    'skewedGaussians',skewedGaussians);

data.trialdata.nBig     = nBig;
data.trialdata.nSmall   = nSmall;
data.trialdata.nTrials  = nTrials;

% choiceInds =    cell2mat(data.trialdata.choices)';
% cueInd =        cell2mat(data.trialdata.cues);
% response =      cell2mat(data.trialdata.chosen);
% stimCols =      data.trialdata.stimCols;

if exist('plotSimilarityMatrix','file') % Part of the MCC repo
    plotSimilarityMatrix(data.trialdata.similarityMatrix)
    title('x Similarity Matrix')
else
    figure,
    imagesc(data.trialdata.similarityMatrix') % transposed so cue on x-axis
    axis equal tight
    colormap('gray')
    colorbar
    caxis([0 1])
    xlabel('Choice')
    ylabel('Cue')
    title('x Similarity Matrix')
end

end
