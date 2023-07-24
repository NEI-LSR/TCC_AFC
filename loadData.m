function loadData


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
