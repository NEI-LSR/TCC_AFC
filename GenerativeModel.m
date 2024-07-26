function [nll, cleandata] = GenerativeModel(optimisationParams,varargin)

% Generates simulated data based on passed parameters
% ... using the TCC model [ref 1]

% [ref 1] Schurgin, M.W., Wixted, J.T. and Brady, T.F., 2020. Psychophysical scaling reveals a unified theory of visual memory strength. Nature Human Behaviour, 4(11), pp.1156–1172. https://doi.org/10.1038/s41562-020-00938-0.
% [ref 2] McGraw, K.O. and Wong, S.P., 1992. A common language effect size statistic. Psychological Bulletin, 111(2), p.361. https://doi.org/10.1037/0033-2909.111.2.361.

% Math I don't fully understand but which is relevent
% - https://stats.stackexchange.com/questions/50501/probability-of-one-random-variable-being-greater-than-another
% - https://stats.stackexchange.com/questions/243293/x-is-greater-than-y-z-x-n-when-all-normal-distributions?rq=1

%% Pre-flight
% The things that we need to set up before anything else happens

%rng(42); % fixes random number generator for reproducibility

default_similarityMatrix      = NaN;  % generate below if not passed
default_dPrime                = 1;    % Eyeballed from [ref 1], fig 3
default_stimulusRemappingCart = NaN;
default_stimulusRemappingPol  = NaN;
default_nBig                  = 64;   % the large n-AFC that we are trying to emulate
default_nSmall                = 4;    % the actual n we are using
default_nTrials               = 10000;% number of trials
default_CorrectOpt            = 1;    % always having a correct option presented as a response option
default_stimCols              = NaN;  % Specify stimulus. Default is set to NaN (Not a Number) so that it calls an additional function to generate stimuli.
default_choiceInds            = NaN;  % Specify response option indices. Default is set to NaN (Not a Number) so that it calls an additional function to generate stimuli.
default_cueInd                = NaN;  % Hard-coded stimulus. Default is set to NaN (Not a Number) so that it will normally generate stimuli options internally.
default_response              = NaN;  % If you're generating data, this function will generate chosen values, if you're not, it will return probabilities instead
default_attractorPoints       = NaN;  % NaN = no attractor points
default_attractorWeights      = NaN;  % NaN = no attractor weights
default_optimisationMeta      = NaN;  % optimisation metadata (What is contained in `optimisationParams`?)
default_forceNumerical        = false;
default_skewedGaussians       = NaN;
default_gaussianWidth         = 60;
default_pltSimFigs            = false; % plot similarity function graphs
default_offsetGaussians       = false;

errorMsg = 'Value must be positive, scalar, and numeric.';
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) && (x > 0),errorMsg);

ip = inputParser;
addParameter(ip,'similarityMatrix',default_similarityMatrix);
addParameter(ip,'dPrime',default_dPrime);
addParameter(ip,'stimulusRemappingCart',default_stimulusRemappingCart);
addParameter(ip,'stimulusRemappingPol',default_stimulusRemappingPol);
addParameter(ip,'nBig',default_nBig, validationFcn);
addParameter(ip,'nSmall',default_nSmall, validationFcn);
addParameter(ip,'nTrials',default_nTrials, validationFcn);
addParameter(ip,'CorrectOpt',default_CorrectOpt); %needs to be either 0 or 1, haven't coded in any validation check (yet).
addParameter(ip,'stimCols',default_stimCols); %Needs to be a list of numbers or NaN, haven't coded in any validation check (yet).
addParameter(ip,'choiceInds',default_choiceInds); %Needs to be a list of numbers or NaN, haven't coded in any validation check (yet).
addParameter(ip,'cueInd',default_cueInd);
addParameter(ip,'response',default_response);
addParameter(ip,'attractorPoints',default_attractorPoints)
addParameter(ip,'attractorWeights',default_attractorWeights)
addParameter(ip,'optimisationMeta',default_optimisationMeta)
addParameter(ip,'forceNumerical',default_forceNumerical);
addParameter(ip,'skewedGaussians',default_skewedGaussians)
addParameter(ip,'gaussianWidth',default_gaussianWidth)
addParameter(ip,'pltSimFigs',default_pltSimFigs)
addParameter(ip,'offsetGaussians',default_offsetGaussians)

parse(ip,varargin{:});

similarityMatrix      = ip.Results.similarityMatrix;
dPrime                = ip.Results.dPrime;
stimulusRemappingCart = ip.Results.stimulusRemappingCart;
stimulusRemappingPol  = ip.Results.stimulusRemappingPol;
nBig                  = ip.Results.nBig;
nSmall                = ip.Results.nSmall;
nTrials               = ip.Results.nTrials;
CorrectOpt            = ip.Results.CorrectOpt;
stimCols              = ip.Results.stimCols;
choiceInds            = ip.Results.choiceInds;
cueInd                = ip.Results.cueInd;
response              = ip.Results.response;
attractorPoints       = ip.Results.attractorPoints;
attractorWeights      = ip.Results.attractorWeights;
om                    = ip.Results.optimisationMeta;
forceNumerical        = ip.Results.forceNumerical;
skewedGaussians       = ip.Results.skewedGaussians;
gaussianWidth         = ip.Results.gaussianWidth;
pltSimFigs            = ip.Results.pltSimFigs;
offsetGaussians       = ip.Results.offsetGaussians;

% Overwrite if passed within optimisationMeta (used for model fitting)
if ~isnan(om)
    if om(1,1)
        similarityVector       = optimisationParams(1 : om(1,2));
        similarityMatrix       = reshape(similarityVector,[nBig, nBig]);
    end
    if om(2,1)
        dPrime                 = optimisationParams(sum(prod(om(1:1,:),2)) + 1 : sum(prod(om(1:2,:),2)));
    end
    if om(3,1)
        stimulusRemappingCart  = optimisationParams(sum(prod(om(1:2,:),2)) + 1 : sum(prod(om(1:3,:),2)));
    end
    if om(4,1)
        stimulusRemappingPol   = optimisationParams(sum(prod(om(1:3,:),2)) + 1 : sum(prod(om(1:4,:),2)));
    end
    if om(5,1)
        gaussianWidth          = optimisationParams(sum(prod(om(1:4,:),2)) + 1 : sum(prod(om(1:5,:),2)));
    end
    if om(6,1)
        skewedGaussians        = optimisationParams(sum(prod(om(1:5,:),2)) + 1 : sum(prod(om(1:6,:),2)));
    end
    if om(7,1)
        offsetGaussians       = optimisationParams(sum(prod(om(1:6,:),2)) + 1 : sum(prod(om(1:7,:),2)));
    end
end

if exist('stimulusRemappingPol','var') && all(~isnan(stimulusRemappingPol))
    if length(stimulusRemappingPol) ~= nBig
        currentInterval = 360/length(stimulusRemappingPol);
        % figure, hold on, plot(0:currentInterval:360,stimulusRemappingPol([1:end,1]),'*-');
        goalInterval = 360/nBig;
        stimulusRemappingPol = interp1(...
            0:currentInterval:360,...
            stimulusRemappingPol([1:end,1]),... % for loop around
            0:goalInterval:360);
        stimulusRemappingPol = stimulusRemappingPol(1:end-1);
        % plot(0:goalInterval:360-goalInterval,stimulusRemappingPol,'*-')
    end
end

if exist('skewedGaussians','var') && all(~isnan(skewedGaussians))
    if length(skewedGaussians) ~= nBig
        currentInterval = 360/length(skewedGaussians);
        % figure, hold on, plot(0:currentInterval:360,stimulusRemappingPol([1:end,1]),'*-');
        goalInterval = 360/nBig;
        skewedGaussians = interp1(...
            0:currentInterval:360,...
            skewedGaussians([1:end,1]),... % for loop around
            0:goalInterval:360);
        skewedGaussians = skewedGaussians(1:end-1);
        % plot(0:goalInterval:360-goalInterval,stimulusRemappingPol,'*-')
    end
end

if exist('offsetGaussians','var') && all(~isnan(offsetGaussians))
    if length(offsetGaussians) ~= nBig
        currentInterval = 360/length(offsetGaussians);
        % figure, hold on, plot(0:currentInterval:360,stimulusRemappingPol([1:end,1]),'*-');
        goalInterval = 360/nBig;
        offsetGaussians = interp1(...
            0:currentInterval:360,...
            offsetGaussians([1:end,1]),... % for loop around
            0:goalInterval:360);
        offsetGaussians = offsetGaussians(1:end-1);
        % plot(0:goalInterval:360-goalInterval,stimulusRemappingPol,'*-')
    end
end
%% Generate cues

if any(isnan(stimCols)) %is stimCols == NaN
    [stimCols, stimCols_pol] = generateStimCols('nBig',nBig);
end

if any(isnan(cueInd))
    cueInd = repmat(1:nBig,1,ceil(nTrials/nBig));
    cueInd = cueInd(1:nTrials); % chop down to size, in case nBig doesn't divide neatly into nTrials
    cueInd = cueInd(randperm(nTrials));
end

%cues = stimCols(:,cueInd);

%figure, hold on, axis equal
%scatter(cues(1,:),cues(2,:),'filled','DisplayName','Nominal cues','MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01); legend;

%% Stimulus remapping

% To account for a perceptually non-uniform colorspace

%disp(stimulusRemapping(1:10))

if ~isnan(stimulusRemappingCart)
    stimulusRemappingCart = reshape(stimulusRemappingCart, [nBig,2])';
    stimCols_pr = stimCols + stimulusRemappingCart; %stimcols post-remap
    [stimCols_pol_pr(1,:),stimCols_pol_pr(2,:)] = cart2pol(stimCols_pr(1,:),stimCols_pr(2,:));
    stimCols_pol_pr(1,:) = rad2deg(stimCols_pol_pr(1,:));
elseif ~isnan(stimulusRemappingPol)
    % normalise input to add to 360
    stimulusRemappingPol = stimulusRemappingPol*(360/sum(stimulusRemappingPol));
    % convert stimCols to polar
    [~,stimCols_pol_pr(2,:)] = cart2pol(stimCols(1,:),stimCols(2,:));
    stimCols_pol_pr(1,:) = cumsum([0,stimulusRemappingPol(1:end-1)]);
    % convert back to cartesian
    [stimCols_pr(1,:),stimCols_pr(2,:)] = pol2cart(deg2rad(stimCols_pol_pr(1,:)),stimCols_pol_pr(2,:));
else
    stimCols_pr = stimCols;
    stimCols_pol_pr = stimCols_pol;
end

%% Generate response options

if any(isnan(choiceInds)) % is respOpts == NaN
    allRespInd = 1:nBig; % all response options
    choiceInds = zeros(nSmall,nTrials); % response options indices
    if CorrectOpt % "correct" is always an option
        for i = 1:nTrials
            potDistsInd = allRespInd(allRespInd ~= cueInd(i));  % potential distractors: all the options apart from the cue
            choiceInds(1:nSmall, i) = [cueInd(i), potDistsInd(randperm(nBig-1, nSmall-1))];
        end
    else % if we're allowing the possibility that there might not always be the correct option displayed
        for i = 1:nTrials
            choiceInds(1:nSmall, i) = randperm(nBig,nSmall);
        end
    end
end

if size(choiceInds,2) ~= nTrials
    error('The size of the provided response options indices (respOptInds) does not equal the number of trials (nTrials)')
end

%respOpts = reshape(stimCols(:,choiceInds),[2,nSmall,nTrials]);
%scatter(respOpts(1,:),respOpts(2,:),'filled','DisplayName','Nominal response options');

%% Generate similarity matrix


if isnan(similarityMatrix)

    % Inspired by `TCC_Code_InManuscriptOrder\Figure1\F\ShapeOfSimFunction`
    % [ref 1] (https://osf.io/j2h65)

    D = zeros(nBig,nBig); % distance (angular)
    for i = 1:nBig
        D(i,:) = rad2deg(angdiff(deg2rad(stimCols_pol_pr(1,:)),deg2rad(stimCols_pol_pr(1,i))));
    end

    if pltSimFigs
        figure,
        imagesc(D)
        title('Distance')
        axis equal tight
        colormap('gray')
        colorbar
        %caxis([0 1])
        xlabel('Choice')
        ylabel('Cue')
        % saveas(gcf,['distanceMatrix_',datestr(now,'yymmdd'),'.svg'])
    end

    x = -180:0.1:180;

    if isnan(skewedGaussians) % no skew
        skewedGaussians = ones(nBig,1) * 0.5;
    end

    SplitGauss = @(x,sd_left,sd_right) [...
        exp(-((x(x<=0).^2)/(2*sd_left^2))),...
        exp(-((x(x> 0).^2)/(2*sd_right^2)))];

    simFunc = zeros(nBig,length(x));

    if isscalar(gaussianWidth)
        gaussianWidth = ones(nBig,1)*gaussianWidth;
    end
    
    for i = 1:nBig
        simFunc(i,:) = SplitGauss(x,...
            skewedGaussians(i)     * (2 * gaussianWidth(i)),...
            (1-skewedGaussians(i)) * (2 * gaussianWidth(i)));
    end

    if any(offsetGaussians)
        for i = 1:nBig
            simFunc(i,:) = circshift(simFunc(i,:),round(offsetGaussians(i)*abs(1/(x(1)-x(2)))),2); % TODO: Current limitation, actually only shifts to the nearest 0.1 degrees
        end
    end

    if pltSimFigs

        figure,
        if all(all(simFunc(1,:) == simFunc(:,:))) % if all the rows are the same
            plot(x,simFunc(1,:),'k', 'LineWidth',3)
            title('Similarity Function')
        else
            plot(x,simFunc(32,:),'k', 'LineWidth',3)
            title('Similarity function for cue 32')
        end
        axis tight
        xline(0,'k:')
        xlim([-180,180]);
        xticks(-180:45:180);
        ylim([0,1]);
        yticks([0,1]);
        box off
        xlabel('Degrees')
        ylabel('Similarity')
        % saveas(gcf,['SimilarityFunction_',datestr(now,'yymmdd'),'.svg'])

        figure, hold on, axis tight
        surf(simFunc,simFunc,'edgecolor','none')                            % TODO Make it so that the x-axis shows x
        title('Similarity Function Set')
        colormap('gray')
        xline(find(x==0),'--r')
        view(2)
        colorbar
        set(gca, 'YDir','reverse')
        % saveas(gcf,['SimilarityFunctionSet_',datestr(now,'yymmdd'),'.svg'])

    end

    similarityMatrix = NaN(nBig,nBig);
    for row = 1:nBig
        similarityMatrix(row,:) = simFunc(row,round(D(row,:)*abs(1/(x(1)-x(2)))) + find(x==0));
    end

    if pltSimFigs
        figure,
        imagesc(similarityMatrix)
        title('Similarity Matrix')
        axis equal tight
        colormap('gray')
        colorbar
        %caxis([0 1])
        xlabel('Choice')
        ylabel('Cue')
        % saveas(gcf,['SimilarityMatrix_',datestr(now,'yymmdd'),'.svg'])

        figure, hold on
        plot(1:nBig,similarityMatrix(32,:),'Color',[120,202,188]/255, 'LineWidth',2,'LineStyle','--')
        xline(32,'k:','LineWidth',2)
        % title('Similarity function for cue 32')
        axis tight
        xlim([1,nBig]);
        xticks(32);
        ylim([0,1]);
        yticks([0,1]);
        box off
        xlabel('Stimulus Index')
        ylabel('Similarity')
        % saveas(gcf,['SimilarityFunction_',datestr(now,'yymmdd'),'.svg'])
    end

end

if any(~isnan(attractorPoints(:))) && any(~isnan(attractorWeights(:)))

    if size(attractorPoints,1) == 2
        attractorPoints_u = attractorPoints(1,:);
        attractorPoints_v = attractorPoints(2,:);
    elseif size(attractorPoints,1) == 1
        attractorPoints_u = attractorPoints(1:length(attractorPoints)/2);
        attractorPoints_v = attractorPoints(length(attractorPoints)/2+1:end);
    end
    nAttactors = length(attractorPoints_u);

    % Create attractor similarity matrix

    sm_a = zeros(nBig,nBig,nAttactors); % similiarity matrix - attractors
    for i = 1:nAttactors

        % if stim A is within attractor,
        % stim A's similarity to *all the other stimuli*
        % that are also within radius i,
        % and that are closer to the stim A than the attractor point is
        % get max value (1)

        D_a(i,:) = sqrt(((attractorPoints_u(i) - stimCols_pr(1,:)).^2) + ((attractorPoints_v(i) - stimCols_pr(2,:)).^2)); % distances, attractor point
        isWithinAttractorRadius(i,:) = (stimCols_pr(1,:)-attractorPoints_u(i)).^2 + (stimCols_pr(2,:)-attractorPoints_v(i)).^2 < attractorWeights(i)^2; % Is each stimulus within the radius of the ith attractor point

        for stim = 1:nBig
            if isWithinAttractorRadius(i,stim)
                isOtherStimBetweenStimAndAttractor(i,stim,:) = (stimCols_pr(1,:)-attractorPoints_u(i)).^2 + (stimCols_pr(2,:)-attractorPoints_v(i)).^2 ...
                    < (stimCols_pr(1,stim)-attractorPoints_u(i)).^2 + (stimCols_pr(2,stim)-attractorPoints_v(i)).^2;
                for otherStim = 1:nBig
                    if isOtherStimBetweenStimAndAttractor(i,stim,otherStim)
                        distanceBetweenStimAndOtherStim = D(stim,otherStim);
                        distanceBetweenStimAndAttractorPoint = D_a(i,stim);
                        if distanceBetweenStimAndOtherStim ...
                                <= distanceBetweenStimAndAttractorPoint + 0.01 % SCARY MAGIC NUMBER BE HEREEE
                            sm_a(stim,otherStim,i) = 1;
                        end
                    end
                end
            end
        end

        %         figure,
        %         imagesc(sm_a(:,:,i))
        %         colorbar
        %         axis equal tight
        %         colormap('gray')
    end

    % Combine with standard similarity matrix

    similarityMatrix = max(cat(3,similarityMatrix,sm_a),[],3);

    %         figure,
    %         imagesc(similarityMatrix)
    %         colorbar
    %         axis equal tight
    %         colormap('gray')

end


%% Compute probabilty of each choice being selected

if and(nSmall<=5,~forceNumerical)
    if size(cueInd,1) ~= 1
        cueInd = cueInd';
    end

    cf = zeros(nSmall,nTrials); % choice familiarities
    cueInd_rep = repmat(cueInd,nSmall,1);
    idx_Lin = sub2ind([nBig,nBig],cueInd_rep(:),choiceInds(:));
    cf(:) = similarityMatrix(idx_Lin);
    % H/T to https://stackoverflow.com/a/72857254/6464224

    cf = cf*dPrime;

    zmat = (- permute(shiftdim(cf,-1),[2,1,3]) + shiftdim(cf,-1))/sqrt(2);

    pmat = normcdf(zmat); % this is the probability that one choice will be more familiar than another

    if nSmall == 2
        pmat_star = pmat;
    elseif nSmall == 3 % correction values from [ref 2]. Ideally, I'd find a better way to do this, because this is only rough, and not great
        pmat_star = 0.88*pmat + 0.11;
    elseif nSmall == 4
        pmat_star = 0.90*pmat + 0.15;
    elseif nSmall == 5
        pmat_star = 0.93*pmat + 0.17;
    else
        warning('No pstar correction listed for this value of nSmall. Continuing with uncorrected values.')
        pmat_star = pmat;
    end

    p_star = squeeze(prod(pmat_star)*(1/pmat_star(1))); % second part to balance out the '0.5' (times pstar factor) factor that occurs when comparing to self

    p_star_norm = p_star./sum(p_star);

    responseInd =       zeros(1,nTrials);
    probOfSelected =    zeros(1,nTrials);

    if isnan(response)
        % Actually apply noise and pick a choice (only if no response is provided as input)
        for i = 1:nTrials
            responseInd(i) = randsample(nSmall,1,true,p_star_norm(:,i));
            response(i) = choiceInds(responseInd(i),i);
            probOfSelected(i) = p_star_norm(responseInd(i),i);
        end
    else % if the response is passed, compute the probability of that response having been selected
        try
            responseInd = find(repmat(response',nSmall,1) == choiceInds);
            probOfSelected = p_star_norm(responseInd);
        catch % Juuuuust in case the above method doesn't work - The above is faster, but: there is a single trial in Pollux's data where there are two correct matches and that would break the above method (I think - untested)
            warning('using slow method')
            for i = 1:nTrials
                responseInd(i) = find(response(i) == choiceInds(:,i),1); % the '1' at the end of this means that it will only pick out the first match.
                probOfSelected(i) = p_star_norm(responseInd(i),i);
            end
        end
    end

    % Compute negative log likelihood
    nll = sum(-(log(probOfSelected)));

    %fprintf('%f\n',nll)

else
    xvals = -2.8:0.4:8; % This covers the required range for dPrimes up to about 5. E.g: figure, hold on; plot(xvals,normcdf(xvals,0,1)); plot(xvals,normcdf(xvals,5,1))
    cf = zeros(nSmall,nTrials); % choice familiarities

    for i = 1:nTrials
        cf(:,i) = similarityMatrix(cueInd(i),choiceInds(:,i));
    end

    % How likely is an item of every possible psych. distance to generate
    % every possible memory strength signal:
    normals = normcdf(repmat(xvals,nSmall,1,nTrials),...
        permute(repmat(cf,1,1,size(xvals,2)),[1,3,2]).*dPrime,...
        1);

    %normals(normals < 1e-10) = 1e-10; % to avoid a divide-by-0 issue

    % figure, hold on
    % plot(xvals,normals(:,:,1)')

    % How likely is each psych. distance to generate the maximum memory
    % match signal? Distribution of the maximums is just prod(CDFs):
    maxDist = prod(normals);

    % plot(xvals,maxDist(1,:,1),'k')
    probs = zeros(nSmall,nTrials);

    for i = 1:nTrials
        pdf = zeros(nSmall,1);
        for k = 1:nSmall
            % Remove the current signal from the max distribution(we want to know
            % the chance that ANOTHER color was larger than this one) and turn
            % this distribution into a PDF instead of a cdf:
            maxDistCur = maxDist(:,:,i) ./ normals(k,:,i);
            maxPdf = [0 diff(maxDistCur)];
            maxPdf = maxPdf ./ sum(maxPdf);

            % Now ask the chance this value was greater than all of the other
            % ones:
            chanceOverThisVal = (1-normals(k,:,i));
            pdf(k) = sum(maxPdf.*chanceOverThisVal);
        end
        probs(:,i) = pdf./sum(pdf);
    end

    responseInd =       zeros(1,nTrials);
    probOfSelected =    zeros(1,nTrials);

    if isnan(response)
        % Actually apply noise and pick a choice (only if no response is provided as input)
        for i = 1:nTrials
            responseInd(i) = randsample(nSmall,1,true,probs(:,i));
            response(i) = choiceInds(responseInd(i),i);
            probOfSelected(i) = probs(responseInd(i),i);
        end
    else % if the response is passed, compute the probability of that response having been selected
        for i = 1:nTrials
            responseInd(i) = find(response(i) == choiceInds(:,i),1); % the '1' at the end of this means that it will only pick out the first match. There is a single trial in Pollux's data where there are two correct matches...
            probOfSelected(i) = probs(responseInd(i),i);
        end
    end

    % Compute negative log likelihood
    nll = sum(-(log(probOfSelected)));
end


%% Convert into cleandata format to match real data format

if nargout > 1 % only run this if it is requested (it's slowwww)
    cleandata.trialdata.choices             = num2cell(choiceInds',2);
    cleandata.trialdata.cues                = num2cell(cueInd');
    cleandata.trialdata.chosen              = num2cell(response);
    cleandata.trialdata.stimCols            = stimCols;
    cleandata.trialdata.similarityMatrix    = similarityMatrix;
    cleandata.trialdata.dPrime              = dPrime;
end

end
