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
default_dprime                = 1;    % Eyeballed from [ref 1], fig 3
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
default_lambda                = -0.1; % similarity function
default_sigma                 = 7;    % perceptual function

errorMsg = 'Value must be positive, scalar, and numeric.';
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) && (x > 0),errorMsg);

ip = inputParser;
addParameter(ip,'similarityMatrix',default_similarityMatrix);
addParameter(ip,'dprime',default_dprime);
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
addParameter(ip,'lambda',default_lambda);
addParameter(ip,'sigma',default_sigma);

parse(ip,varargin{:});

similarityMatrix      = ip.Results.similarityMatrix;
dprime                = ip.Results.dprime;
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
lambda                = ip.Results.lambda;
sigma                 = ip.Results.sigma;

if ~isnan(om)
    if om(1,1)
        similarityVector       = optimisationParams(1 : om(1,2));
        similarityMatrix       = reshape(similarityVector,[nBig, nBig]);
    end
    if om(2,1)
        dprime                 = optimisationParams(sum(prod(om(1:1,:),2)) + 1 : sum(prod(om(1:2,:),2)));
    end
    if om(3,1)
        stimulusRemappingCart  = optimisationParams(sum(prod(om(1:2,:),2)) + 1 : sum(prod(om(1:3,:),2)));
    end
    if om(4,1)
        stimulusRemappingPol   = optimisationParams(sum(prod(om(1:3,:),2)) + 1 : sum(prod(om(1:4,:),2)));
    end
    if om(5,1)
        attractorPoints        = optimisationParams(sum(prod(om(1:4,:),2)) + 1 : sum(prod(om(1:5,:),2)));
    end
    if om(6,1)
        attractorWeights       = optimisationParams(sum(prod(om(1:5,:),2)) + 1 : sum(prod(om(1:6,:),2)));
    end
    if om(7,1)
        lambda                 = optimisationParams(sum(prod(om(1:6,:),2)) + 1 : sum(prod(om(1:7,:),2)));
    end
    if om(8,1)
        sigma                  = optimisationParams(sum(prod(om(1:7,:),2)) + 1 : sum(prod(om(1:8,:),2)));
    end
end


%% Generate cues

if any(isnan(stimCols)) %is stimCols == NaN
    [stimCols, ~] = generateStimCols('nBig',nBig);
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
    %disp(stimCols_pr(1,1:10))
elseif ~isnan(stimulusRemappingPol)
    % convert stimCols to polar
    [stimCols_pol(1,:),stimCols_pol(2,:)] = cart2pol(stimCols(1,:),stimCols(2,:));
    stimCols_pol(1,:) = rad2deg(stimCols_pol(1,:));
    stimCols_pol(1,stimCols_pol(1,:) < 0) = stimCols_pol(1,stimCols_pol(1,:) < 0) + 360;
    % apply stimulus remapping
    stimCols_pol(1,:) = stimCols_pol(1,:) + stimulusRemappingPol;
    stimCols_pol(stimCols_pol < 0)    = stimCols_pol(stimCols_pol < 0)    + 360;
    stimCols_pol(stimCols_pol >= 360) = stimCols_pol(stimCols_pol >= 360) - 360;
    % convert back to cartesian
    [stimCols_pr(1,:),stimCols_pr(2,:)] = pol2cart(deg2rad(stimCols_pol(1,:)),stimCols_pol(2,:));
else
    stimCols_pr = stimCols;
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
    %
    % There's probably a slicker way of doing this, but I like that it is
    % set up so that you can seperately control the similarity function and
    % the perceptual function seperately, should we choose/need to do that
    % at some point.

    x = 0:0.1:300; % be careful if changing this - it will modify the convolution below, and also the scaling used to generate the similarity matrix (which depends on the interval, and the fact that this range starts at zero)

    %lambda = -0.1;
    similarityFunc = exp(x * lambda);

    %sigma = 7; %
    perceptualFunc = normpdf(x,0,sigma);

    combinedSim = conv(...
        padarray(similarityFunc,[0,length(similarityFunc)],'replicate'),...
        perceptualFunc,...
        'valid');
    combinedSim = combinedSim(2:length(perceptualFunc)+1);
    combinedSim = combinedSim/max(combinedSim);


    D = zeros(nBig,nBig);
    for i = 1:nBig
        D(i,:) = sqrt(((stimCols_pr(1,i) - stimCols_pr(1,:)).^2) + ((stimCols_pr(2,i) - stimCols_pr(2,:)).^2));
    end

%     figure, hold on,
%     legend
%     axis square
%     plot(x,similarityFunc/max(similarityFunc),':', 'LineWidth', 3, 'DisplayName','Similarity Function')
%     %figure,
%     plot(x,perceptualFunc/max(perceptualFunc),':', 'LineWidth', 3, 'DisplayName','Perceptual Function')
%     %figure,
%     plot(x,combinedSim, '-', 'Color', [.1 .7 .1], 'LineWidth', 3, 'DisplayName', 'Combined Function')
%     set(gca, 'XLim', [0 max(D(:))], 'XTick', [0,max(D(:))/2,max(D(:))], 'XTickLabels', {0,90,180}, 'YTick', [0 0.5 1], 'FontSize', 16);
%     xlabel('Distance between colors')
%     ylabel('Psychophysical Similarity')

%     figure, hold on,
%     axis square
%     %plot([-(x(end:-1:1)),x],[combinedSim(end:-1:1),combinedSim], '-', 'Color', [.1 .7 .1], 'LineWidth', 3, 'DisplayName', 'Combined Function')
%     plot([-(x(end:-1:1)),x],[combinedSim(end:-1:1),combinedSim]*3, '-', 'Color', [.1 .7 .1], 'LineWidth', 3, 'DisplayName', 'Combined Function')
%     set(gca, 'XLim', [-max(D(:)) max(D(:))], 'XTick', [-max(D(:)),0,max(D(:))], 'XTickLabels', {-180,0,180}, 'Ylim', [-5,5], 'YTick', [-5 -3 -1 1 3 5], 'FontSize', 16, 'YAxisLocation','right');
%     xlabel('Distance between colors')
%     ylabel('Familiarity')
% 
%     plot([-(x(end:-1:1)),x],[combinedSim(end:-1:1),combinedSim]*3+randn(1,size(combinedSim,2)*2), 'ok', 'DisplayName', 'Combined Function')

%     figure,
%     imagesc(D)
%     axis equal tight
%     colormap('gray')
%     colorbar
%     %caxis([0 1])
%     xlabel('Choice')
%     ylabel('Cue')

    similarityMatrix = combinedSim(round(D*abs(1/(x(1)-x(2))))+1);

%     figure,
%     imagesc(similarityMatrix)
%     axis equal tight
%     colormap('gray')
%     colorbar
%     caxis([0 1])
%     xlabel('Choice')
%     ylabel('Cue')

%     figure, plot(similarityMatrix(:,1))

end

if any(~isnan(attractorPoints)) && any(~isnan(attractorWeights))

    attractorPoints_u = attractorPoints(1:length(attractorPoints)/2);
    attractorPoints_v = attractorPoints(length(attractorPoints)/2+1:end);

       % Create attractor similarity matrix

    sm_a = zeros(nBig,nBig,length(attractorPoints)/2); % similiarity matrix - attractors
    for i = 1:length(attractorPoints)/2

        % if stim A is within attractor, 
        % stim A's similarity to *all the other stimuli* 
        % that are also within radius i, 
        % and that are closer to the stim A than the attractor point is
        % get max value (1)

        D_a(i,:) = sqrt(((attractorPoints_u(i) - stimCols_pr(1,:)).^2) + ((attractorPoints_v(i) - stimCols_pr(2,:)).^2)); % distances, attractor point
        isWithinAttractorRadius(i,:) = (stimCols_pr(1,:)-attractorPoints_u(i)).^2 + (stimCols_pr(2,:)-attractorPoints_v(i)).^2 < attractorWeights(i)^2; % Is each stimulus within the radius of the ith attractor point

        for stim = 1:nBig
            if isWithinAttractorRadius(i,stim)
                isOtherStimBetweenStimAndAttractor(i,stim,:) = (stimCols_pr(1,:)-attractorPoints_u(i)).^2 + (stimCols_pr(2,:)-attractorPoints_v(i)).^2 < (stimCols_pr(1,stim)-attractorPoints_u(i)).^2 + (stimCols_pr(2,stim)-attractorPoints_v(i)).^2;
                for otherStim = 1:nBig
                    if isOtherStimBetweenStimAndAttractor(i,stim,otherStim)
                        distanceBetweenStimAndOtherStim = D(stim,otherStim);
                        distanceBetweenStimAndAttractorPoint = D_a(i,stim);
                        if distanceBetweenStimAndOtherStim <= distanceBetweenStimAndAttractorPoint + 0.01 % SCARY MAGIC NUMBER BE HEREEE
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

    cf = cf*dprime;

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
    xvals = -2.8:0.4:8; % This covers the required range for dprimes up to about 5. E.g: figure, hold on; plot(xvals,normcdf(xvals,0,1)); plot(xvals,normcdf(xvals,5,1))
    cf = zeros(nSmall,nTrials); % choice familiarities

    for i = 1:nTrials
        cf(:,i) = similarityMatrix(cueInd(i),choiceInds(:,i));
    end

    % How likely is an item of every possible psych. distance to generate
    % every possible memory strength signal:
    normals = normcdf(repmat(xvals,nSmall,1,nTrials),...
        permute(repmat(cf,1,1,size(xvals,2)),[1,3,2]).*dprime,... 
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
    cleandata.trialdata.chosen              = num2cell(response');
    cleandata.trialdata.stimCols            = stimCols;
    cleandata.trialdata.similarityMatrix    = similarityMatrix;
    cleandata.trialdata.dprime              = dprime;
end

end
