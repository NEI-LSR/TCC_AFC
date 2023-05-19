    function stop = optimplotx_SimilarityMatrix(x,optimValues,state,varargin)
% OPTIMPLOTX Plot current point at each iteration.
%
%   STOP = OPTIMPLOTX(X,OPTIMVALUES,STATE) plots the current point, X, as a
%   bar plot of its elements at the current iteration.
%
%   Example:
%   Create an options structure that will use OPTIMPLOTX
%   as the plot function
%       options = optimset('PlotFcns',@optimplotx);
%
%   Pass the options into an optimization problem to view the plot
%       fminbnd(@sin,3,10,options)

%   Copyright 2006-2010 The MathWorks, Inc.

stop = false;

% nBig = sqrt(length(x)-1);
% similarityMatrix = reshape(x(1:end-1),[nBig, nBig]);
% 
% dprime = x(end);

nBig = sqrt(length(x));
similarityMatrix = reshape(x,[nBig, nBig]);

% [~, data] = GenerativeModel([],...
%     'attractorPoints',generateStimCols('nBig',64),...
%     'attractorWeights',x,...
%     'nBig',64);
% 
% similarityMatrix = data.trialdata.similarityMatrix;

% c = forcezerominandonemax(x);
% c(1) = c(1) + 0.0001;
% c(2) = 1-(c(2) + 0.0001);

switch state
    case 'iter'
        if optimValues.iteration == 0
            % The 'iter' case is  called during the zeroth iteration,
            % but it now has values that were empty during the 'init' case

            plotx = imagesc(similarityMatrix);
            set(gca,'YDir','reverse')
            axis equal tight
            colormap('gray')
            colorbar
            caxis([0 1]) 

%             title(['Current similarity matrix', newline, 'dprime:', num2str(dprime), newline, 'mean min/maxes:', num2str(c)]);
%             title(['Current similarity matrix', newline, 'mean min/maxes:', num2str(c)]);
            xlabel('Choice');
            ylabel('Cue');

            set(plotx,'Tag','optimplotx');
        else
            plotx = findobj(get(gca,'Children'),'Tag','optimplotx');
            set(plotx,'Cdata',similarityMatrix);

%             title(['Current similarity matrix', newline, 'dprime:', num2str(dprime), newline, 'mean min/maxes:', num2str(c)]);
%            title(['Current similarity matrix', newline, 'mean min/maxes:', num2str(c)]);
        end
        try
            saveas(gcf,[pwd,filesep,'plots',filesep,datestr(now,'yyyymmddTHHMMSS'), '.png'])
        catch
        end
end

