function [cart,pol] = generateStimCols(varargin)

% Generates stimulus values 

%% Run as script instead of function (good for bugfixing)

if ~exist('varargin','var')
    clc, clear, close all
    varargin = {}; %sets varargin to empty so that the defaults are used
end

%% Pre-flight
% The things that we need to set up before anything else happens

default_nBig        = 28;      % the large n-AFC that we are trying to emulate
default_origin      = [0,0];   % origin around which the stimuli are specified (might not be required/useful)
default_fullSatOnly = 1;       % Full saturation only? (a ring around the extent of colourspace or a grid across colourspace?)
default_sat         = 38;      % Saturation (38, taken from Bae+ 2015)
default_showFig     = false;

errorMsg = 'Value must be positive, scalar, and numeric.';
validationFcn = @(x) assert(isnumeric(x) && isscalar(x) && (x > 0),errorMsg);

p = inputParser;
addParameter(p,'nBig',default_nBig, validationFcn);
addParameter(p,'origin',default_origin); %needs to be two values, haven't coded in any validation check (yet).
addParameter(p,'fullSatOnly',default_fullSatOnly); %needs to be either 0 or 1, haven't coded in any validation check (yet).
addParameter(p,'sat',default_sat,validationFcn);
addParameter(p,'showFig',default_showFig);

parse(p,varargin{:});

nBig        = p.Results.nBig;
origin      = p.Results.origin;
fullSatOnly = p.Results.fullSatOnly;
sat         = p.Results.sat;
showFig     = p.Results.showFig;

%disp(char(['generateStimCols is using defaults for:' p.UsingDefaults])) % Show warning saying which defaults are being used. Would be nice if it could list what the values for these parameters are.

if ~isequal(origin, [0,0])
    error('The code for generating stimuli not centred on [0,0] has not yet been writtten')
end

%% Generate stimulus colours

if fullSatOnly == 1
    % could be done with 'nsidedpoly' instead, then pulling out Vertices.
    
    polarAngs = linspace(0,360-(360/nBig),nBig); %Polar Angles
    %figure, polarplot(deg2rad(polarAngs),ones(nBig,1)*sat)
    
    [a,b] = pol2cart(deg2rad(polarAngs),ones(1,nBig)*sat); %Convert polar co-ordinates to cartesian co-ordinates
        
    cart = [a;b];
    pol  = [polarAngs;repmat(sat,size(polarAngs))];
    
elseif strcmp(fullSatOnly,'tritess') %triangular tesselation
    if nBig == 19
        % This gives close to what we want, but not perfect.
        % https://www.mathworks.com/matlabcentral/answers/24614-circle-packed-with-circles
        %
        % Would be nice to let this happen automatically, for different
        % nBig values, but this seems like it would be fairly tricky so it
        % is hard coded and slightly hacky, at least for now.
        
        polarAngs = linspace(0,360-(360/6),6); %Compute 6 exterior points
        [a,b] = pol2cart(deg2rad(polarAngs),ones(1,6)*sat); %Convert polar co-ordinates to cartesian co-ordinates
        cart = [a;b];
        for i = 1:6
            inbetween = [cart(1,i)/2;cart(2,i)/2];
            cart = [cart,inbetween];
        end
        for i = 1:5
            inbetween = [(cart(1,i)+cart(1,i+1))/2;(cart(2,i)+cart(2,i+1))/2];
            cart = [cart,inbetween];
        end
        inbetween = [(cart(1,6)+cart(1,1))/2;(cart(2,6)+cart(2,1))/2];
        cart = [cart,inbetween];
        cart = [cart,[0;0]]; %add origin
        cart = cart(:,[1,13,2,14,3,15,4,16,5,17,6,18,7,8,9,10,11,12,19]); % I hate this.
        
        [pol(1,:),pol(2,:)] = cart2pol(cart(1,:),cart(2,:));
        pol(1,:) = round(rad2deg(pol(1,:)));
        pol(1,pol(1,:)<0) = 360+pol(1,pol(1,:)<0); % the rad2deg remapping puts 190deg as -170deg (for example) - I want everything positive from 0 -> 360.
        
    else
        error('The code for triangular tesselation for nBig other than 19 has not been written (yet?)')
    end
else
    error('Unrecognised value for fullSatOnly')
end

if showFig
    figure(1000), hold on, axis equal
    stimCols_sRGB = LuvTosRGB([repelem(76.0693, nBig); cart(1,:);cart(2,:)]);
    colvals = im2double(stimCols_sRGB);
    scatter(cart(1,:),cart(2,:),[],colvals,"filled")
    % un-comment if you want labels
%     for t = 1:nBig
%         text(cart(1,t)+1.5, cart(2,t), num2str(t));
%     end
end

end

