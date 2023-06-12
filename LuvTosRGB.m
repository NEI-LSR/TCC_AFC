function sRGB = LuvTosRGB(Luv)

% Converts CIELUV values to sRGB.
%
% Uses mainly PTB functions with the exception of 'whitepoint' which is an
% internal MATLAB function.
%
% Assumes D65.
%
% Usage example:
% sRGB = LuvTosRGB([ones(1,64)*70;generateStimCols('nBig',64)])
% image(permute(sRGB,[1,3,2]))
%
% To see it throw an error when something is out of gamut, try:
% sRGB = LuvTosRGB([ones(1,64)*70;generateStimCols('nBig',64,'sat',300)]);
% image(permute(sRGB,[1,3,2]))

if size(Luv,1) ~= 3 % transpose if needed
    Luv = Luv';
    if size(Luv,1) ~= 3 % if transposition didn't work, throw error
        error('expected 3 values')
    end
end

XYZ_D65 = [0.95047;1.0;1.08883];
XYZ = LuvToXYZ(Luv,XYZ_D65);

if any(XYZ < 0,'all')
    warning('gamut boundary issue')
end

sRGBlin = XYZToSRGBPrimary(XYZ);
sRGB = uint8(SRGBGammaCorrect(sRGBlin,0)');

