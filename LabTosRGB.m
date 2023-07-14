function sRGB = LabTosRGB(Lab)

% See LuvTosRGB for documentation

if size(Lab,1) ~= 3 % transpose if needed
    Lab = Lab';
    if size(Lab,1) ~= 3 % if transposition didn't work, throw error
        error('expected 3 values')
    end
end

XYZ_D65 = [0.95047;1.0;1.08883];
XYZ = LabToXYZ(Lab,XYZ_D65);

if any(XYZ < 0,'all')
    warning('gamut boundary issue')
end

sRGBlin = XYZToSRGBPrimary(XYZ);
sRGB = uint8(SRGBGammaCorrect(sRGBlin,0)');

