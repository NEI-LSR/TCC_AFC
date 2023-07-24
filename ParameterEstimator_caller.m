clear, clc, close all

data = generateSimulatedData;

%%

params.dPrime               = true;
params.gaussianWidth        = true;

rn = 0;

[x,aic,bic] = ParameterEstimator(data,params,rn);

disp(x)