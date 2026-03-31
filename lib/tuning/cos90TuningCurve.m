function [ cosTuningCurve ] = cos90TuningCurve(neuronNum, probScale)

preferedHDstep = 2*pi/(neuronNum);
preferedHD(:,1) = 0:preferedHDstep:(2*pi-preferedHDstep);

A1 = 0.5*probScale;
cosTuningCurve = -A1*cos(preferedHD - pi/2)+A1;

end

