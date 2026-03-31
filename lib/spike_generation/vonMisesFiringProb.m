function [ vonMisesSpikeMatrix ] = vonMisesFiringProb(interHD, neuronNum, probScale, sigma, seed)

preferedHDstep = 2*pi/(neuronNum);
preferedHD(:,1) = 0:preferedHDstep:(2*pi-preferedHDstep);

l1 = length(interHD);
repPreferedHD = repmat(preferedHD, 1, l1)'; % Create matrix of preffered HDs to use in circ_vmpdf
repActualHD = repmat(interHD, 1, neuronNum); % Create matrix of actual HDs to use in circ_vmpdf

tau = 1 / (sigma*sigma);
K = circ_vmpdf(0, 0, tau )/probScale;

firingProbability = circ_vmpdf(repActualHD, repPreferedHD, tau);
minFiringProbability = min(firingProbability);
repMinFiringProbability = repmat(minFiringProbability, l1, 1);

firingProbability = (firingProbability - repMinFiringProbability)./K;

% Changing firingProbability to logicals
rng(seed);
randNum = rand(l1,neuronNum);
vonMisesSpikeMatrix = firingProbability >= randNum;

end

