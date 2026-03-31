function [ CV0spikeMatrix ] = CV0spikeGenerator(interHD, itpos, neuronNum, probScale, sigma)

    preferedHDstep = 2*pi/(neuronNum);
    preferedHD(:,1) = 0:preferedHDstep:(2*pi-preferedHDstep);

    l1 = length(interHD);
    repPreferedHD = repmat(preferedHD, 1, l1)';
    repIndexInterHD = repmat(interHD, 1, neuronNum);

    tau = 1 / (sigma*sigma);
    K = circ_vmpdf(0, 0, tau )/probScale;

    firingProbability = circ_vmpdf(repIndexInterHD, repPreferedHD, tau);
    minFiringProbability = min(firingProbability);
    repMinFiringProbability = repmat(minFiringProbability', 1, l1)';

    firingRate = 0.001./((firingProbability - repMinFiringProbability)./K);
    

    CV0spikeMatrix = zeros(size(itpos,1), neuronNum);
    lastSpike = ones(neuronNum, 1)*0.001;
    timeSinceLast = zeros(neuronNum, 1);
    
    for i = 1:size(itpos, 1)
        timeSinceLast = abs(itpos(i) - lastSpike);
        CV0spikeMatrix(i,:) = timeSinceLast > firingRate(i,:)';
        lastSpike(CV0spikeMatrix(i,:) == 1) = itpos(i);
    end

end

