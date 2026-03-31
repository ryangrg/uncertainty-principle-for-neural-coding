function [ ThetaMatrix ] = PlacePhaseBurstGenerator(interHD, itpos, neuronNum, cyclesperlap, basefreq,peakrate)

    
    baseinterval=1/basefreq;
    
%     preferedHDstep = 2*pi/(neuronNum);
%     preferedHD(:,1) = 0:preferedHDstep:(2*pi-preferedHDstep);
% 
%     l1 = length(interHD);
%     repPreferedHD = repmat(preferedHD, 1, l1)';
%     repIndexInterHD = repmat(interHD, 1, neuronNum);
% 
%     tau = 1 / (sigma*sigma);
%     K = circ_vmpdf(0, 0, tau )/probScale;
% 
%     firingProbability = circ_vmpdf(repIndexInterHD, repPreferedHD, tau);
%     minFiringProbability = min(firingProbability);
%     repMinFiringProbability = repmat(minFiringProbability', 1, l1)';
% 
%     firingRate = 0.001./((firingProbability - repMinFiringProbability)./K);
    

    SpikeMatrix = zeros(size(itpos,1), neuronNum);
    SpikeMatrix(1,1)=1; %reference cell spikes on first timestep
    lastRefSpike=itpos(1);
    lastSpike = ones(neuronNum, 1)*-1000;
    timeSinceLastThis = zeros(neuronNum, 1)-1000;
    
    for i = 1:size(itpos, 1)
        %generate reference cell spike train
        timeSinceLastRef = abs(itpos(i) - lastRefSpike);       
        if timeSinceLastRef>=baseinterval;
            timeSinceLastRef=0;
            lastRefSpike=itpos(i);
            SpikeMatrix(i,1) = 1;
        end
        %generate non-reference cell spike trains
        for j=2:neuronNum
            targetInterval = mod((j-1)*baseinterval/neuronNum + (cyclesperlap(j-1)*interHD(i)/(2*pi))*baseinterval,baseinterval);
                             %   |----- stagger phases ------| 
            timeSinceLastRef = abs(itpos(i) - lastRefSpike);
            timeSinceLastThis(j-1) = abs(itpos(i) - lastSpike(j-1));
            if (timeSinceLastRef > targetInterval) & (timeSinceLastThis(j-1) >= timeSinceLastRef);
                SpikeMatrix(i,j) = 1;
                lastSpike(j-1) = itpos(i);
            end
        end
    end

    for i=1:32
        temp(:,i)=conv(SpikeMatrix(:,i),gausswin(125),'same');
    end
    [ ThetaMatrix ] = InhomogeneousPoisson(temp*peakrate, 1);

end

