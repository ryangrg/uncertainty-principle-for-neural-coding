function [ ThetaMatrix ] = PlaceMultiPhaseBurstGenerator(interHD, itpos, neuronNum, cyclesperlap, basefreq,peakrate)

    
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
    

    numphases=8;
    SpikeMatrix = zeros(size(itpos,1), neuronNum);
    SpikeMatrix(1,1)=1; %reference cell spikes on first timestep
    lastRefSpike=itpos(1)+[0:(numphases-1)]*baseinterval/numphases;
    lastSpike = ones(neuronNum, numphases)*-1000;
    timeSinceLastThis = zeros(neuronNum, numphases)-1000;
    
    %cycdex=[2 7 10 14 20 21 30];
    cycdex=[1 2 3];
    
    for i = 1:size(itpos, 1)
      for ph=1:numphases  
        %generate reference cell spike train
        timeSinceLastRef(ph) = abs(itpos(i) - lastRefSpike(ph));
        if timeSinceLastRef(ph)>=baseinterval;
            timeSinceLastRef(ph)=0;
            lastRefSpike(ph)=itpos(i);
            SpikeMatrix(i,1+(ph-1)*(neuronNum/numphases)) = 1;
        end
        %generate non-reference cell spike trains
        for j=2:4%2:numphases:neuronNum
            targetInterval = mod((ph-1)*baseinterval/numphases+(cyclesperlap(cycdex(j-1))*interHD(i)/(2*pi))*baseinterval,baseinterval);
            timeSinceLastRef(ph) = abs(itpos(i) - lastRefSpike(ph));
            timeSinceLastThis(j-1,ph) = abs(itpos(i) - lastSpike(j-1,ph));
            if (timeSinceLastRef > targetInterval) & (timeSinceLastThis(j-1,ph) >= timeSinceLastRef(ph));
                SpikeMatrix(i,j+(ph-1)*(neuronNum/numphases)) = 1;
                lastSpike(j-1,ph) = itpos(i);
            end
        end
      end
    end

    for i=1:32
        temp(:,i)=conv(SpikeMatrix(:,i),gausswin(125),'same');
    end
    [ ThetaMatrix ] = InhomogeneousPoisson(temp*peakrate, 1);

end

