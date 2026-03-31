function [ prunedSpikeMatrix ] = SpikePrune(tpos, networkSize, spikeMatrix, tauInh, inhRadius)

% initialize vector of inhibitory neurons
inhibitoryNetwork(1:length(tpos),1:networkSize) = 1;

% vector of initial time points for inhibitory decay function
t0_inh = zeros(1, networkSize) - 1000000;
spikePruneProb = rand(length(tpos),networkSize);

for i = 1:length(tpos)
    
    % Remove any spikes that do not fire due to fast inhibition from the
    % previous run through the loop
    if i ~= 1
        inhibitoryNetwork(i,:) = 1 - exp(-(tpos(i)-t0_inh)/tauInh);
        spikeMatrix(i,:) = spikeMatrix(i,:).*(inhibitoryNetwork(i,:) >= spikePruneProb(i,:));
    end
    
    % Get index of spikes in ith row spikeMatrix
    spikeIndex = find(spikeMatrix(i,:) == 1);
    
    % if there are spikes in the current ith row of spikeMatrix, then
    % inhibit the current cell that recieved a spike and adjacent cells per
    % the value given by inhRadius (inhibition radius), in the next time
    % step.
    if ~isempty(spikeIndex)
        inhibitoryNetwork(i,:) = 1 - exp(-(tpos(i)-t0_inh)/tauInh);

        tempSpikeIndex = [];
        for j = 1:length(spikeIndex)
            tempSpikeIndex = [tempSpikeIndex, ( (spikeIndex(j)-inhRadius):(spikeIndex(j)+ inhRadius) ) ];
        end
        % find indexes that exceed the network size and mod them so they
        % wrap back around.
        plusNetworkSize = find(tempSpikeIndex > networkSize);
        tempSpikeIndex(plusNetworkSize) = mod(tempSpikeIndex(plusNetworkSize), networkSize);
        
        % find indexes that are less than 1 and adjust there value so they
        % wrap back around.
        neg0 = find(tempSpikeIndex <= 0);
        tempSpikeIndex(neg0) = networkSize + tempSpikeIndex(neg0);
        
        % remove non-unque indexes from tempSpikeIndex
        tempSpikeIndex = unique(tempSpikeIndex);
        
        spikeIndex = tempSpikeIndex;
    end
    
    % t0_inh gets current time time stamp 
    t0_inh(spikeIndex) = tpos(i);
    inhibitoryNetwork(i,:) = 1 - exp(-(tpos(i)-t0_inh)/tauInh);
    
end
prunedSpikeMatrix = spikeMatrix;
%totSpikes = sum(sum(spikeMatrix(cutTime:end,:)));

end % end function

%{
  if i ~= 1
        spikePruneProb = rand(1,networkSize );
        spikeMatrix(i,:) = spikeMatrix(i,:).*(inhibitoryNetwork(i,:) >= spikePruneProb);
    end
%}
