function [gMatrixBetween] = Generate_Between_STD(itpos, spikeMatrix, gMatrixWithin, pSTD, tau)



% itpos = timestamp data for behavior
% spikeMatrix = matrix of spike trains (one per neuron)
% gMatrixWithin = conductance matrix for within cell firing rates (one per neuron)

networkSize = size(gMatrixWithin,2);


%%%%%%%%%%%%%% COMPUTE CONDUCTANCES FOR BETWEEN CELL FIRING RATES

sep_ang(:,1) = 1:networkSize;                                           %separation angles to process
spike_B = zeros(length(gMatrixWithin), networkSize, length(sep_ang));
S_B = zeros(length(gMatrixWithin), networkSize, length(sep_ang));
gMatrixBetween = zeros(length(gMatrixWithin),length(sep_ang));          %conductance matrix for within cell firing rates (one per separation angle)
index_i(:,1) = 1:networkSize;


    for j = 1:length(sep_ang) %loop through separation angles

        index_j(:,1) = mod(index_i + sep_ang(j), networkSize);  %index for 2nd cell in the pair
        index_j(index_j == 0) = networkSize; 

        for i = 1:networkSize %loop through HD cells in the ring
stdspk=(spikeMatrix(:, index_i(i))-pSTD(:, index_i(i)));
stdspk(find(stdspk<0))=0;
            spike_B(:, index_i(i), j) = gMatrixWithin(:, index_j(i)).*stdspk;
            S_B(:, index_i(i), j) = Convolved_Conductance(itpos, 1, spike_B(:, index_i(i), j), tau);

        end
        gMatrixBetween(:,j) =  sum(S_B(:,:, j), 2);
    end
    
