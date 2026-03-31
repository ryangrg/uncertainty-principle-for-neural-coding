
clear;
load('noThetaSpks.mat');  %load file with spike matrix and behavior data

cf_Theta_testerr=[]; cf_Theta_testshift= [];
cf_testerr=[]; cf_testshift= [];
cf_Comp_testerr=[]; cf_Comp_testshift= [];

ncells=size(spikeMatrix,2); 

 [ cosTuningCurve ] = cosPiTuningCurve(ncells, 1);
 cosTuningCurve = [cosTuningCurve((round(ncells/2)+1):end); cosTuningCurve(1:round(ncells/2))];
 [ cosTuningCurve2 ] = cos90TuningCurve(ncells, 1);
 cosTuningCurve2  = cosTuningCurve2(end:-1:1);

dt = itpos(2)-itpos(1);  %simulaton time step


for ps=15
    
    %clearvars -except ti hd tpos vscale tau* ps centerfreq cf cf_* ncells basefrequency peakrate tau cosTuningCurve* seeds t1 t2 dt sigma networkType test* train* Count_*;
    tau_within=.1; %decay time constant for firing rates
    tau_between=tau_within; %decay time constant for co-firing rates
    tau_std=0.5; %recovery time constant for vesicle depletion (STD)
    
    perc_std=0;%.8*ps.^(-1.3);%percentage of vesicles used up per presynaptic spike
    
    % balance clk an ccw turning in the training set for symmetric velocity weights
    trainstart=18421;
    trainend=265930; %+(splitpoint-trainstart);
    trainrange=trainstart:trainend;
        
    teststart=trainend+1;
    testend=length(ixpos)-3000;
    testrange=teststart:testend;

    %%% TRAIN THE FIRING AND CO-FIRING RATE DECODERS (BUT FIRING RATE DECODER WILL USE COSINE WEIGHTS!!!)
    Spikes = spikeMatrix(trainrange,:); %get spikes for training segment
    Gstd = Generate_STD(itpos(trainrange), ncells, Spikes, tau_std, perc_std); %short term synaptic depression
    stdspk=Spikes-Gstd; stdspk(find(stdspk<0))=0;
    withinG = Generate_Conductance(itpos(trainrange), ncells, Spikes, tau_within); %firing rate traces
    betweenG = Generate_Between_STD(itpos(trainrange), Spikes, withinG, Gstd, tau_between); %co-firing rate traces
    [BW_weights, BetWith_weights, trainerr, trainshift] = Universal_Train(ixpos(trainrange), ispd(trainrange), withinG, betweenG)    
    
    %%% TEST THE FIRING AND CO-FIRING RATE DECODERS (BUT FIRING RATE DECODER WILL USE COSINE WEIGHTS!!!)
    Spikes = spikeMatrix(testrange,:); %get spikes for testing segment
    Gstd = Generate_STD(itpos(testrange), ncells, Spikes, tau_std, perc_std); %short term synaptic depression
    stdspk=Spikes-Gstd; stdspk(find(stdspk<0))=0;  
    withinG = Generate_Conductance(itpos(testrange), ncells, Spikes, tau_within); %firing rate traces
    betweenG = Generate_Between_STD(itpos(testrange), Spikes, withinG, Gstd, tau_between); %co-firing rate traces
    %re-scale the cosine weights to match the scale of fitted weights:
    BetWith_weights(:,3)=mean(abs([max(BetWith_weights(:,3)) min(BetWith_weights(:,3))]))*(cosTuningCurve-.5)*2;
    BetWith_weights(:,4)=mean(abs([max(BetWith_weights(:,4)) min(BetWith_weights(:,4))]))*(cosTuningCurve2-.5)*2;
    [testerr, testshift] = Universal_Test(ixpos(testrange), ispd(testrange), withinG, betweenG, BW_weights, BetWith_weights);
    cf_testerr=[cf_testerr; testerr(:)'];
    cf_testshift= [cf_testshift; testshift(:)'];
    

drawnow;

end
