
clear;
load('HeadDirectionData.mat');

cf_Theta_testerr=[]; cf_Theta_testshift= [];
cf_CV1_testerr=[]; cf_CV1_testshift= [];
cf_CV0_testerr=[]; cf_CV0_testshift= [];
cf_Comp_testerr=[]; cf_Comp_testshift= [];

figure(300); clf;
figure(301); clf;
figure(302); clf;
ncells=32; basefrequency=7; peakrate=.15;
tau=.1;

[ cosTuningCurve ] = cosPiTuningCurve(32, 1);
cosTuningCurve = [cosTuningCurve((round(ncells/2)+1):end); cosTuningCurve(1:round(ncells/2))];
[ cosTuningCurve2 ] = cos90TuningCurve(32, 1);
cosTuningCurve2  = cosTuningCurve2(end:-1:1);

load('SeedSet_10x120.mat');

sigma = 0.5; %standard deviation of vin Mises HD tuning curve

t1 = 30;   %simulation start time (s)
t2 = 930;   %simulation end time (s)
dt = .001;  %simulaton time step

networkType = 'UC'; %uncompressed (for file naming)

%tauInh = [1.3 0.565 0.305 0.2 0.138 0.099 0.077 0.0575 0.046 0.038 0.0325 0.0275 0.023 0.02 0.01655 0.015];

pmod=[2.6 1.27 .825 .66 .625 .66 .74 .84];
tauInh = [.004 .0172 .044 .07 .0812 .080 .072 .0625];
%for centerfreq=1:15
for ps=1
%for vscale=1:4
    
 %   ps=8;
 sigma=2*.2;
%    peakrate = .005+(ps-1)*.0025;
    peakrate = .05*pmod(2);
    ti=tauInh(2);
    
    cf=.1;%*centerfreq;
    
    clearvars -except hd_real ti vscale tauInh ps pmod centerfreq cf cf_* ncells basefrequency peakrate tau cosTuningCurve* seeds t1 t2 dt sigma networkType test* train* Count_*;
    
    %duration=9000
    %trunc=2500;
    duration=15000
    trunc=5000;
    tpos=(1/30):(1/30):duration;
    av=rand(30*duration,1);
%    dfh = designfilt('bandpassiir', 'StopbandFrequency1', cf*.1, 'PassbandFrequency1', cf*.2, 'PassbandFrequency2', cf*1.8, 'StopbandFrequency2', cf*1.9, 'StopbandAttenuation1', 30, 'PassbandRipple', 1, 'StopbandAttenuation2', 30, 'SampleRate', 30);
    dfh = designfilt('bandpassiir', 'StopbandFrequency1', .01, 'PassbandFrequency1', .02, 'PassbandFrequency2', .2*1.8, 'StopbandFrequency2', .2*1.9, 'StopbandAttenuation1', 30, 'PassbandRipple', 1, 'StopbandAttenuation2', 30, 'SampleRate', 30);
    bpav=filtfilt(dfh,av*.66);
    bpav=bpav-.25*min(bpav(trunc*30:(duration-trunc)*30));
    bpav(find(bpav<.25*max(bpav(trunc*30:(duration-trunc)*30))))=0;
    gw=90;
    bpav=3*conv(bpav,gausswin(gw)/sum(gausswin(gw)),'same');
    figure(1); clf;
    subplot(2,1,1);
    plot(tpos((trunc*30:(duration-trunc)*30)),360*av((trunc*30:(duration-trunc)*30))); 
    hold on; plot(tpos((trunc*30:(duration-trunc)*30)),360*bpav((trunc*30:(duration-trunc)*30))/(2*pi/30));
    subplot(2,1,2);
    hd=mod(cumsum(bpav(trunc*30:(duration-trunc)*30)),2*pi);
    plot(tpos((trunc*30:(duration-trunc)*30)),hd);
    tpos=tpos((trunc*30:(duration-trunc)*30));
    tpos=tpos-tpos(1);
    %figure(2); clf; scatter(bpav(trunc*30:(duration-trunc)*30),hd);
    
    spd=bpav(trunc*30:(duration-trunc)*30);
    figure(1); clf; histogram(spd*30*40/(2*pi),0:2:50);
    drawnow;
    tic;
    %networkSize = 128; %number of neurons in the HD ring
    
    % Upsampled head direction data to compare with predicted head direction
    %
    ihd(:,1) = interpolateHDData(hd, tpos, dt, t1, t2);
    % balance clk an ccw turning in the training set for symmetric velocity weights
    trainstart=4569;
    splitpoint=313513;
    trainend=splitpoint+(splitpoint-trainstart);
    trainrange=trainstart:trainend;
        
    %ihd(splitpoint:trainend)=2*pi-ihd(trainstart:splitpoint);
    %ihd(647564)=mean([ihd(647563) ihd(647565)]); %AV(636744)=AV(636743);
    %trainend=13334+66667*9;
    teststart=trainend+1;
    testend=teststart+66667*3;
    testrange=teststart:testend;
    
    % interpolate tpos using interp1 takes orginal signal and interpolates into steps
    % determined by dt.
    %
    itpos(:,1) = interp1(tpos(1:end), tpos(1:end), 0:dt:((t2 - t1) - dt))';
    
    % variable to cut out conductance warmup from predictedHD and meanSpikes,
    % time will be cut from the start of itpos. Enter varialbe in seconds
    % desired to be cut.
    cutTime = 1;
    cutTime = cutTime*1000+1;
    
    % Compute angular velocity as circular difference of HD
    AV = circ_dist(ihd(2:end), ihd(1:end-1)); 
    %AV(splitpoint-1)=mean([AV(splitpoint-2) AV(splitpoint)]);
    %AV(trainend)=mean([AV(trainend-1) AV(trainend+1)]);
%    AV(325039:(325039+311705))=-AV(325039:-1:trainstart);
    
    
    numseeds = 1;
    probScale(1, :) = [25:2.5:150];
    
    
    seedSpikes = zeros(numseeds, 1);        %total number of spikes, M, for a given seed
    seedMSE = zeros(numseeds, 1);           %mean squared error for predicted AV
    M_Spikes = zeros(numseeds, length(probScale));
    M_MSE = zeros(numseeds, length(probScale));
    
    %cyclesperlap, basefreq,peakrate,tau)
    
    %spikeMatrix = PlacePhaseSpikeGenerator(behdata, itpos, networkSize, cyclesperlap, basefreq, peakrate); %create spike trains
    %spikeMatrix = PlacePhaseBurstGenerator(behdata, itpos, networkSize, cyclesperlap, basefreq, peakrate); %create spike trains
    
%     seed=1;
%     CV0_Spikes = CV0spikeGenerator(ihd(trainrange), itpos(trainrange), ncells, peakrate, sigma); %create spike trains
%     CV0_G = Generate_Conductance(itpos(trainrange), ncells, CV0_Spikes,tau); %convolve with exponential decay kernel
%     [CV0_BW_weights, CV0_BetWith_weights, CV0_trainerr, CV0_trainshift] = Universal_Train(1, itpos(trainrange), ihd(trainrange), AV(trainrange), CV0_Spikes, CV0_G);
%     
%     CV0_Spikes = CV0spikeGenerator(ihd(testrange), itpos(testrange), ncells, peakrate, sigma); %create spike trains
%     Count_CV0(ps)=sum(CV0_Spikes(:));
%     CV0_G = Generate_Conductance(itpos(testrange), ncells, CV0_Spikes,tau); %convolve with exponential decay kernel
%     %scale pure cosine weights to match scale of fitted weights
%     CV0_BetWith_weights(:,3)=mean(abs([max(CV0_BetWith_weights(:,3)) min(CV0_BetWith_weights(:,3))]))*(cosTuningCurve-.5)*2;
%     CV0_BetWith_weights(:,4)=mean(abs([max(CV0_BetWith_weights(:,4)) min(CV0_BetWith_weights(:,4))]))*(cosTuningCurve2-.5)*2;
%     [CV0_testerr, CV0_testshift] = Universal_Test(1, itpos(testrange), ihd(testrange), AV(testrange), CV0_Spikes, CV0_G, CV0_BW_weights, CV0_BetWith_weights);
%     cf_CV0_testerr=[cf_CV0_testerr; CV0_testerr(:)'];
%     cf_CV0_testshift= [cf_CV0_testshift; CV0_testshift(:)'];
%     
% %    figure(300); subplot(4,4,centerfreq); plot(-31:2:31,CV0_BetWith_weights(:,5)); axis tight;
%     figure(300); subplot(4,4,ps); plot(-31:2:31,CV0_BetWith_weights(:,5)); axis tight;
%     drawnow;
%     
%         seed=1;
%         CV1_Spikes = vonMisesFiringProb(ihd(trainrange), ncells, peakrate, sigma, seed); %create spike trains
%         CV1_G = Generate_Conductance(itpos(trainrange), ncells, CV1_Spikes,tau); %convolve with exponential decay kernel
%         [CV1_BW_weights, CV1_BetWith_weights, CV1_trainerr, CV1_trainshift] = Universal_Train(1, itpos(trainrange), ihd(trainrange), AV(trainrange), CV1_Spikes, CV1_G);
%     
%         CV1_Spikes = vonMisesFiringProb(ihd(testrange), ncells, peakrate, sigma, seed); %create spike tests
%         Count_CV1(ps)=sum(CV1_Spikes(:));
%         CV1_G = Generate_Conductance(itpos(testrange), ncells, CV1_Spikes,tau); %convolve with exponential decay kernel
%         %scale pure cosine weights to match scale of fitted weights
%         CV1_BetWith_weights(:,3)=mean(abs([max(CV1_BetWith_weights(:,3)) min(CV1_BetWith_weights(:,3))]))*(cosTuningCurve-.5)*2;
%         CV1_BetWith_weights(:,4)=mean(abs([max(CV1_BetWith_weights(:,4)) min(CV1_BetWith_weights(:,4))]))*(cosTuningCurve2-.5)*2;
%         [CV1_testerr, CV1_testshift] = Universal_Test(1, itpos(testrange), ihd(testrange), AV(testrange), CV1_Spikes, CV1_G, CV1_BW_weights, CV1_BetWith_weights);
%         cf_CV1_testerr=[cf_CV1_testerr; CV1_testerr(:)'];
%         cf_CV1_testshift= [cf_CV1_testshift; CV1_testshift(:)'];
%         
%     figure(301); subplot(4,4,ps); plot(-31:2:31,CV1_BetWith_weights(:,5)); axis tight;
%     drawnow;
% 
%         seed=1;
%         CV1_Spikes = vonMisesFiringProb(ihd(trainrange), ncells, .26, sigma, seed); %create spike trains
%         Comp_Spikes = SpikePrune(itpos(trainrange), ncells, CV1_Spikes, ti, 2);
%         Comp_G = Generate_Conductance(itpos(trainrange), ncells, Comp_Spikes,tau); %convolve with exponential decay kernel
%         [Comp_BW_weights, Comp_BetWith_weights, Comp_trainerr, Comp_trainshift] = Universal_Train(1, itpos(trainrange), ihd(trainrange), AV(trainrange), Comp_Spikes, Comp_G);
%     
%         CV1_Spikes = vonMisesFiringProb(ihd(testrange), ncells, .26, sigma, seed); %create spike tests
%         Comp_Spikes = SpikePrune(itpos(testrange), ncells, CV1_Spikes, ti, 2);
%         Count_Comp(ps)=sum(Comp_Spikes(:));
%         Comp_G = Generate_Conductance(itpos(testrange), ncells, Comp_Spikes,tau); %convolve with exponential decay kernel
%         %scale pure cosine weights to match scale of fitted weights
%         Comp_BetWith_weights(:,3)=mean(abs([max(Comp_BetWith_weights(:,3)) min(Comp_BetWith_weights(:,3))]))*(cosTuningCurve-.5)*2;
%         Comp_BetWith_weights(:,4)=mean(abs([max(Comp_BetWith_weights(:,4)) min(Comp_BetWith_weights(:,4))]))*(cosTuningCurve2-.5)*2;
%         [Comp_testerr, Comp_testshift] = Universal_Test(1, itpos(testrange), ihd(testrange), AV(testrange), Comp_Spikes, Comp_G, Comp_BW_weights, Comp_BetWith_weights);
%         cf_Comp_testerr=[cf_Comp_testerr; Comp_testerr(:)'];
%         cf_Comp_testshift= [cf_Comp_testshift; Comp_testshift(:)'];
%         
%     figure(302); subplot(4,4,ps); plot(-31:2:31,Comp_BetWith_weights(:,5)); axis tight;
%     drawnow;
    
    
        shiftspercycle=[1:4 1:4 1:4 1:4 1:4 1:4 1:4 1:3];
        %shiftspercycle=[1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 ];
        ThetaSpikes = PlacePhaseBurstGenerator(ihd(trainrange), itpos(trainrange), ncells, shiftspercycle, basefrequency, .12); %create spike trains
        ThetaG = Generate_Conductance(itpos(trainrange), ncells, ThetaSpikes,tau); %convolve with exponential decay kernel
        [Theta_BW_weights, Theta_BetWith_weights, Theta_trainerr, Theta_trainshift] = Universal_Train(1, itpos(trainrange), ihd(trainrange), AV(trainrange), ThetaSpikes, ThetaG);
    
        ThetaSpikes = PlacePhaseBurstGenerator(ihd(testrange), itpos(testrange), ncells, shiftspercycle, basefrequency, .12); %create spike trains
        ThetaG = Generate_Conductance(itpos(testrange), ncells, ThetaSpikes,tau); %convolve with exponential decay kernel
        [Theta_testerr, Theta_testshift] = Universal_Test(1, itpos(testrange), ihd(testrange), AV(testrange), ThetaSpikes, ThetaG, Theta_BW_weights, Theta_BetWith_weights);
        cf_Theta_testerr=[cf_Theta_testerr; Theta_testerr];
        cf_Theta_testshift= [cf_Theta_testshift; Theta_testshift];
    
end

% figure(1); clf; 
% subplot(1,2,2);
% plot([.33 .66 1 1.33 1.66]',cf_CV1_testerr(:,3));
% hold on; plot([.33 .66 1 1.33 1.66]',cf_CV0_testerr(:,3));
% hold on; plot([.33 .66 1 1.33 1.66]',cf_Comp_testerr(:,3));
% subplot(1,2,2);
% plot(Count_CV1',cf_CV1_testshift(:,3));
% hold on; plot(Count_CV0',cf_CV0_testshift(:,3));
% hold on; plot(Count_Comp',cf_Comp_testshift(:,3));


% figure(2); clf; 
% subplot(1,2,2);
% plot([.2:.2:1.6]',cf_CV1_testerr(:,5));
% hold on; plot([.2:.2:1.6]',cf_CV0_testerr(:,5));
% hold on; plot([.2:.2:1.6]',cf_Comp_testerr(:,5));
% subplot(1,2,2);
% plot([.33 .66 1 1.33 1.66]',cf_CV1_testshift(:,5));
% hold on; plot([.33 .66 1 1.33 1.66]',cf_CV0_testshift(:,5));
% hold on; plot([.33 .66 1 1.33 1.66]',cf_Comp_testshift(:,5));

% figure(1); clf;
% subplot(1,2,2);
% plot([.2:.2:1.6]',cf_CV1_testerr(:,3));
% hold on; plot([.2:.2:1.6]',cf_CV0_testerr(:,3));
% hold on; plot([.2:.2:1.6]',cf_Comp_testerr(:,3));
% figure(2); 
% subplot(1,2,2);
% plot(.33*[.2 .4 .6 .8 1]',cf_CV1_testerr(:,5));
% hold on; plot(.33*[.2 .4 .6 .8 1]',cf_CV0_testerr(:,5));
% hold on; plot(.33*[.2 .4 .6 .8 1]',cf_Comp_testerr(:,5));

%[Theta_testerr, Theta_testshift] = Universal_Test(1, itpos((teststart-length(testrange)):testend), ihd((teststart-length(testrange)):testend), AV((teststart-length(testrange)):testend), ThetaSpikes, ThetaG, Theta_BW_weights, Theta_BetWith_weights);