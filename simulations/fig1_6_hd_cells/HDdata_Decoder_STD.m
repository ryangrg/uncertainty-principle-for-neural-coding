
clear;
load('HeadDirectionData.mat');
%hd=mod((2*pi)/60:(2*pi)/60:9000*2*pi/60,2*pi)';

cf_Theta_testerr=[]; cf_Theta_testshift= [];
cf_CV1_testerr=[]; cf_CV1_testshift= [];
cf_CV0_testerr=[]; cf_CV0_testshift= [];
cf_Comp_testerr=[]; cf_Comp_testshift= [];

figure(300); clf;
figure(301); clf;
figure(302); clf;
ncells=12; basefrequency=7; peakrate=.15;



[ cosTuningCurve ] = cosPiTuningCurve(ncells, 1);
cosTuningCurve = [cosTuningCurve((round(ncells/2)+1):end); cosTuningCurve(1:round(ncells/2))];
[ cosTuningCurve2 ] = cos90TuningCurve(ncells, 1);
cosTuningCurve2  = cosTuningCurve2(end:-1:1);

load('SeedSet_10x120.mat');

sigma = 0.5; %standard deviation of vin Mises HD tuning curve

t1 = 10;   %simulation start time (s)
t2 = 300;   %simulation end time (s)
dt = .001;  %simulaton time step

networkType = 'UC'; %uncompressed (for file naming)

tauInh = [1.3 0.565 0.305 0.2 0.138 0.099 0.077 0.0575 0.046 0.038 0.0325 0.0275 0.023 0.02 0.01655 0.015];

%for centerfreq=1:5
centerfreq=0;
for ps=1:15
%for vscale=1:4
    
 %   ps=8;
    sigma=0.5;
    peakrate = .01+(ps-1)*.025;%.005+(ps-1)*.0025;
%    peakrate = .005+(6-1)*.0025;
    ti=tauInh(ps);
    
    cf=.1;%*centerfreq;
    
    clearvars -except ti hd tpos vscale tau* ps centerfreq cf cf_* ncells basefrequency peakrate tau cosTuningCurve* seeds t1 t2 dt sigma networkType test* train* Count_*;
    tau_within=.1%*2^(centerfreq-1); 
    tau_between=tau_within;
    tau_std=.5;
    
    perc_std=.8*ps.^(-1.3);%.9-.0607*(ps-1)*2.33;
    
    %duration=9000
    %trunc=2500;
%     duration=15000
%     trunc=5000;
%     tpos=(1/30):(1/30):duration;
%     av=rand(30*duration,1)*2-1;
% %    dfh = designfilt('bandpassiir', 'StopbandFrequency1', cf*.1, 'PassbandFrequency1', cf*.2, 'PassbandFrequency2', cf*1.8, 'StopbandFrequency2', cf*1.9, 'StopbandAttenuation1', 30, 'PassbandRipple', 1, 'StopbandAttenuation2', 30, 'SampleRate', 30);
%     dfh = designfilt('bandpassiir', 'StopbandFrequency1', .01, 'PassbandFrequency1', .02, 'PassbandFrequency2', 2*1.8, 'StopbandFrequency2', 2*1.9, 'StopbandAttenuation1', 30, 'PassbandRipple', 1, 'StopbandAttenuation2', 30, 'SampleRate', 30);
%     bpav=filtfilt(dfh,av);
%     figure(1); clf;
%     subplot(2,1,1);
%     plot(tpos((trunc*30:(duration-trunc)*30)),av((trunc*30:(duration-trunc)*30))); hold on; plot(tpos((trunc*30:(duration-trunc)*30)),bpav((trunc*30:(duration-trunc)*30)));
%     subplot(2,1,2);
%     hd=mod(cumsum(bpav(trunc*30:(duration-trunc)*30)),2*pi);
%     plot(tpos((trunc*30:(duration-trunc)*30)),hd);
%     tpos=tpos((trunc*30:(duration-trunc)*30));
%     tpos=tpos-tpos(1);
    %figure(2); clf; scatter(bpav(trunc*30:(duration-trunc)*30),hd);
    
    drawnow;
    tic;
    %networkSize = 128; %number of neurons in the HD ring
    
    % Upsampled head direction data to compare with predicted head direction
    %
    ihd(:,1) = interpolateHDData(hd, tpos, dt, t1, t2);
    % balance clk an ccw turning in the training set for symmetric velocity weights
    trainstart=18421;
    splitpoint=165930;
    trainend=splitpoint+(splitpoint-trainstart);
    trainrange=trainstart:trainend;
        
    ihd=[ihd(1:(trainstart-1)); ihd(trainstart:splitpoint); 2*pi-ihd(trainstart:splitpoint); ihd((splitpoint+1):end)];
    %ihd(647564)=mean([ihd(647563) ihd(647565)]); %AV(636744)=AV(636743);
    %trainend=13334+66667*9;
    teststart=trainend+1;
    testend=length(ihd)-3000;
    testrange=teststart:testend;
    
    % interpolate tpos using interp1 takes orginal signal and interpolates into steps
    % determined by dt.
    %
    itpos(:,1) = dt*[0:length(ihd)-1]'; %(interp1(tpos(1:end), tpos(1:end), 0:dt:((t2 - t1) - dt))';
    
    % variable to cut out conductance warmup from predictedHD and meanSpikes,
    % time will be cut from the start of itpos. Enter varialbe in seconds
    % desired to be cut.
    %cutTime = 1;
    %cutTime = cutTime*1000+1;
    
    % Compute angular velocity as circular difference of HD
    AV = circ_dist(ihd(2:end), ihd(1:end-1))*1000; 
    AV(splitpoint)=mean([AV(splitpoint-1) AV(splitpoint+1)]);
    AV(trainend+1)=mean([AV(trainend) AV(trainend+2)]);
%    AV(325039:(325039+311705))=-AV(325039:-1:trainstart);
    
    
    numseeds = 1;
    probScale(1, :) = [25:2.5:150];
    
    
    seedSpikes = zeros(numseeds, 1);        %total number of spikes, M, for a given seed
    seedMSE = zeros(numseeds, 1);           %mean squared error for predicted AV
    M_Spikes = zeros(numseeds, length(probScale));
    M_MSE = zeros(numseeds, length(probScale));
      
    
    seed=1;
    Spikes = CV0spikeGenerator(ihd(trainrange), itpos(trainrange), ncells, peakrate, sigma); %create spike trains
    Gstd = Generate_STD(itpos(trainrange), ncells, Spikes, tau_std, perc_std); %convolve with exponential decay kernel
    stdspk=Spikes-Gstd;
    stdspk(find(stdspk<0))=0;
    withinG = Generate_Conductance(itpos(trainrange), ncells, Spikes, tau_within); %convolve with exponential decay kernel
    betweenG = Generate_Between_STD(itpos(trainrange), Spikes, withinG, Gstd, tau_between);
    %betweenG = Generate_Between(itpos(trainrange), Spikes, withinG, tau_between);
    %[CV0_BW_weights, CV0_BetWith_weights, CV0_trainerr, CV0_trainshift] = Universal_Train(1, itpos(trainrange), ihd(trainrange), AV(trainrange), CV0_Spikes, CV0_withinG);
%    withinGstd = Generate_Conductance(itpos(trainrange), ncells, Spikes-Gstd, tau_within); %convolve with exponential decay kernel
    [CV0_BW_weights, CV0_BetWith_weights, CV0_trainerr, CV0_trainshift] = Universal_Train(ihd(trainrange), AV(trainrange), withinG, betweenG)    
    
    Spikes = CV0spikeGenerator(ihd(testrange), itpos(testrange), ncells, peakrate, sigma); %create spike trains
    Count_CV0(ps)=sum(Spikes(:));
    Gstd = Generate_STD(itpos(testrange), ncells, Spikes, tau_std, perc_std); %convolve with exponential decay kernel
    stdspk=Spikes-Gstd;
    stdspk(find(stdspk<0))=0;
    withinG = Generate_Conductance(itpos(testrange), ncells, Spikes, tau_within); %convolve with exponential decay kernel
    betweenG = Generate_Between_STD(itpos(testrange), Spikes, withinG, Gstd, tau_between);
%    withinGstd = Generate_Conductance(itpos(testrange), ncells, Spikes-Gstd, tau_within); %convolve with exponential decay kernel
%    betweenG = Generate_Between(itpos(testrange), Spikes, withinG, tau_between);
    %scale pure cosine weights to match scale of fitted weights
    CV0_BetWith_weights(:,3)=mean(abs([max(CV0_BetWith_weights(:,3)) min(CV0_BetWith_weights(:,3))]))*(cosTuningCurve-.5)*2;
    CV0_BetWith_weights(:,4)=mean(abs([max(CV0_BetWith_weights(:,4)) min(CV0_BetWith_weights(:,4))]))*(cosTuningCurve2-.5)*2;
    [CV0_testerr, CV0_testshift] = Universal_Test(ihd(testrange), AV(testrange), withinG, betweenG, CV0_BW_weights, CV0_BetWith_weights);
    cf_CV0_testerr=[cf_CV0_testerr; CV0_testerr(:)'];
    cf_CV0_testshift= [cf_CV0_testshift; CV0_testshift(:)'];
    
%    figure(300); subplot(4,4,centerfreq); plot(-31:2:31,CV0_BetWith_weights(:,5)); axis tight;
%    figure(300); subplot(4,4,ps); plot(-31:2:31,CV0_BetWith_weights(:,5)); axis tight;
%     drawnow;
    
        seed=1;
        Spikes = vonMisesFiringProb(ihd(trainrange), ncells, peakrate, sigma, seed); %create spike trains
    Gstd = Generate_STD(itpos(trainrange), ncells, Spikes, tau_std, perc_std); %convolve with exponential decay kernel
    stdspk=Spikes-Gstd;
    stdspk(find(stdspk<0))=0;
%        withinG = Generate_Conductance(itpos(trainrange), ncells, stdspk, tau_within); %convolve with exponential decay kernel
        withinG = Generate_Conductance_STD(itpos(trainrange), ncells, Spikes, tau_within); %convolve with exponential decay kernel
    betweenG = Generate_Between_STD(itpos(trainrange), Spikes, withinG, Gstd, tau_between);
%        betweenG = Generate_Between(itpos(trainrange), Spikes, withinG, tau_between);
    %withinGstd = Generate_Conductance(itpos(trainrange), ncells, Spikes-Gstd, tau_within); %convolve with exponential decay kernel
        [CV1_BW_weights, CV1_BetWith_weights, CV1_trainerr, CV1_trainshift] = Universal_Train(ihd(trainrange), AV(trainrange), withinG, betweenG)    
    
        Spikes = vonMisesFiringProb(ihd(testrange), ncells, peakrate, sigma, seed); %create spike tests
        Count_CV1(ps)=sum(Spikes(:));
    Gstd = Generate_STD(itpos(testrange), ncells, Spikes, tau_std, perc_std); %convolve with exponential decay kernel
    stdspk=Spikes-Gstd;
    stdspk(find(stdspk<0))=0;
        withinG = Generate_Conductance_STD(itpos(testrange), ncells, Spikes, tau_within); %convolve with exponential decay kernel
%        withinGstd = Generate_Conductance_STD(itpos(testrange), ncells, Spikes, tau_within); %convolve with exponential decay kernel
    betweenG = Generate_Between_STD(itpos(testrange), Spikes, withinG, Gstd, tau_between);
%    withinGstd = Generate_Conductance(itpos(testrange), ncells, Spikes, tau_within); %convolve with exponential decay kernel
%        betweenG = Generate_Between(itpos(testrange), Spikes, withinG, tau_between);
        %scale pure cosine weights to match scale of fitted weights
        CV1_BetWith_weights(:,3)=mean(abs([max(CV1_BetWith_weights(:,3)) min(CV1_BetWith_weights(:,3))]))*(cosTuningCurve-.5)*2;
        CV1_BetWith_weights(:,4)=mean(abs([max(CV1_BetWith_weights(:,4)) min(CV1_BetWith_weights(:,4))]))*(cosTuningCurve2-.5)*2;
        [CV1_testerr, CV1_testshift] = Universal_Test(ihd(testrange), AV(testrange), withinG, betweenG, CV1_BW_weights, CV1_BetWith_weights);
        cf_CV1_testerr=[cf_CV1_testerr; CV1_testerr(:)'];
        cf_CV1_testshift= [cf_CV1_testshift; CV1_testshift(:)'];
        
%     figure(301); subplot(4,4,ps); plot(-31:2:31,CV1_BetWith_weights(:,5)); axis tight;
%     drawnow;

%         seed=1;
%         Spikes = vonMisesFiringProb(ihd(trainrange), ncells, .15, sigma, seed); %create spike trains
%         Spikes = SpikePrune(itpos(trainrange), ncells, Spikes, ti, 2);
%         withinG = Generate_Conductance(itpos(trainrange), ncells, Spikes,tau_within); %convolve with exponential decay kernel
%         betweenG = Generate_Between(itpos(trainrange), Spikes, withinG, tau_between);
%         [Comp_BW_weights, Comp_BetWith_weights, Comp_trainerr, Comp_trainshift] = Universal_Train(ihd(trainrange), AV(trainrange), withinG, betweenG)    
%     
%         Spikes = vonMisesFiringProb(ihd(testrange), ncells, .15, sigma, seed); %create spike tests
%         Spikes = SpikePrune(itpos(testrange), ncells, Spikes, ti, 2);
%         Count_Comp(ps)=sum(Spikes(:));
%         withinG = Generate_Conductance(itpos(testrange), ncells, Spikes,tau_within); %convolve with exponential decay kernel
%         betweenG = Generate_Between(itpos(testrange), Spikes, withinG, tau_between);
%         %scale pure cosine weights to match scale of fitted weights
%         Comp_BetWith_weights(:,3)=mean(abs([max(Comp_BetWith_weights(:,3)) min(Comp_BetWith_weights(:,3))]))*(cosTuningCurve-.5)*2;
%         Comp_BetWith_weights(:,4)=mean(abs([max(Comp_BetWith_weights(:,4)) min(Comp_BetWith_weights(:,4))]))*(cosTuningCurve2-.5)*2;
%         [Comp_testerr, Comp_testshift] = Universal_Test(ihd(testrange), AV(testrange), withinG, betweenG, Comp_BW_weights, Comp_BetWith_weights);
%         cf_Comp_testerr=[cf_Comp_testerr; Comp_testerr(:)'];
%         cf_Comp_testshift= [cf_Comp_testshift; Comp_testshift(:)'];
%         
%     figure(302); subplot(4,4,ps); plot(-31:2:31,Comp_BetWith_weights(:,5)); axis tight;
%     drawnow;
    




% tau_z=1600;
%     ckern=gausswin(tau_z); ckern=ckern-min(ckern); 
%     ckern=ckern.*cos(2*pi*basefrequency*[(-tau_z/2):(tau_z/2-1)]/1000)'; 
%     ckern(find(ckern<0))=0;
%     expkern=exp(-[.001:.001:.25]/.02); expkern=[expkern(end:-1:1) expkern];
%     ckern=conv(ckern,expkern,'same');
%     ckern=ckern((tau_z/2+1):tau_z); ckern=ckern/sum(ckern);
%     figure(1); clf; plot(ckern);
    
%         shiftspercycle=[1:4 1:4 1:4 1:4 1:4 1:4 1:4 1:3];
%         %shiftspercycle=[1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 5 1 1.7 2.9 ];
%         Spikes = PlacePhaseBurstGenerator(ihd(trainrange), itpos(trainrange), ncells, shiftspercycle, basefrequency, peakrate/5); %create spike trains
%         Count_Theta(ps)=sum(Spikes(:));
%         withinG = Generate_Conductance(itpos(trainrange), ncells, Spikes,tau_within); %convolve with exponential decay kernel
% %     clear withinGx;
% %     for i=1:size(withinG,2)        
% %        withinGx(:,i) = conv(double(Spikes(:,i)),ckern);        
% %     end
% %    withinGx = withinGx(1:end-(tau_z/2-1),:);        
% %        betweenG = Generate_Between(itpos(trainrange), Spikes, withinG, tau_between);
%         betweenGx = Generate_Between(itpos(trainrange), Spikes, withinGx, tau_between);
% %        [Theta_BW_weights, Theta_BetWith_weights, Theta_trainerr, Theta_trainshift] = Universal_Train(ihd(trainrange), AV(trainrange), withinG, betweenG);    
%         [Theta_BW_weights_x, Theta_BetWith_weights_x, Theta_trainerr_x, Theta_trainshift_x] = Universal_Train(ihd(trainrange), AV(trainrange), withinG, betweenGx);    
%         
%         Spikes = PlacePhaseBurstGenerator(ihd(testrange), itpos(testrange), ncells, shiftspercycle, basefrequency, peakrate/5); %create spike trains
%         withinG = Generate_Conductance(itpos(testrange), ncells, Spikes,tau_within); %convolve with exponential decay kernel
% %     clear withinGx;
% %     for i=1:size(withinG,2)        
% %        withinGx(:,i) = conv(double(Spikes(:,i)),ckern);        
% %     end
% %    withinGx = withinGx(1:end-(tau_z/2-1),:);        
% %        betweenG = Generate_Between(itpos(testrange), Spikes, withinG, tau_between);
%         betweenGx = Generate_Between(itpos(testrange), Spikes, withinGx, tau_between);
% %        [Theta_testerr, Theta_testshift] = Universal_Test(ihd(testrange), AV(testrange), withinG, betweenG, Theta_BW_weights, Theta_BetWith_weights);
%         [Theta_testerr_x, Theta_testshift_x] = Universal_Test(ihd(testrange), AV(testrange), withinG, betweenGx, Theta_BW_weights_x, Theta_BetWith_weights_x);
%         cf_Theta_testerr=[cf_Theta_testerr; Theta_testerr_x(:)'];
%         cf_Theta_testshift= [cf_Theta_testshift; Theta_testshift_x(:)'];
    

% %%%% figure 1 plots position decoded from within rates
% figure(110+centerfreq); clf; 
% subplot(1,2,1);
% plot(Count_CV1',cf_CV1_testerr(:,3),'o-');
% hold on; plot(Count_CV0',cf_CV0_testerr(:,3),'o-');
% %hold on; plot(Count_Comp',cf_Comp_testerr(:,3),'o-');
% hold on; plot(Count_Theta',cf_Theta_testerr(:,2),'o-'); %theta decoded from between rather than within
% subplot(1,2,2);
% plot(Count_CV1',cf_CV1_testshift(:,3),'o-');
% hold on; plot(Count_CV0',cf_CV0_testshift(:,3),'o-');
% %hold on; plot(Count_Comp',cf_Comp_testshift(:,3),'o-');
% hold on; plot(Count_Theta',cf_Theta_testshift(:,2),'o-'); %theta decoded from between rather than within
% 
% 
% %%%% figure 2 plots velocity decoded from between rates
% figure(210+centerfreq); clf; 
% subplot(1,2,1);
% plot(Count_CV1',cf_CV1_testerr(:,5),'o-');
% hold on; plot(Count_CV0',cf_CV0_testerr(:,5),'o-');
% %hold on; plot(Count_Comp',cf_Comp_testerr(:,5),'o-');
% subplot(1,2,2);
% plot(Count_CV1',cf_CV1_testshift(:,5),'o-');
% hold on; plot(Count_CV0',cf_CV0_testshift(:,5),'o-');
% %hold on; plot(Count_Comp',cf_Comp_testshift(:,5),'o-');

drawnow;

end
%end