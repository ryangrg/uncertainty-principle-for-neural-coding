function [fweights, fweights2, besterr, bestshift] = Universal_Train(behdata, veldata, gMatrixWithin, gMatrixBetween)



% behdata = behavior data to decode, sampled at 1 KHz (e.g., HD data)
% veldata = time deriative of behavior data (e.g., AV data)
% gMatrixWithin = conductance matrix for within cell firing rates (one per neuron)
% gMatrixBetween = conductance matrix for between cell firing rates (one per neuron)

%%%%%%%%%%%%%% LOOP THROUGH TIMS SHIFTS TO FIGURE OUT OPTIMAL PREDICTION DELAY


minerror_betwith=10000000;
minerror_between=10000000;
minerror_within=10000000;

minerrorav_betwith=10000000;
minerrorav_between=10000000;
minerrorav_within=10000000;

edex=1;
figure(3); clf;
figure(4); clf;

lowtime=-250;
inctime=25;
hitime=250;

for sh=lowtime:inctime:hitime %
    
    allvel=1001:(length(veldata)-1000); %indices for the conductances with 1st and last second truncated
    
    [x, y]=pol2cart(behdata,1); %separate sin & cos components of angular behdata (x is sin, y is cos)
    
    %compute weights to predict from both between & within
    weights(:,1) = pinv([gMatrixBetween(allvel,:) gMatrixWithin(allvel,:)])*x(allvel-sh);
    weights(:,2) = pinv([gMatrixBetween(allvel,:) gMatrixWithin(allvel,:)])*y(allvel-sh);
    weights(:,3) = pinv([gMatrixBetween(allvel,:) gMatrixWithin(allvel,:)])*veldata(allvel-sh);
    %use weights to predict from both between & within
    predx_betwith = sum( ([gMatrixBetween(allvel,:) gMatrixWithin(allvel,:)].*repmat(weights(:,1), 1, size([gMatrixBetween(allvel,:) gMatrixWithin(allvel,:)],1))')')';
    predy_betwith = sum( ([gMatrixBetween(allvel,:) gMatrixWithin(allvel,:)].*repmat(weights(:,2), 1, size([gMatrixBetween(allvel,:) gMatrixWithin(allvel,:)],1))')')';
    [predpos_betwith,r]=cart2pol(predx_betwith,predy_betwith);
    predpos_betwith(find(predpos_betwith<0))=predpos_betwith(find(predpos_betwith<0))+2*pi;
    predAV_betwith = sum( ([gMatrixBetween(allvel,:) gMatrixWithin(allvel,:)].*repmat(weights(:,3), 1, size([gMatrixBetween(allvel,:) gMatrixWithin(allvel,:)],1))')')';
    
    MSE_betwith(edex) = sum(circ_dist(predpos_betwith, behdata(allvel-sh)).^2)/length(allvel);
    if MSE_betwith(edex)<minerror_betwith
        minerror_betwith=MSE_betwith(edex);
        fweights(:,1:2)=weights(:,1:2);
        bestshift(1,1)=sh;
        besterr(1,1)=MSE_betwith(edex);
        figure(3);
        subplot(3,1,1);
        hold off; plot(behdata); hold on; plot(allvel,predpos_betwith);
    end

    MSEav_betwith(edex) = sum(circ_dist(predAV_betwith, veldata(allvel-sh)).^2)/length(allvel);
    if MSEav_betwith(edex)<minerrorav_betwith
        minerrorav_betwith=MSEav_betwith(edex);
        fweights(:,3)=weights(:,3);
        bestshift(1,2)=sh;
        besterr(1,2)=MSEav_betwith(edex);
        figure(4);
        subplot(3,1,1);
        hold off; plot(veldata); hold on; plot(allvel,predAV_betwith);
    end
    
    %compute weights to predict from between only
    weights2(:,1) = pinv(gMatrixBetween(allvel,:))*x(allvel-sh);
    weights2(:,2) = pinv(gMatrixBetween(allvel,:))*y(allvel-sh);
    weights2(:,5) = pinv(gMatrixBetween(allvel,:))*veldata(allvel-sh);
    %use weights to predict from between only
    predx_between = sum( (gMatrixBetween(allvel, :).*repmat(weights2(:,1), 1, size(gMatrixBetween(allvel,:),1))')')';
    predy_between = sum( (gMatrixBetween(allvel, :).*repmat(weights2(:,2), 1, size(gMatrixBetween(allvel,:),1))')')';
    [predpos_between,r]=cart2pol(predx_between,predy_between);
    predpos_between(find(predpos_between<0))=predpos_between(find(predpos_between<0))+2*pi;
    predAV_between = sum( (gMatrixBetween(allvel, :).*repmat(weights2(:,5), 1, size(gMatrixBetween(allvel,:),1))')')';
    
    MSE_between(edex) = sum(circ_dist(predpos_between, behdata(allvel-sh)).^2)/length(allvel);
    if MSE_between(edex)<minerror_between
        minerror_between=MSE_between(edex);
        fweights2(:,1:2)=weights2(:,1:2);
        bestshift(2,1)=sh;
        besterr(2,1)=MSE_between(edex);
        figure(3);
        subplot(3,1,2);
        hold off; plot(behdata); hold on; plot(allvel,predpos_between);
    end

    MSEav_between(edex) = sum(circ_dist(predAV_between, veldata(allvel-sh)).^2)/length(allvel);
    if MSEav_between(edex)<minerrorav_between
        minerrorav_between=MSEav_between(edex);
        fweights2(:,5)=weights2(:,5);
        bestshift(2,2)=sh;
        besterr(2,2)=MSEav_between(edex);
        figure(4);
        subplot(3,1,2);
        hold off; plot(veldata); hold on; plot(allvel,predAV_between);
    end
    
    %compute weights to predict from within only
    weights2(:,3) = pinv(gMatrixWithin(allvel,:))*x(allvel-sh);
    weights2(:,4) = pinv(gMatrixWithin(allvel,:))*y(allvel-sh);
    weights2(:,6) = pinv(gMatrixWithin(allvel,:))*veldata(allvel-sh);
    %use weights to predict from within only
    predx_within = sum( (gMatrixWithin(allvel, :).*repmat(weights2(:,3), 1, size(gMatrixWithin(allvel,:),1))')')';
    predy_within = sum( (gMatrixWithin(allvel, :).*repmat(weights2(:,4), 1, size(gMatrixWithin(allvel,:),1))')')';
    [predpos_within,r]=cart2pol(predx_within,predy_within);
    predpos_within(find(predpos_within<0))=predpos_within(find(predpos_within<0))+2*pi;
    predAV_within = sum( (gMatrixWithin(allvel, :).*repmat(weights2(:,6), 1, size(gMatrixWithin(allvel,:),1))')')';
    
    MSE_within(edex) = sum(circ_dist(predpos_within, behdata(allvel-sh)).^2)/length(allvel);
    if MSE_within(edex)<minerror_within
        minerror_within=MSE_within(edex);
        fweights2(:,3:4)=weights2(:,3:4);
        bestshift(3,1)=sh;
        besterr(3,1)=MSE_within(edex);
        figure(3);
        subplot(3,1,3);
        hold off; plot(behdata); hold on; plot(allvel,predpos_within);
    end
 
    MSEav_within(edex) = sum(circ_dist(predAV_within, veldata(allvel-sh)).^2)/length(allvel);
    if MSEav_within(edex)<minerrorav_within
        minerrorav_within=MSEav_within(edex);
        fweights2(:,6)=weights2(:,6);
        bestshift(3,2)=sh;
        besterr(3,2)=MSEav_within(edex);
        figure(4);
        subplot(3,1,3);
        hold off; plot(veldata); hold on; plot(allvel,predAV_within);
    end

 edex=edex+1;
 drawnow;


end

figure(101); clf; 
subplot(2,1,1);
hold on; plot(lowtime:inctime:hitime,MSE_betwith);
hold on; plot(lowtime:inctime:hitime,MSE_between);
plot(lowtime:inctime:hitime,MSE_within);
subplot(2,1,2);
hold on; plot(lowtime:inctime:hitime,MSEav_betwith);
hold on; plot(lowtime:inctime:hitime,MSEav_between);
plot(lowtime:inctime:hitime,MSEav_within);


