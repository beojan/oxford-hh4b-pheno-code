% Estimation of reducible background for the RESOLVED channel

close all
clear all

load('sampleHighStats') % Need to load HH 4-tag, and 4b 4-tag

% Global variables
luminosityScaling = 4000; %Integrated luminosity
numBins = 17; % For histogramming data
Plots = 1; % To produce the distribution plots of the samples at the beginning

% Channel Sort
[resHH4tag,intHH4tag,boostHH4tag] = channelSort(HH4tag);
[res4b4tag,int4b4tag,boost4b4tag] = channelSort(qcd4b4tag);

%Split into training and Samples
res4b4tagTrainWeight = res4b4tag.weight(1:length(res4b4tag.weight)/2,1); %Training (first half)
res4b4tagTrainmH0 = res4b4tag.m_H0(1:length(res4b4tag.weight)/2,1); %mH0 
res4b4tagTrainmH1 = res4b4tag.m_H1(1:length(res4b4tag.weight)/2,1); %mH1
res4b4tagSampleWeight = res4b4tag.weight(length(res4b4tag.weight)/2:length(res4b4tag.weight),1); %Sample to fit (testing) (second half)
res4b4tagSamplemH0 = res4b4tag.m_H0(length(res4b4tag.weight)/2:length(res4b4tag.weight),1); %mH0 
res4b4tagSamplemH1 = res4b4tag.m_H1(length(res4b4tag.weight)/2:length(res4b4tag.weight),1); %mH1

resHH4tagTrainWeight = resHH4tag.weight(1:length(resHH4tag.weight)/2,1); %Training (first half) %Weight
resHH4tagTrainmH0 = resHH4tag.m_H0(1:length(resHH4tag.weight)/2,1); %mH0 
resHH4tagTrainmH1 = resHH4tag.m_H1(1:length(resHH4tag.weight)/2,1); %mH1
resHH4tagSampleWeight = resHH4tag.weight(round(length(resHH4tag.weight)/2):length(resHH4tag.weight),1); %Sample to fit (testing) (second half)
resHH4tagSamplemH0 = resHH4tag.m_H0(round(length(resHH4tag.weight)/2):length(resHH4tag.weight),1); %mH0 
resHH4tagSamplemH1 = resHH4tag.m_H1(round(length(resHH4tag.weight)/2):length(resHH4tag.weight),1); %mH1

% Histograms (mH0, mH1) - 4b and HH (Training and Testing halves)
[bin4b4res,bin4b4Varres,binValue14b4res,binValue24b4res] = HistoBins(res4b4tagTrainWeight,res4b4tagTrainmH0,res4b4tagTrainmH1,luminosityScaling,numBins,Plots,'4b 4-tag Resolved Training Sample','m_H0','m_H1');
[binHH4res,binHH4Varres,binValue1HH4res,binValue2HH4res] = HistoBins(resHH4tagTrainWeight,resHH4tagTrainmH0,resHH4tagTrainmH1,luminosityScaling,numBins,Plots,'HH 4-tag Resolved Training Sample','m_H0','m_H1');
[binHH4res32,binHH4Varres32,binValue1HH4res32,binValue2HH4res32] = HistoBins(resHH4tagTrainWeight,resHH4tagTrainmH0,resHH4tagTrainmH1,luminosityScaling,33,Plots,'HH 4-tag Resolved Training Sample','m_H0','m_H1');

[bin4b4resSample,bin4b4VarresSample,binValue14b4resSample,binValue24b4resSample] = HistoBins(res4b4tagSampleWeight,res4b4tagSamplemH0,res4b4tagSamplemH1,luminosityScaling,numBins,Plots,'4b 4-tag Resolved Testing Sample','m_H0','m_H1');
[binHH4resSample,binHH4VarresSample,binValue1HH4resSample,binValue2HH4resSample] = HistoBins(resHH4tagSampleWeight,resHH4tagSamplemH0,resHH4tagSamplemH1,luminosityScaling,numBins,Plots,'HH 4-tag Resolved Testing Sample','m_H0','m_H1');

% 1D plots - no fit, to overlay the fit to later
histogram1Dboth(bin4b4res,sqrt(bin4b4Varres),binValue14b4res,binValue24b4res,'4b-4-tag-Resolved','m_H0','m_H1');
histogram1Dboth(binHH4res32,sqrt(binHH4Varres32),binValue1HH4res32,binValue2HH4res32,'HH-4-tag-Resolved','m_H0','m_H1');

% 1D fit:
% Preallocate
binValue1mid4b4res = zeros(numBins-1,1);
binValue2mid4b4res = zeros(numBins-1,1);
binValue1midHH4res = zeros(numBins-1,1);
binValue2midHH4res = zeros(numBins-1,1);
binValue1midHH4res32 = zeros(numBins-1,1);
binValue2midHH4res32 = zeros(numBins-1,1);

%Inverse Variance for weights in fits
invVar4b = 1./(bin4b4Varres); %4b
invVar4bmH0 = 1./(sum(bin4b4Varres,2)); % variance 4b mH0
invVar4bmH1 = 1./(sum(bin4b4Varres,1)); % variance 4b mH1
invVarHH = 1./(binHH4Varres); % variance HH 
invVarHHmH0 = 1./(sum(binHH4Varres32,2)); % HH mH0
invVarHHmH1 = 1./(sum(binHH4Varres32,1)); % HH mH1
invVarSum = 1./(binHH4VarresSample + bin4b4VarresSample); %Fit is for the sum of testing samples

% Sum of testing samples for final likelihood fit
binSumBothSample = bin4b4resSample + binHH4resSample; % For the final fit

%1D fits:
for i=2:length(binValue14b4res) %4b - 16 bins
    binValue1mid4b4res(i-1) = 0.5*(binValue14b4res(i-1) + binValue14b4res(i));
    binValue2mid4b4res(i-1) = 0.5*(binValue24b4res(i-1) + binValue24b4res(i));
end
for i=2:length(binValue1HH4res) %HH - 16 bins
    binValue1midHH4res(i-1) = 0.5*(binValue1HH4res(i-1) + binValue1HH4res(i));
    binValue2midHH4res(i-1) = 0.5*(binValue2HH4res(i-1) + binValue2HH4res(i));
end
for i=2:length(binValue1HH4res32) %HH - 32 bins
    binValue1midHH4res32(i-1) = 0.5*(binValue1HH4res32(i-1) + binValue1HH4res32(i));
    binValue2midHH4res32(i-1) = 0.5*(binValue2HH4res32(i-1) + binValue2HH4res32(i));
end

% HH 1D fits:
x0_mH0 = [1,5,125,7,100,1,2]; %Start point
optimal_x_HHmH0 = nlinfit(binValue1midHH4res32,sum(binHH4res32,2),@CrystalBall2,x0_mH0,'Weights',invVarHHmH0); %Fit to Crystal Ball function
output_mH0 = CrystalBall2(optimal_x_HHmH0,85:0.1:165); % Best Fit
x0_mH1 = [1,10,120,7,100,0.75,1.5]; % Start Point
optimal_x_HHmH1 = nlinfit(binValue2midHH4res32,sum(binHH4res32,1).',@CrystalBall2,x0_mH1,'Weights',invVarHHmH1.'); %Fit to Crystal Ball function
output_mH1 = CrystalBall2(optimal_x_HHmH1,85:0.1:165); % Best Fit

%Next 2 figures overlay the best fits to the histogrammed data
figure(13) %m_H1 HH fit
hold on
plot(85:0.1:165,output_mH1)
hold off
xlabel('m_{H1} / GeV')
ylabel('Number of Events')
title('HH 4-tag Resolved m_{H1} distribution')
legend('Sample','Error','Crystal Ball Fit')

figure(14) %m_H0 HH fit
hold on
plot(85:0.1:165,output_mH0)
hold off
xlabel('m_{H0} / GeV')
ylabel('Number of Events')
title('HH 4-tag Resolved m_{H0} distribution')
legend('Sample','Error','Crystal Ball Fit')

binHH4resFitmH0 = CrystalBall2(optimal_x_HHmH0,binValue1midHH4res); %16 bins, HH fit in mH0
binHH4resFitmH1 = CrystalBall2(optimal_x_HHmH1,binValue1midHH4res); %16 bins, HH fit in mH1
binHH4resFitmH032 = CrystalBall2(optimal_x_HHmH0,binValue1midHH4res32); %32 bins, HH fit in mH0
binHH4resFitmH132 = CrystalBall2(optimal_x_HHmH1,binValue1midHH4res32); %32 bins, HH fit in mH1

%Normalise for fit later:
N = sum(sum(binSumBothSample)); %Total number of events in testing sample
%mH0 HH
optimal_x_HHmH0(5) = optimal_x_HHmH0(5)/sum(binHH4resFitmH0); %scaled to 1
%mH1 HH
optimal_x_HHmH1(5) = optimal_x_HHmH1(5)/sum(binHH4resFitmH1); %scaled to 1

% Combining the 1D fits to 2D distribution for HH signal
fitHH4res = zeros(numBins-1,numBins-1);
for i=1:(numBins-1)
    for j=1:(numBins-1)
        fitHH4res(i,j) = CrystalBall2(optimal_x_HHmH0,binValue1midHH4res(i))*CrystalBall2(optimal_x_HHmH1,binValue2midHH4res(j));
    end
end

% 2D fit (4b and HH included):
% 2D Gaussian fit:
ftGauss2 = fittype('a*exp(-(b*(mH0-c)^2 + 2*d*(mH0-c)*(mH1-e) + f*(mH1-e)^2)) + g','independent',{'mH0','mH1'},'coefficients',{'a','b','c','d','e','f','g'});

%Set up variables for fits:
n=0;
fitValues4b=zeros((numBins-1)^2,4); % 4b background
fitValuesHH=zeros((numBins-1)^2,4); % HH signal
fitValuesBothSample=zeros((numBins-1)^2,4); % Testing sample (for a chi2 fit, not viable for current size)
for i=1:(numBins-1)
    for j=1:(numBins-1)
        n=n+1;
        fitValues4b(n,:) = [binValue1mid4b4res(i),binValue2mid4b4res(j),bin4b4res(i,j),invVar4b(i,j)];
        fitValuesHH(n,:) = [binValue1midHH4res(i),binValue2midHH4res(j),binHH4res(i,j),invVarHH(i,j)];
        fitValuesBothSample(n,:) = [binValue1mid4b4res(i),binValue2mid4b4res(j),binSumBothSample(i,j),invVarSum(i,j)];
    end
end

% The 2D fits for 4b and HH:
[fitGauss2D4b4res,out4b2e1,out4b2d2] = fit([fitValues4b(:,1),fitValues4b(:,2)],fitValues4b(:,3),ftGauss2,'weight',fitValues4b(:,4),'startpoint',[5e3,1e-3,125,-1e-3,125,1e-4,1e2]); %[1e4,1e-3,125,-1e-3,125,1e-4,1e3] [1e-2,1e-3,125,-1e-3,125,1e-4,1e-3]
[fitGauss2DHH4res,outhh2d1,outhh2d2] = fit([fitValuesHH(:,1),fitValuesHH(:,2)],fitValuesHH(:,3),ftGauss2,'weight',fitValuesHH(:,4),'startpoint',[100,1e-3,125,1e-3,125,1e-3,1] ); %[100,1e-3,125,1e-3,125,1e-3,1] [5e-2,5e-3,125,-1e-3,125,5e-3,1e-5]

% Histogram of fits to compare to sample (for a chi2 value)
fit4b4res2D = zeros(numBins-1,numBins-1);
fitHH4res2D = zeros(numBins-1,numBins-1);
for i=1:(numBins-1)
    for j=1:(numBins-1)
        fit4b4res2D(i,j) = fitGauss2D4b4res(binValue1mid4b4res(i),binValue2mid4b4res(j)); %4b 4-tag resolved fitted dist
        fitHH4res2D(i,j) = fitGauss2DHH4res(binValue1midHH4res(i),binValue2midHH4res(j)); %HH 4-tag resolved fitted dist 2D gauss
    end
end

%This figure compares the 2D templates for signal and background
figure()
subplot(3,2,1)
histogram2('XBinEdges',binValue1HH4res,'YBinEdges',binValue2HH4res,'BinCounts',binHH4res,'FaceColor','flat')
title('HH Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,2,2)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',bin4b4res,'FaceColor','flat')
title('4b Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,2,3)
histogram2('XBinEdges',binValue1HH4res,'YBinEdges',binValue2HH4res,'BinCounts',fitHH4res,'FaceColor','flat')
title('HH fit - Crystal Ball Functions')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,2,5)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',fitHH4res2D,'FaceColor','flat')
title('HH fit - 2D gaussian with constant')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,2,6)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',fit4b4res2D,'FaceColor','flat')
title('4b fit - 2D gaussian with constant')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')

% Fit the total events:
a = fitGauss2D4b4res.a/sum(sum(fit4b4res2D));
b = fitGauss2D4b4res.b;
c = fitGauss2D4b4res.c;
d = fitGauss2D4b4res.d;
e = fitGauss2D4b4res.e;
f = fitGauss2D4b4res.f;
g = fitGauss2D4b4res.g/sum(sum(fit4b4res2D));

%Chi2 fit of templates to a testing sample
ftSumBoth = fittype( @(alpha,mH0,mH1) N*((1-alpha).*(a.*exp(-(b.*(mH0-c).^2 + 2*d.*(mH0-c).*(mH1-e) + f.*(mH1-e).^2)) + g) + alpha.*(  optimal_x_HHmH0(5).*(heaviside((optimal_x_HHmH0(3)-optimal_x_HHmH0(1).*optimal_x_HHmH0(4)) - mH0).*(( optimal_x_HHmH0(2) ./ abs(optimal_x_HHmH0(1))).^(optimal_x_HHmH0(2)).*exp(-(abs(optimal_x_HHmH0(1)).^2)/2) .* (optimal_x_HHmH0(2)./abs(optimal_x_HHmH0(1)) - abs(optimal_x_HHmH0(1)) - (mH0 - optimal_x_HHmH0(3))./optimal_x_HHmH0(4)).^-optimal_x_HHmH0(2)) + (heaviside(mH0 - (optimal_x_HHmH0(3)-optimal_x_HHmH0(1).*optimal_x_HHmH0(4))).*heaviside((optimal_x_HHmH0(3)+optimal_x_HHmH0(7).*optimal_x_HHmH0(4)) - mH0).*(exp(-(mH0-optimal_x_HHmH0(3)).^2/(2.*optimal_x_HHmH0(4).^2)))) + (heaviside(mH0 - (optimal_x_HHmH0(3)+optimal_x_HHmH0(7).*optimal_x_HHmH0(4))).*(( optimal_x_HHmH0(6) ./ abs(optimal_x_HHmH0(7))).^(optimal_x_HHmH0(6)).*exp(-(abs(optimal_x_HHmH0(7)).^2)/2) .* (optimal_x_HHmH0(6)./abs(optimal_x_HHmH0(7)) - abs(optimal_x_HHmH0(7)) - (optimal_x_HHmH0(3) - mH0)./optimal_x_HHmH0(4)).^-optimal_x_HHmH0(6)))) ).*(  optimal_x_HHmH1(5).*(heaviside((optimal_x_HHmH1(3)-optimal_x_HHmH1(1).*optimal_x_HHmH1(4)) - mH1).*(( optimal_x_HHmH1(2) ./ abs(optimal_x_HHmH1(1))).^(optimal_x_HHmH1(2)).*exp(-(abs(optimal_x_HHmH1(1)).^2)/2) .* (optimal_x_HHmH1(2)./abs(optimal_x_HHmH1(1)) - abs(optimal_x_HHmH1(1)) - (mH1 - optimal_x_HHmH1(3))./optimal_x_HHmH1(4)).^-optimal_x_HHmH1(2)) + (heaviside(mH1 - (optimal_x_HHmH1(3)-optimal_x_HHmH1(1).*optimal_x_HHmH1(4))).*heaviside((optimal_x_HHmH1(3)+optimal_x_HHmH1(7).*optimal_x_HHmH1(4)) - mH1).*(exp(-(mH1-optimal_x_HHmH1(3)).^2/(2.*optimal_x_HHmH1(4).^2)))) + (heaviside(mH1 - (optimal_x_HHmH1(3)+optimal_x_HHmH1(7).*optimal_x_HHmH1(4))).*(( optimal_x_HHmH1(6) ./ abs(optimal_x_HHmH1(7))).^(optimal_x_HHmH1(6)).*exp(-(abs(optimal_x_HHmH1(7)).^2)/2) .* (optimal_x_HHmH1(6)./abs(optimal_x_HHmH1(7)) - abs(optimal_x_HHmH1(7)) - (optimal_x_HHmH1(3) - mH1)./optimal_x_HHmH1(4)).^-optimal_x_HHmH1(6)))))),'independent',{'mH0','mH1'});
[fittedBoth,gofitBoth,outputFitBoth] = fit([fitValuesBothSample(:,1),fitValuesBothSample(:,2)],fitValuesBothSample(:,3),ftSumBoth,'weight',fitValuesBothSample(:,4),'startpoint',1e-4,'upper',1,'lower',0);

%Histogram of chi2 template fit
binFittedValuesSample = zeros(numBins-1,numBins-1);
for i=1:(numBins-1)
    for j=1:(numBins-1)
        binFittedValuesSample(i,j) = fittedBoth(binValue1mid4b4res(i),binValue2mid4b4res(j));
    end
end

%Figure contains: Testing sample, Fit to testing sample with templates, and
%the two templates used (4b and HH)
figure()
subplot(2,2,1)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',binSumBothSample,'FaceColor','flat')
title('Sum of HH and 4B Resolved Testing Sample')
colorbar
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,2)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',binFittedValuesSample,'FaceColor','flat')
title('Fitted Distribution of Resolved Testing Sample')
colorbar
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,3)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',fit4b4res2D,'FaceColor','flat')
title('Fitted 4b training sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,4)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',fitHH4res2D,'FaceColor','flat')
title('Fitted HH training sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')

%Confidence intervals of chi2 fit at 1-4 sigma intervals on fraction of
%signal to background (alpha)
confidence = [confint(fittedBoth,0.6827),confint(fittedBoth,0.9545),confint(fittedBoth,0.9973),confint(fittedBoth,0.99993)];

%Scale HH fit:
fitHH4resScale = fitHH4res*sum(sum(binHH4resSample));
% Residual plots:
% HH, 4b, Sum
% ( ( Sample - Fit )/(SampleErr) )^2
residualsHH2D = (binHH4resSample -fitHH4resScale).^2./(binHH4VarresSample); % 1D fits
residuals4b2D = (bin4b4resSample -fit4b4res2D).^2./(bin4b4VarresSample); % 2D fit
residualsSum2D = (binSumBothSample -binFittedValuesSample).^2./(binHH4VarresSample + bin4b4VarresSample);
% Residuals for 1D fits:
residualsHHmH016 = (sum(binHH4res,2).' - 2*binHH4resFitmH0.').^2./(sum(binHH4Varres,2).');
residualsHHmH116 = (sum(binHH4res,1) - 2*binHH4resFitmH1.').^2./(sum(binHH4Varres,1));
residualsHHmH032 = (sum(binHH4res32,2).' - binHH4resFitmH032.').^2./(sum(binHH4Varres32,2).');
residualsHHmH132 = (sum(binHH4res32,1) - binHH4resFitmH132.').^2./(sum(binHH4Varres32,1));

% Figure contains: 2D residuals (HH, 4b and Testing Sample)
figure()
subplot(2,2,1)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',residualsHH2D,'FaceColor','flat')
title('Residuals HH 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')
subplot(2,2,2)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',residuals4b2D,'FaceColor','flat')
title('Residuals 4b 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')
subplot(2,2,3)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',residualsSum2D,'FaceColor','flat')
title('Residuals Signal+Background 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')


% Figure contains: each row: (HH, 4b, Testing), and then the Sample, fit
% and residuals for each respectively.
figure()
subplot(3,3,1)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',binHH4resSample,'FaceColor','flat')
title('HH 4-tag Resolved Testing Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,2)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',fitHH4resScale,'FaceColor','flat')
title('HH 4-tag Resolved 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,3)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',residualsHH2D,'FaceColor','flat')
title('Residuals HH 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')

subplot(3,3,4)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',bin4b4resSample,'FaceColor','flat')
title('4b 4-tag Resolved Testing Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,5)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',fit4b4res2D,'FaceColor','flat')
title('4b 4-tag Resolved 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,6)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',residuals4b2D,'FaceColor','flat')
title('Residuals 4b 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')

subplot(3,3,7)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',binSumBothSample,'FaceColor','flat')
title('Signal+Background 4-tag Resolved Testing Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,8)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',binFittedValuesSample,'FaceColor','flat')
title('Signal+Background 4-tag Resolved 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,9)
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',residualsSum2D,'FaceColor','flat')
title('Residuals Signal+Background 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')


%Chi2 values of each 2D fit
chi24b = sum(sum(residuals4b2D));
chi2HH = sum(sum(residualsHH2D));
chi2Sum = sum(sum(residualsSum2D));

fprintf('\nManual Chi2 Values (2D):\n 4b: %s\n HH: %s\n Sum: %s\n',chi24b,chi2HH,chi2Sum)

%Chi2 values of each 1D fit (and scale for 16 bins)
chi2HHmH016 = sum(sum(residualsHHmH016));
chi2HHmH116 = sum(sum(residualsHHmH116));
chi2HHmH032 = sum(sum(residualsHHmH032));
chi2HHmH132 = sum(sum(residualsHHmH132));

fprintf('\nManual Chi2 Values (1D):\n HH mH0 (32 bins): %s\n HH mH1 (32 bins): %s\n HH mH0 (16 bins): %s\n HH mH1 (16 bins): %s\n',chi2HHmH032,chi2HHmH132,chi2HHmH016,chi2HHmH116)


%Manual calculate chi^2 for different alpha (fraction of HH): (Manual chi2
%fit)
alpha = 0;
fit4b4res2DNorm = fit4b4res2D/sum(sum(fit4b4res2D));
bin4b4resNorm = bin4b4res/sum(sum(bin4b4res));
binHH4resNorm = binHH4res/sum(sum(binHH4res));
for i=1:10000
    alpha = alpha + 0.00001;
    for j=1:(numBins-1)
        for k=1:(numBins-1)
            manualBinsSimple(j,k) = N*( (1-alpha)*fit4b4res2DNorm(j,k) + alpha*fitHH4res(j,k));
            alphaZeroCalc(j,k) = N*( (1-alpha)*bin4b4resNorm(j,k) + alpha*binHH4resNorm(j,k));
        end
    end
    residualsmanual2DSimple = (binSumBothSample - manualBinsSimple).^2./(binHH4VarresSample + bin4b4VarresSample);
    residualsmanualAlphaZero = (binSumBothSample - alphaZeroCalc).^2./(binHH4VarresSample + bin4b4VarresSample);
    chi2Manual(i,1) = sum(sum(residualsmanual2DSimple));
    chi2AlphaZero(i,1) = sum(sum(residualsmanualAlphaZero));
    alphaManual(i,1) = alpha;
end

figure()
plot(alphaManual,chi2Manual)
title('Chi2 v. Fraction diHiggs')
xlabel('alpha - fraction diHiggs')
ylabel('chi^{2}')
figure()
plot(alphaManual,chi2AlphaZero)
title('Chi2 v. Fraction diHiggs - Alpha Zero')
xlabel('alpha - fraction diHiggs')
ylabel('chi^{2}')

%cross-sections
fprintf('%s\n %s\n %s\n %s\n %s\n',sum(sum(res4b4tagTrainWeight)),sum(sum(resHH4tagTrainWeight)),sum(sum(fit4b4res2D))/4000,sum(sum(binHH4resFitmH032))/4000,sum(sum(binHH4resFitmH132))/4000)
xsec_HH_test = sum(sum(resHH4tagSampleWeight));
xsec_4b_test = sum(sum(res4b4tagSampleWeight));

%S/root(B)
radiusSignal = [7,14,21,28];
signalCount = zeros(1,length(radiusSignal));
backgroundCount = zeros(1,length(radiusSignal));

for k = 1:length(radiusSignal)
    for i=1:length(res4b4tagSampleWeight)
        if (res4b4tagSamplemH0(i) - optimal_x_HHmH0(3))^2 + (res4b4tagSamplemH1(i) - optimal_x_HHmH1(3))^2 < radiusSignal(k)^2 %Centered on each fit peak
            backgroundCount(k) = backgroundCount(k) + res4b4tagSampleWeight(i)*luminosityScaling;
        end
    end
    for i=1:length(resHH4tagSampleWeight)
        if (resHH4tagSamplemH0(i) - optimal_x_HHmH0(3))^2 + (resHH4tagSamplemH1(i) - optimal_x_HHmH1(3))^2 < radiusSignal(k)^2 %Centered on each fit peak
            signalCount(k) = signalCount(k) + resHH4tagSampleWeight(i)*luminosityScaling;
        end
    end
end

signalSig = signalCount./sqrt(backgroundCount);

%S/root(B) full sample (not just testing)
radiusSignalFull = [7,14,21,28];
signalCountFull = zeros(1,length(radiusSignalFull));
backgroundCountFull = zeros(1,length(radiusSignalFull));

for k = 1:length(radiusSignalFull)
    for i=1:length(res4b4tag.weight)
        if (res4b4tag.m_H0(i) - optimal_x_HHmH0(3))^2 + (res4b4tag.m_H1(i) - optimal_x_HHmH1(3))^2 < radiusSignalFull(k)^2 %Centered on each fit peak
            backgroundCountFull(k) = backgroundCountFull(k) + res4b4tag.weight(i)*luminosityScaling;
        end
    end
    for i=1:length(resHH4tag.weight)
        if (resHH4tag.m_H0(i) - optimal_x_HHmH0(3))^2 + (resHH4tag.m_H1(i) - optimal_x_HHmH1(3))^2 < radiusSignalFull(k)^2 %Centered on each fit peak
            signalCountFull(k) = signalCountFull(k) + resHH4tag.weight(i)*luminosityScaling;
        end
    end
end

signalSigFull = signalCountFull./sqrt(backgroundCountFull);


figure()
sp1 = subplot(1,2,1);
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',bin4b4res,'FaceColor','flat','DisplayStyle','tile')
title('4b Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
caxis([0,6e4])
sp2 = subplot(1,2,2);
histogram2('XBinEdges',binValue14b4res,'YBinEdges',binValue24b4res,'BinCounts',fit4b4res2D,'FaceColor','flat','DisplayStyle','tile')
title('4b 4-tag Resolved 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
caxis([0,6e4])

set(sp1,'Units','normalized');
set(sp2,'Units','normalized');
set(sp1,'position',[0.075,0.125,0.35,0.8]);
set(sp2,'position',[0.625,0.125,0.35,0.8]);