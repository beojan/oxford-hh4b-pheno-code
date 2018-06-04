% Estimation of irreducible background: for the BOOSTED channel

close all
clear all

load('sampleHighStats')  % Need to load HH 4-tag, and 4b 4-tag

% Global variables
luminosityScaling = 4000; %Integrated Luminosity
numBins = 17; % For the histograms
Plots = 1; %To produce the distribution plots of the samples at the beginning

% Channel Sort
[resHH4tag,intHH4tag,boostHH4tag] = channelSort(HH4tag);
[res4b4tag,int4b4tag,boost4b4tag] = channelSort(qcd4b4tag);

boost4b4tag(39473,:) = [];

%Split into training and Samples
boost4b4tagTrainWeight = boost4b4tag.weight(1:length(boost4b4tag.weight)/2,1); %Training (first half)
boost4b4tagTrainmH0 = boost4b4tag.m_H0(1:length(boost4b4tag.weight)/2,1); %mH0
boost4b4tagTrainmH1 = boost4b4tag.m_H1(1:length(boost4b4tag.weight)/2,1); %mH1
boost4b4tagSampleWeight = boost4b4tag.weight(round(length(boost4b4tag.weight)/2):round(length(boost4b4tag.weight)),1); %Sample to fit (second half)
boost4b4tagSamplemH0 = boost4b4tag.m_H0(round(length(boost4b4tag.weight)/2):round(length(boost4b4tag.weight)),1); %mH0
boost4b4tagSamplemH1 = boost4b4tag.m_H1(round(length(boost4b4tag.weight)/2):round(length(boost4b4tag.weight)),1); %mH1

boostHH4tagTrainWeight = boostHH4tag.weight(1:length(boostHH4tag.weight)/2,1); %Training (first half)
boostHH4tagTrainmH0 = boostHH4tag.m_H0(1:length(boostHH4tag.weight)/2,1); %mH0
boostHH4tagTrainmH1 = boostHH4tag.m_H1(1:length(boostHH4tag.weight)/2,1); %mH1
boostHH4tagSampleWeight = boostHH4tag.weight(round(length(boostHH4tag.weight)/2):round(length(boostHH4tag.weight)),1); %Sample to fit (second half)
boostHH4tagSamplemH0 = boostHH4tag.m_H0(round(length(boostHH4tag.weight)/2):round(length(boostHH4tag.weight)),1); %mH0
boostHH4tagSamplemH1 = boostHH4tag.m_H1(round(length(boostHH4tag.weight)/2):round(length(boostHH4tag.weight)),1); %mH1

% Histograms (mH0, mH1) (4b and HH, training and testing halves).
[bin4b4boost,bin4b4Varboost,binValue14b4boost,binValue24b4boost] = HistoBins(boost4b4tagTrainWeight,boost4b4tagTrainmH0,boost4b4tagTrainmH1,luminosityScaling,numBins,Plots,'4b 4-tag Boosted Training Sample','m_H0','m_H1');
[binHH4boost,binHH4Varboost,binValue1HH4boost,binValue2HH4boost] = HistoBins(boostHH4tagTrainWeight,boostHH4tagTrainmH0,boostHH4tagTrainmH1,luminosityScaling,numBins,Plots,'HH 4-tag Boosted Training Sample','m_H0','m_H1');
[binHH4boost32,binHH4Varboost32,binValue1HH4boost32,binValue2HH4boost32] = HistoBins(boostHH4tagTrainWeight,boostHH4tagTrainmH0,boostHH4tagTrainmH1,luminosityScaling,33,Plots,'HH 4-tag Boosted Training Sample','m_H0','m_H1');

[bin4b4boostSample,bin4b4VarboostSample,binValue14b4boostSample,binValue24b4boostSample] = HistoBins(boost4b4tagSampleWeight,boost4b4tagSamplemH0,boost4b4tagSamplemH1,luminosityScaling,numBins,Plots,'4b 4-tag Boosted Testing Sample','m_H0','m_H1');
[binHH4boostSample,binHH4VarboostSample,binValue1HH4boostSample,binValue2HH4boostSample] = HistoBins(boostHH4tagSampleWeight,boostHH4tagSamplemH0,boostHH4tagSamplemH1,luminosityScaling,numBins,Plots,'HH 4-tag Boosted Testing Sample','m_H0','m_H1');

% 1D plots - no fit
histogram1Dboth(bin4b4boost,sqrt(bin4b4Varboost),binValue14b4boost,binValue24b4boost,'4b-4-tag-Boosted','m_H0','m_H1');
histogram1Dboth(binHH4boost32,sqrt(binHH4Varboost32),binValue1HH4boost32,binValue2HH4boost32,'HH-4-tag-Boosted','m_H0','m_H1');

% 1D fit:
% Preallocate
binValue1mid4b4boost = zeros(numBins-1,1);
binValue2mid4b4boost = zeros(numBins-1,1);
binValue1midHH4boost = zeros(numBins-1,1);
binValue2midHH4boost = zeros(numBins-1,1);
binValue1midHH4boost32 = zeros(numBins-1,1);
binValue2midHH4boost32 = zeros(numBins-1,1);

%Inverse Frac errors for weights in fits
invVar4b = 1./(bin4b4Varboost); %4b
invVar4bmH0 = 1./(sum(bin4b4Varboost,2)); %4b mH0
invVar4bmH1 = 1./(sum(bin4b4Varboost,1)); %4b mH1
invVarHH = 1./(binHH4Varboost); %HH
invVarHHmH0 = 1./(sum(binHH4Varboost32,2)); %HH mH0
invVarHHmH1 = 1./(sum(binHH4Varboost32,1)); %HH mH1
invVarSum = 1./(binHH4VarboostSample + bin4b4VarboostSample); %Fit is for the sum of testing samples

% Normalise to 1
binSumBothSample = bin4b4boostSample + binHH4boostSample; % For the final fit

%1D fits:
ftExp = fittype('a*exp(-b*(x)) + c','coefficients',{'a','b','c'},'independent','x'); %Exponential fit, not degenerate
ftExp1a = fittype('a*exp(-b*(x))','coefficients',{'a','b'},'independent','x'); %No constant, exponential
ftinvLogistic = fittype('a/(1+b*exp(c*(x)))','coefficients',{'a','b','c'},'independent','x'); %4 coeffs, Logistic
% Making histogram bins:
for i=2:length(binValue14b4boost)
    binValue1mid4b4boost(i-1) = 0.5*(binValue14b4boost(i-1) + binValue14b4boost(i));
    binValue2mid4b4boost(i-1) = 0.5*(binValue24b4boost(i-1) + binValue24b4boost(i));
end
for i=2:length(binValue1HH4boost)
    binValue1midHH4boost(i-1) = 0.5*(binValue1HH4boost(i-1) + binValue1HH4boost(i));
    binValue2midHH4boost(i-1) = 0.5*(binValue2HH4boost(i-1) + binValue2HH4boost(i));
end
for i=2:length(binValue1HH4boost32)
    binValue1midHH4boost32(i-1) = 0.5*(binValue1HH4boost32(i-1) + binValue1HH4boost32(i));
    binValue2midHH4boost32(i-1) = 0.5*(binValue2HH4boost32(i-1) + binValue2HH4boost32(i));
end

%4b 1D fits:
fitmH04b4boost = fit(binValue1mid4b4boost,sum(bin4b4boost,2),ftinvLogistic,'weight',invVar4bmH0,'startpoint',[1e5,0.075,0.05]); %Gauss: [0.01,80,25,1e-3] |Exp: [1e4,1,80,1e3] |Linear: []
fitmH14b4boost = fit(binValue2mid4b4boost,sum(bin4b4boost,1).',ftExp1a,'weight',invVar4bmH1,'startpoint',[2e5,0.03]); %Gauss: [1e5,80,25,1e3] |Exp: [1e4,1,80,1e3] |Linear: []

%HH 1D fits:
x0_mH0 = [1.5,1.5,125,7,20,2.5,1.5]; %Start Point
optimal_x_HHmH0 = nlinfit(binValue1midHH4boost32,sum(binHH4boost32,2),@CrystalBall2,x0_mH0,'Weights',invVarHHmH0); %Fit to Crystal Ball function
output_mH0 = CrystalBall2(optimal_x_HHmH0,85:0.1:165); %Best Fit
x0_mH1 = [1,3,125,7,20,1.5,1.5]; %Start Point
optimal_x_HHmH1 = nlinfit(binValue2midHH4boost32,sum(binHH4boost32,1).',@CrystalBall2,x0_mH1,'Weights',invVarHHmH1.');  %Fit to Crystal Ball function
output_mH1 = CrystalBall2(optimal_x_HHmH1,85:0.1:165); %Best Fit

% Overlaying 1D fits to the 1D histograms (mH0 or mH1) for (4b or HH)
figure(11)
hold on
plot(fitmH14b4boost)
hold off
xlabel('m_{H1} / GeV')
ylabel('Number of Events')
title('4b 4-tag Boosted m_{H1} distribution')
legend('Sample','Error','Exponential Fit')

figure(12)
hold on
plot(fitmH04b4boost)
hold off
xlabel('m_{H0} / GeV')
ylabel('Number of Events')
title('4b 4-tag Boosted m_{H0} distribution')
legend('Sample','Error','Logistic Curve Fit')

figure(13)
hold on
plot(85:0.1:165,output_mH1)
hold off
xlabel('m_{H1} / GeV')
ylabel('Number of Events')
title('HH 4-tag Boosted m_{H1} distribution')
legend('Sample','Error','Crystal Ball Fit')

figure(14)
hold on
plot(85:0.1:165,output_mH0)
hold off
xlabel('m_{H0} / GeV')
ylabel('Number of Events')
title('HH 4-tag Boosted m_{H0} distribution')
legend('Sample','Error','Crystal Ball Fit')

%Making the 1D fit histograms
fitHH4boost=zeros(numBins-1,numBins-1);
fit4b4boost=zeros(numBins-1,numBins-1);
bin4b4boostFitmH0 = zeros(numBins-1,1);
bin4b4boostFitmH1 = zeros(numBins-1,1);

for i=1:(numBins-1)
    bin4b4boostFitmH0(i) = fitmH04b4boost(binValue1mid4b4boost(i));
    bin4b4boostFitmH1(i) = fitmH14b4boost(binValue2mid4b4boost(i));
end

%Histogram of fits for chi2:
binHH4boostFitmH0 = CrystalBall2(optimal_x_HHmH0,binValue1midHH4boost); %16 bins
binHH4boostFitmH1 = CrystalBall2(optimal_x_HHmH1,binValue2midHH4boost); 
binHH4boostFitmH032 = CrystalBall2(optimal_x_HHmH0,binValue1midHH4boost32); %32 bins
binHH4boostFitmH132 = CrystalBall2(optimal_x_HHmH1,binValue2midHH4boost32); 

%Normalise (For fit later) - Scale each fit to 1 (for distribution function
%instead of differential cross section)
N = sum(sum(binSumBothSample));
% mH0 4b
a = fitmH04b4boost.a/sum(bin4b4boostFitmH0);
b = fitmH04b4boost.b;
c = fitmH04b4boost.c;
% mH1 4b
a2 = fitmH14b4boost.a/sum(bin4b4boostFitmH1);
b2 = fitmH14b4boost.b;
% mH0 HH
optimal_x_HHmH0(5) = optimal_x_HHmH0(5)/sum(binHH4boostFitmH0); %scaled to PDF
% mH1 HH
optimal_x_HHmH1(5) = optimal_x_HHmH1(5)/sum(binHH4boostFitmH1); %So can then use with n bins (ie n=16)


%Combining 1D fits into 2D distribution (template)
for i=1:(numBins-1)
    for j=1:(numBins-1)
        fit4b4boost(i,j) = (a/(1+b*exp(c*(binValue1mid4b4boost(i)))))*(a2*exp(-b2*(binValue2mid4b4boost(j))) );% + c2);
        fitHH4boost(i,j) = CrystalBall2(optimal_x_HHmH0,binValue1midHH4boost(i))*CrystalBall2(optimal_x_HHmH1,binValue2midHH4boost(j));
    end
end

% Figure contains: Samples and fitted templates for signal and background
figure()
subplot(2,2,1)
histogram2('XBinEdges',binValue1HH4boost,'YBinEdges',binValue2HH4boost,'BinCounts',binHH4boost,'FaceColor','flat')
title('HH 4-tag Boosted Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,2)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',bin4b4boost,'FaceColor','flat')
title('4b 4-tag Boosted Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
colorbar
caxis([0,1500])
subplot(2,2,3)
histogram2('XBinEdges',binValue1HH4boost,'YBinEdges',binValue2HH4boost,'BinCounts',fitHH4boost,'FaceColor','flat')
title('HH fit - 1D gaussians')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,4)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',fit4b4boost,'FaceColor','flat')
title('4b fit - 1D gaussians')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')

% 2D Linear Combination set up:
n=0;
fitValuesBoth=zeros((numBins-1)^2,4);
for i=1:(numBins-1)
    for j=1:(numBins-1)
        n=n+1;
        fitValuesBoth(n,:) = [binValue1mid4b4boost(i),binValue2mid4b4boost(j),binSumBothSample(i,j),invVarSum(i,j)];
    end
end

% Fit the total events using a chi2 fit:
ftSumBoth = fittype( @(alpha,mH0,mH1) N.*((1-alpha).*(a./(1+b.*exp(c.*(mH0)))).*(a2*exp(-b2*(mH1)) ) + alpha.*(  optimal_x_HHmH0(5).*(heaviside((optimal_x_HHmH0(3)-optimal_x_HHmH0(1).*optimal_x_HHmH0(4)) - mH0).*(( optimal_x_HHmH0(2) ./ abs(optimal_x_HHmH0(1))).^(optimal_x_HHmH0(2)).*exp(-(abs(optimal_x_HHmH0(1)).^2)/2) .* (optimal_x_HHmH0(2)./abs(optimal_x_HHmH0(1)) - abs(optimal_x_HHmH0(1)) - (mH0 - optimal_x_HHmH0(3))./optimal_x_HHmH0(4)).^-optimal_x_HHmH0(2)) + (heaviside(mH0 - (optimal_x_HHmH0(3)-optimal_x_HHmH0(1).*optimal_x_HHmH0(4))).*heaviside((optimal_x_HHmH0(3)+optimal_x_HHmH0(7).*optimal_x_HHmH0(4)) - mH0).*(exp(-(mH0-optimal_x_HHmH0(3)).^2/(2.*optimal_x_HHmH0(4).^2)))) + (heaviside(mH0 - (optimal_x_HHmH0(3)+optimal_x_HHmH0(7).*optimal_x_HHmH0(4))).*(( optimal_x_HHmH0(6) ./ abs(optimal_x_HHmH0(7))).^(optimal_x_HHmH0(6)).*exp(-(abs(optimal_x_HHmH0(7)).^2)/2) .* (optimal_x_HHmH0(6)./abs(optimal_x_HHmH0(7)) - abs(optimal_x_HHmH0(7)) - (optimal_x_HHmH0(3) - mH0)./optimal_x_HHmH0(4)).^-optimal_x_HHmH0(6)))) ).*(  optimal_x_HHmH1(5).*(heaviside((optimal_x_HHmH1(3)-optimal_x_HHmH1(1).*optimal_x_HHmH1(4)) - mH1).*(( optimal_x_HHmH1(2) ./ abs(optimal_x_HHmH1(1))).^(optimal_x_HHmH1(2)).*exp(-(abs(optimal_x_HHmH1(1)).^2)/2) .* (optimal_x_HHmH1(2)./abs(optimal_x_HHmH1(1)) - abs(optimal_x_HHmH1(1)) - (mH1 - optimal_x_HHmH1(3))./optimal_x_HHmH1(4)).^-optimal_x_HHmH1(2)) + (heaviside(mH1 - (optimal_x_HHmH1(3)-optimal_x_HHmH1(1).*optimal_x_HHmH1(4))).*heaviside((optimal_x_HHmH1(3)+optimal_x_HHmH1(7).*optimal_x_HHmH1(4)) - mH1).*(exp(-(mH1-optimal_x_HHmH1(3)).^2/(2.*optimal_x_HHmH1(4).^2)))) + (heaviside(mH1 - (optimal_x_HHmH1(3)+optimal_x_HHmH1(7).*optimal_x_HHmH1(4))).*(( optimal_x_HHmH1(6) ./ abs(optimal_x_HHmH1(7))).^(optimal_x_HHmH1(6)).*exp(-(abs(optimal_x_HHmH1(7)).^2)/2) .* (optimal_x_HHmH1(6)./abs(optimal_x_HHmH1(7)) - abs(optimal_x_HHmH1(7)) - (optimal_x_HHmH1(3) - mH1)./optimal_x_HHmH1(4)).^-optimal_x_HHmH1(6)))))),'independent',{'mH0','mH1'});
[fittedBoth,gofitBoth,outputFitBoth] = fit([fitValuesBoth(:,1),fitValuesBoth(:,2)],fitValuesBoth(:,3),ftSumBoth,'weight',fitValuesBoth(:,4),'startpoint',0.001,'upper',1,'lower',0);

% Histogram of chi2 template fit:
binFittedValuesSample = zeros(numBins-1,numBins-1);
for i=1:(numBins-1)
    for j=1:(numBins-1)
        binFittedValuesSample(i,j) = fittedBoth(binValue1mid4b4boost(i),binValue2mid4b4boost(j));
    end
end

% Figure contains: Testing sample, fit to testing sample with templtates,
% and the two templates used (4b and HH).
figure()
subplot(2,2,1)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',binSumBothSample,'FaceColor','flat')
title('Sum of HH and 4B Boosted Testing Sample')
colorbar
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,2)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',binFittedValuesSample,'FaceColor','flat')
title('Fitted Distribution of Boosted Testing Sample')
colorbar
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,3)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',fit4b4boost,'FaceColor','flat')
title('Fitted 4b 4-tag Boosted')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,4)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',fitHH4boost,'FaceColor','flat')
title('Fitted HH 4-tag Boosted')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')

%Confidence intervals of chi2 fit at 1-4 sigma intervals on fraction of
%signal to background (alpha)
confidence = [confint(fittedBoth,0.6827),confint(fittedBoth,0.9545),confint(fittedBoth,0.9973),confint(fittedBoth,0.99993)];

%Chi2 for each fit:
% Scaling for chi2:
fitHH4boostScale = sum(sum(binHH4boostSample))*fitHH4boost;
fit4b4boostScale = sum(sum(bin4b4boostSample))*fit4b4boost;

% Residual plots:
% HH, 4b, Sum
% ( ( Sample - Fit )/(SampleErr) )^2
residualsHH2D = (binHH4boostSample - fitHH4boostScale).^2./(binHH4VarboostSample); %chi2
residuals4b2D = (bin4b4boostSample - fit4b4boostScale).^2./(bin4b4VarboostSample);
residualsSum2D = (binSumBothSample - binFittedValuesSample).^2./(binHH4VarboostSample + bin4b4VarboostSample);
residuals4bmH0 = (sum(bin4b4boost,2).' - bin4b4boostFitmH0).^2./(sum(bin4b4Varboost,2).'); %1D fits
residuals4bmH1 = (sum(bin4b4boost,1) - bin4b4boostFitmH1).^2./(sum(bin4b4Varboost,1));
residualsHHmH016 = (sum(binHH4boost,2).' - 2*binHH4boostFitmH0.').^2./(sum(binHH4Varboost,2).');
residualsHHmH116 = (sum(binHH4boost,1) - 2*binHH4boostFitmH1.').^2./(sum(binHH4Varboost,1));
residualsHHmH032 = (sum(binHH4boost32,2).' - binHH4boostFitmH032.').^2./(sum(binHH4Varboost32,2).');
residualsHHmH132 = (sum(binHH4boost32,1) - binHH4boostFitmH132.').^2./(sum(binHH4Varboost32,1));


% Figure contains: each row: (HH, 4b, Testing), and then the Sample, fit
% and residuals for each respectively.
figure()
subplot(3,3,1)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',binHH4boostSample,'FaceColor','flat')
title('HH 4-tag Boosted Testing Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,2)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',fitHH4boostScale,'FaceColor','flat')
title('HH 4-tag Boosted 1D fits')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,3)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',residualsHH2D,'FaceColor','flat')
title('Residuals HH 1D fits')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')

subplot(3,3,4)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',bin4b4boostSample,'FaceColor','flat')
title('4b 4-tag Boosted Testing Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,5)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',fit4b4boostScale,'FaceColor','flat')
title('4b 4-tag Boosted 1D fits')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,6)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',residuals4b2D,'FaceColor','flat')
title('Residuals 4b 1D fits')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')

subplot(3,3,7)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',binSumBothSample,'FaceColor','flat')
title('Signal+Background 4-tag Boosted Testing Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,8)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',binFittedValuesSample,'FaceColor','flat')
title('Signal+Background 4-tag Boosted 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,9)
histogram2('XBinEdges',binValue14b4boost,'YBinEdges',binValue24b4boost,'BinCounts',residualsSum2D,'FaceColor','flat')
title('Residuals Signal+Background 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')


% Chi2: 2D fit values:
chi24b = sum(sum(residuals4b2D));
chi2HH = sum(sum(residualsHH2D));
chi2Sum = sum(sum(residualsSum2D));

fprintf('\nManual Chi2 Values (2D):\n 4b: %s\n HH: %s\n Sum: %s\n',chi24b,chi2HH,chi2Sum)

% Chi2: 1D fit values:
chi24bmH0 = sum(sum(residuals4bmH0));
chi24bmH1 = sum(sum(residuals4bmH1));
chi2HHmH016 = sum(sum(residualsHHmH016));
chi2HHmH116 = sum(sum(residualsHHmH116));
chi2HHmH032 = sum(sum(residualsHHmH032));
chi2HHmH132 = sum(sum(residualsHHmH132));

fprintf('\nManual Chi2 Values (1D):\n 4b mH0 (16 bins): %s\n 4b mH1 (16 bins): %s\n HH mH0 (32 bins): %s\n HH mH1 (32 bins): %s\n HH mH0 (16 bins): %s\n HH mH1 (16 bins): %s\n',chi24bmH0,chi24bmH1,chi2HHmH032,chi2HHmH132,chi2HHmH016,chi2HHmH116)

testRatio = binHH4boostSample./bin4b4boostSample;


% Manually calculating chi2 for different alpha (fraction of HH) in a
% minimising chi2 fit:
alpha = 0;
bin4b4boostNorm = bin4b4boost/sum(sum(bin4b4boost));
binHH4boostNorm = binHH4boost/sum(sum(binHH4boost));
for i=1:10000
    alpha = alpha + 0.00001;
    for j=1:(numBins-1)
        for k=1:(numBins-1)
            manualBinsSimple(j,k) = N.*( (1-alpha)*fit4b4boost(j,k) + alpha*fitHH4boost(j,k));
            AlphaZeroCalc(j,k) = N.*( (1-alpha)*bin4b4boostNorm(j,k) + alpha*binHH4boostNorm(j,k));
        end    
    end
    residualsmanual2DSimple = (binSumBothSample - manualBinsSimple).^2./(binHH4VarboostSample + bin4b4VarboostSample);
    residualsAlphaZero = (binSumBothSample - AlphaZeroCalc).^2./(binHH4VarboostSample + bin4b4VarboostSample);
    chi2ManualSimple(i,1) = sum(sum(residualsmanual2DSimple));
    chi2AlphaZero(i,1) = sum(sum(residualsAlphaZero));
    alphaManual(i,1) = alpha;
end

figure()
plot(alphaManual,chi2ManualSimple)
title('Simple calculation')
xlabel('alpha - fraction diHiggs')
ylabel('chi^{2}')
figure()
plot(alphaManual,chi2AlphaZero)
title('Alpha Zero')
xlabel('alpha - fraction diHiggs')
ylabel('chi^{2}')

%Cross-sections:
fprintf('%s\n %s\n %s\n %s\n %s\n %s\n',sum(sum(boost4b4tagTrainWeight)),sum(sum(boostHH4tagTrainWeight)),sum(sum(bin4b4boostFitmH0))/4000,sum(sum(bin4b4boostFitmH1))/4000,sum(sum(binHH4boostFitmH032))/4000,sum(sum(binHH4boostFitmH132))/4000)
xsec_HH_test = sum(sum(boostHH4tagSampleWeight));
xsec_4b_test = sum(sum(boost4b4tagSampleWeight));


%S/root(B)
radiusSignal = [7,14,21,28];
signalCount = zeros(1,length(radiusSignal));
backgroundCount = zeros(1,length(radiusSignal));

for k = 1:length(radiusSignal)
    for i=1:length(boost4b4tagSampleWeight)
        if (boost4b4tagSamplemH0(i) - optimal_x_HHmH0(3))^2 + (boost4b4tagSamplemH1(i) - optimal_x_HHmH1(3))^2 < radiusSignal(k)^2 %Centered on each fit peak
            backgroundCount(k) = backgroundCount(k) + boost4b4tagSampleWeight(i)*luminosityScaling;
        end
    end
    for i=1:length(boostHH4tagSampleWeight)
        if (boostHH4tagSamplemH0(i) - optimal_x_HHmH0(3))^2 + (boostHH4tagSamplemH1(i) - optimal_x_HHmH1(3))^2 < radiusSignal(k)^2 %Centered on each fit peak
            signalCount(k) = signalCount(k) + boostHH4tagSampleWeight(i)*luminosityScaling;
        end
    end
end

signalSig = signalCount./sqrt(backgroundCount);

%S/root(B) full sample (not just testing)
radiusSignalFull = [7,14,21,28];
signalCountFull = zeros(1,length(radiusSignalFull));
backgroundCountFull = zeros(1,length(radiusSignalFull));

for k = 1:length(radiusSignalFull)
    for i=1:length(boost4b4tag.weight)
        if (boost4b4tag.m_H0(i) - optimal_x_HHmH0(3))^2 + (boost4b4tag.m_H1(i) - optimal_x_HHmH1(3))^2 < radiusSignalFull(k)^2 %Centered on each fit peak
            backgroundCountFull(k) = backgroundCountFull(k) + boost4b4tag.weight(i)*luminosityScaling;
        end
    end
    for i=1:length(boostHH4tag.weight)
        if (boostHH4tag.m_H0(i) - optimal_x_HHmH0(3))^2 + (boostHH4tag.m_H1(i) - optimal_x_HHmH1(3))^2 < radiusSignalFull(k)^2 %Centered on each fit peak
            signalCountFull(k) = signalCountFull(k) + boostHH4tag.weight(i)*luminosityScaling;
        end
    end
end

signalSigFull = signalCountFull./sqrt(backgroundCountFull);