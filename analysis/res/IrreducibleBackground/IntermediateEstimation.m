% Estimation of irreducible background: for the INTERMEDIATE channel

close all
clear all

load('sampleHighStats')  % Need to load HH 4-tag, and 4b 4-tag

% Global variables
luminosityScaling = 4000; %Integrated Luminosity
numBins = 17; %For the histogramming data
Plots = 1; %To produce the distribution plots of the samples at the beginning

% Channel Sort
[resHH4tag,intHH4tag,boostHH4tag] = channelSort(HH4tag);
[res4b4tag,int4b4tag,boost4b4tag] = channelSort(qcd4b4tag);

%Split into training and Samples
int4b4tagTrainWeight = int4b4tag.weight(1:length(int4b4tag.weight)/2,1); %Training (first half)
int4b4tagTrainmH0 = int4b4tag.m_H0(1:length(int4b4tag.weight)/2,1); %mH0 
int4b4tagTrainmH1 = int4b4tag.m_H1(1:length(int4b4tag.weight)/2,1); %mH1 
int4b4tagSampleWeight = int4b4tag.weight(round(length(int4b4tag.weight)/2):round(length(int4b4tag.weight)),1); %Sample to fit (second half)
int4b4tagSamplemH0 = int4b4tag.m_H0(round(length(int4b4tag.weight)/2):round(length(int4b4tag.weight)),1); %mH0 
int4b4tagSamplemH1 = int4b4tag.m_H1(round(length(int4b4tag.weight)/2):round(length(int4b4tag.weight)),1); %mH1 

intHH4tagTrainWeight = intHH4tag.weight(1:length(intHH4tag.weight)/2,1); %Training (first half)
intHH4tagTrainmH0 = intHH4tag.m_H0(1:length(intHH4tag.weight)/2,1); %mH0 
intHH4tagTrainmH1 = intHH4tag.m_H1(1:length(intHH4tag.weight)/2,1); %mH1
intHH4tagSampleWeight = intHH4tag.weight(round(length(intHH4tag.weight)/2):length(intHH4tag.weight),1); %Sample to fit (second half)
intHH4tagSamplemH0 = intHH4tag.m_H0(round(length(intHH4tag.weight)/2):length(intHH4tag.weight),1); %mH0 
intHH4tagSamplemH1 = intHH4tag.m_H1(round(length(intHH4tag.weight)/2):length(intHH4tag.weight),1); %mH1

% Histograms (mH0, mH1) - 4b and HH (training and testing halves)
[bin4b4int,bin4b4Varint,binValue14b4int,binValue24b4int] = HistoBins(int4b4tagTrainWeight,int4b4tagTrainmH0,int4b4tagTrainmH1,luminosityScaling,numBins,Plots,'4b 4-tag Intermediate Training Sample','m_H0','m_H1');
[binHH4int,binHH4Varint,binValue1HH4int,binValue2HH4int] = HistoBins(intHH4tagTrainWeight,intHH4tagTrainmH0,intHH4tagTrainmH1,luminosityScaling,numBins,Plots,'HH 4-tag Intermediate Training Sample','m_H0','m_H1'); %16 bins
[binHH4int32,binHH4Varint32,binValue1HH4int32,binValue2HH4int32] = HistoBins(intHH4tagTrainWeight,intHH4tagTrainmH0,intHH4tagTrainmH1,luminosityScaling,33,Plots,'HH 4-tag Intermediate Training Sample','m_H0','m_H1'); %32 bins

[bin4b4intSample,bin4b4VarintSample,binValue14b4intSample,binValue24b4intSample] = HistoBins(int4b4tagSampleWeight,int4b4tagSamplemH0,int4b4tagSamplemH1,luminosityScaling,numBins,Plots,'4b 4-tag Intermediate Testing Sample','m_H0','m_H1');
[binHH4intSample,binHH4VarintSample,binValue1HH4intSample,binValue2HH4intSample] = HistoBins(intHH4tagSampleWeight,intHH4tagSamplemH0,intHH4tagSamplemH1,luminosityScaling,numBins,Plots,'HH 4-tag Intermediate Testing Sample','m_H0','m_H1');

% 1D plots - no fit will overlay fits later:
histogram1Dboth(bin4b4int,sqrt(bin4b4Varint),binValue14b4int,binValue24b4int,'4b-4-tag-Intermediate','m_H0','m_H1');
histogram1Dboth(binHH4int32,sqrt(binHH4Varint32),binValue1HH4int32,binValue2HH4int32,'HH-4-tag-Intermediate','m_H0','m_H1');

% 1D fit:
% Preallocate
binValue1mid4b4int = zeros(numBins-1,1);
binValue2mid4b4int = zeros(numBins-1,1);
binValue1midHH4int = zeros(numBins-1,1);
binValue2midHH4int = zeros(numBins-1,1);
binValue1midHH4int32 = zeros(numBins-1,1);
binValue2midHH4int32 = zeros(numBins-1,1);

%Inverse Frac errors for weights in fits
invVar4b = 1./(bin4b4Varint); %4b
invVar4bmH0 = 1./(sum(bin4b4Varint,2)); %4b mH0
invVar4bmH1 = 1./(sum(bin4b4Varint,1)); %4b mH1
invVarHH = 1./(binHH4Varint); % HH
invVarHHmH0 = 1./(sum(binHH4Varint32,2)); %HH mH0
invVarHHmH1 = 1./(sum(binHH4Varint32,1)); %HH mH1
invVarSum = 1./(binHH4VarintSample + bin4b4VarintSample); %Fit is for the sum of testing samples

%Sum Samples
binSumBothSample = bin4b4intSample + binHH4intSample; % For the final fit

%1D fits:
ftLinear = fittype('a*x + b','coefficients',{'a','b'},'independent','x'); %Linear
ftpoly2 = fittype('a + b*x + c.*x.^2','coefficients',{'a','b','c'},'independent','x'); % Quadratic
ftExp = fittype('a*exp(-b*(x-c)) + d','coefficients',{'a','b','c','d'},'independent','x'); % Exponential + constant
ftExp1a = fittype('a*exp(-b*(x))','coefficients',{'a','b'},'independent','x'); %No constant, and no degenerate terms

for i=2:length(binValue14b4int)
    binValue1mid4b4int(i-1) = 0.5*(binValue14b4int(i-1) + binValue14b4int(i));
    binValue2mid4b4int(i-1) = 0.5*(binValue24b4int(i-1) + binValue24b4int(i));
end
%4b 1D fits:
fitmH04b4int = fit(binValue1mid4b4int,sum(bin4b4int,2),ftExp1a,'weight',invVar4bmH0,'startpoint',[4e5,0.03]); %fit Gauss to 4b 4-tag resolved [0.001,125,50,1e-3]
fitmH14b4int = fit(binValue2mid4b4int,sum(bin4b4int,1).',ftpoly2,'weight',invVar4bmH1,'startpoint',[3e5,100,0]); %fit Gauss to 4b 4-tag resolved [1e4,125,10,5e3]  [2e4,125,50,1e3]

for i=2:length(binValue1HH4int)
    binValue1midHH4int(i-1) = 0.5*(binValue1HH4int(i-1) + binValue1HH4int(i));
    binValue2midHH4int(i-1) = 0.5*(binValue2HH4int(i-1) + binValue2HH4int(i));
end
for i=2:length(binValue1HH4int32)
    binValue1midHH4int32(i-1) = 0.5*(binValue1HH4int32(i-1) + binValue1HH4int32(i));
    binValue2midHH4int32(i-1) = 0.5*(binValue2HH4int32(i-1) + binValue2HH4int32(i));
end

%HH 1D fits:
x0_mH0 = [1.5,1.5,125,7,20,2.5,1.5]; %Start Point
optimal_x_HHmH0 = nlinfit(binValue1midHH4int32,sum(binHH4int32,2),@CrystalBall2,x0_mH0,'Weights',invVarHHmH0); %Fit to Crystal Ball function
output_mH0 = CrystalBall2(optimal_x_HHmH0,85:0.1:165); %Best Fit
x0_mH1 = [1,3,125,7,20,1.5,1.5]; %Start Point
optimal_x_HHmH1 = nlinfit(binValue2midHH4int32,sum(binHH4int32,1).',@CrystalBall2,x0_mH1,'Weights',invVarHHmH1.'); %Fit to Crystal Ball function
output_mH1 = CrystalBall2(optimal_x_HHmH1,85:0.1:165); %Best Fit

% Overlaying 1D fits to the 1D histograms (mH0 or mH1) for (4b or HH)
figure(11)
hold on
plot(fitmH14b4int)
hold off
xlabel('m_{H1} / GeV')
ylabel('Number of Events')
title('4b 4-tag Intermediate m_{H1} distribution')
legend('Sample','Error','Quadratic Fit')

figure(12)
hold on
plot(fitmH04b4int)
hold off
xlabel('m_{H0} / GeV')
ylabel('Number of Events')
title('4b 4-tag Intermediate m_{H0} distribution')
legend('Sample','Error','Exponential Fit')

figure(13)
hold on
plot(85:0.1:165,output_mH1)%Plot fit over entire mass range
hold off
xlabel('m_{H1} / GeV')
ylabel('Number of Events')
title('HH 4-tag Intermediate m_{H1} distribution')
legend('Sample','Error','Crystal Ball Fit')

figure(14)
hold on
plot(85:0.1:165,output_mH0) %Plot fit over entire mass range
hold off 
xlabel('m_{H0} / GeV')
ylabel('Number of Events')
title('HH 4-tag Intermediate m_{H0} distribution')
legend('Sample','Error','Crystal Ball Fit')

% Making the 1D fit histograms
bin4b4intFitmH0 = zeros(numBins-1,1);
bin4b4intFitmH1 = zeros(numBins-1,1);

for i=1:(numBins-1)
    bin4b4intFitmH0(i) = fitmH04b4int(binValue1mid4b4int(i));
    bin4b4intFitmH1(i) = fitmH14b4int(binValue2mid4b4int(i));
end

% Histogram of Fit for chi2 later
binHH4intFitmH0 = CrystalBall2(optimal_x_HHmH0,binValue1midHH4int); %16 bins
binHH4intFitmH1 = CrystalBall2(optimal_x_HHmH1,binValue2midHH4int);
binHH4intFitmH032 = CrystalBall2(optimal_x_HHmH0,binValue1midHH4int32); %32 bins
binHH4intFitmH132 = CrystalBall2(optimal_x_HHmH1,binValue2midHH4int32);

%Normalise (For fit later) - Scale each fit to 1 (for distribution function
%instead of differential cross section)
N = sum(sum(binSumBothSample)); %Total number of events in testing sample
% mH0 4b
a1 = fitmH04b4int.a/sum(bin4b4intFitmH0);
b1 = fitmH04b4int.b;
% mH1 4b
a2 = fitmH14b4int.a/sum(bin4b4intFitmH1);
b2 = fitmH14b4int.b/sum(bin4b4intFitmH1);
c2 = fitmH14b4int.c/sum(bin4b4intFitmH1);
% mH0 HH
optimal_x_HHmH0(5) = optimal_x_HHmH0(5)/sum(binHH4intFitmH0); %scaled to PDF
% mH1 HH
optimal_x_HHmH1(5) = optimal_x_HHmH1(5)/sum(binHH4intFitmH1); %So can then use with n bins (ie n=16)

%Pre-allocate
fit4b4int = zeros(numBins-1,numBins-1);
fitHH4int = zeros(numBins-1,numBins-1);

%Combining 1D fits into 2D distribution 
for i=1:(numBins-1)
    for j=1:(numBins-1)
        fit4b4int(i,j) = (a1*exp(-b1*(binValue1mid4b4int(i))))*(a2 + b2*binValue2mid4b4int(j) + c2*(binValue2mid4b4int(j)^2));
        fitHH4int(i,j) = CrystalBall2(optimal_x_HHmH0,binValue1midHH4int(i))*CrystalBall2(optimal_x_HHmH1,binValue2midHH4int(j));
    end
end

% Figure contains: Samples and fitted templates for signal and background
figure()
subplot(2,2,1)
histogram2('XBinEdges',binValue1HH4int,'YBinEdges',binValue2HH4int,'BinCounts',binHH4int,'FaceColor','flat')
title('HH 4-tag Intermediate Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,2)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',bin4b4int,'FaceColor','flat')
title('4b 4-tag Intermediate Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,3)
histogram2('XBinEdges',binValue1HH4int,'YBinEdges',binValue2HH4int,'BinCounts',fitHH4int,'FaceColor','flat')
title('HH fit - power-Gauss-power')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,4)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',fit4b4int,'FaceColor','flat')
title('4b fit - Exp and Linear')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')

% 2D Linear Combination set up:
n=0;
fitValuesBoth=zeros((numBins-1)^2,4);
fitValues4b2D=zeros((numBins-1)^2,4);
for i=1:(numBins-1)
    for j=1:(numBins-1)
        n=n+1;
        fitValuesBoth(n,:) = [binValue1mid4b4int(i),binValue2mid4b4int(j),binSumBothSample(i,j),invVarSum(i,j)];
        fitValues4b2D(n,:) = [binValue1mid4b4int(i),binValue2mid4b4int(j),bin4b4int(i,j),invVar4b(i,j)];
    end
end

 
% Fit the total events using a chi2 fit:
ftSumBoth = fittype( @(alpha,mH0,mH1) N.*((1-alpha).*(a1.*exp(-b1.*(mH0))).*(a2 + b2.*mH1 + c2.*(mH1.^2)) + alpha.*(  optimal_x_HHmH0(5).*(heaviside((optimal_x_HHmH0(3)-optimal_x_HHmH0(1).*optimal_x_HHmH0(4)) - mH0).*(( optimal_x_HHmH0(2) ./ abs(optimal_x_HHmH0(1))).^(optimal_x_HHmH0(2)).*exp(-(abs(optimal_x_HHmH0(1)).^2)/2) .* (optimal_x_HHmH0(2)./abs(optimal_x_HHmH0(1)) - abs(optimal_x_HHmH0(1)) - (mH0 - optimal_x_HHmH0(3))./optimal_x_HHmH0(4)).^-optimal_x_HHmH0(2)) + (heaviside(mH0 - (optimal_x_HHmH0(3)-optimal_x_HHmH0(1).*optimal_x_HHmH0(4))).*heaviside((optimal_x_HHmH0(3)+optimal_x_HHmH0(7).*optimal_x_HHmH0(4)) - mH0).*(exp(-(mH0-optimal_x_HHmH0(3)).^2/(2.*optimal_x_HHmH0(4).^2)))) + (heaviside(mH0 - (optimal_x_HHmH0(3)+optimal_x_HHmH0(7).*optimal_x_HHmH0(4))).*(( optimal_x_HHmH0(6) ./ abs(optimal_x_HHmH0(7))).^(optimal_x_HHmH0(6)).*exp(-(abs(optimal_x_HHmH0(7)).^2)/2) .* (optimal_x_HHmH0(6)./abs(optimal_x_HHmH0(7)) - abs(optimal_x_HHmH0(7)) - (optimal_x_HHmH0(3) - mH0)./optimal_x_HHmH0(4)).^-optimal_x_HHmH0(6)))) ).*(  optimal_x_HHmH1(5).*(heaviside((optimal_x_HHmH1(3)-optimal_x_HHmH1(1).*optimal_x_HHmH1(4)) - mH1).*(( optimal_x_HHmH1(2) ./ abs(optimal_x_HHmH1(1))).^(optimal_x_HHmH1(2)).*exp(-(abs(optimal_x_HHmH1(1)).^2)/2) .* (optimal_x_HHmH1(2)./abs(optimal_x_HHmH1(1)) - abs(optimal_x_HHmH1(1)) - (mH1 - optimal_x_HHmH1(3))./optimal_x_HHmH1(4)).^-optimal_x_HHmH1(2)) + (heaviside(mH1 - (optimal_x_HHmH1(3)-optimal_x_HHmH1(1).*optimal_x_HHmH1(4))).*heaviside((optimal_x_HHmH1(3)+optimal_x_HHmH1(7).*optimal_x_HHmH1(4)) - mH1).*(exp(-(mH1-optimal_x_HHmH1(3)).^2/(2.*optimal_x_HHmH1(4).^2)))) + (heaviside(mH1 - (optimal_x_HHmH1(3)+optimal_x_HHmH1(7).*optimal_x_HHmH1(4))).*(( optimal_x_HHmH1(6) ./ abs(optimal_x_HHmH1(7))).^(optimal_x_HHmH1(6)).*exp(-(abs(optimal_x_HHmH1(7)).^2)/2) .* (optimal_x_HHmH1(6)./abs(optimal_x_HHmH1(7)) - abs(optimal_x_HHmH1(7)) - (optimal_x_HHmH1(3) - mH1)./optimal_x_HHmH1(4)).^-optimal_x_HHmH1(6))))  )),'independent',{'mH0','mH1'});
[fittedBoth,gofitBoth,outputFitBoth] = fit([fitValuesBoth(:,1),fitValuesBoth(:,2)],fitValuesBoth(:,3),ftSumBoth,'weight',fitValuesBoth(:,4),'startpoint',0.01,'upper',1,'lower',0);

% Histogram of chi2 template fit:
binFittedValuesSample = zeros(numBins-1,numBins-1);
for i=1:(numBins-1)
    for j=1:(numBins-1)
        binFittedValuesSample(i,j) = fittedBoth(binValue1mid4b4int(i),binValue2mid4b4int(j));
    end
end

% Figure contains: Testing sample, fit to testing sample with templtates,
% and the two templates used (4b and HH).
figure()
subplot(2,2,1)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',binSumBothSample,'FaceColor','flat')
title('Sum of HH and 4B Intermediate Testing Sample')
colorbar
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,2)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',binFittedValuesSample,'FaceColor','flat')
title('Fitted Distribution of Intermediate Testing Sample')
colorbar
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,3)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',fit4b4int,'FaceColor','flat')
title('Fitted 4b 4-tag Intermediate')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(2,2,4)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',fitHH4int,'FaceColor','flat')
title('Fitted HH 4-tag Intermediate')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')

%Confidence intervals of chi2 fit at 1-4 sigma intervals on fraction of
%signal to background (alpha)
confidence = [confint(fittedBoth,0.6827),confint(fittedBoth,0.9545),confint(fittedBoth,0.9973),confint(fittedBoth,0.99993)];

%Chi2 for each fit:
% Scaling for chi2:
fitHH4intScale = sum(sum(binHH4intSample))*fitHH4int;
fit4b4intScale = sum(sum(bin4b4intSample))*fit4b4int;
fit4b4trainScale = sum(sum(bin4b4int))*fit4b4int; %for chi2 train value

% Residual plots:
% HH, 4b, Sum
% ( ( Sample - Fit )/(SampleErr) )^2
residualsHH2D = (binHH4intSample - fitHH4intScale).^2./(binHH4VarintSample); %chi2
residuals4b2D = (bin4b4intSample - fit4b4intScale).^2./(bin4b4VarintSample);
residualsSum2D = (binSumBothSample - binFittedValuesSample).^2./(binHH4VarintSample + bin4b4VarintSample);
residuals4bmH0 = (sum(bin4b4int,2).' - bin4b4intFitmH0.').^2./(sum(bin4b4Varint,2).');
residuals4bmH1 = (sum(bin4b4int,1) - bin4b4intFitmH1.').^2./(sum(bin4b4Varint,1));
residualsHHmH016 = (sum(binHH4int,2).' - 2*binHH4intFitmH0.').^2./(sum(binHH4Varint,2).');
residualsHHmH116 = (sum(binHH4int,1) - 2*binHH4intFitmH1.').^2./(sum(binHH4Varint,1));
residualsHHmH032 = (sum(binHH4int32,2).' - binHH4intFitmH032.').^2./(sum(binHH4Varint32,2).');
residualsHHmH132 = (sum(binHH4int32,1) - binHH4intFitmH132.').^2./(sum(binHH4Varint32,1));

residuals4btrain2D = (bin4b4int - fit4b4trainScale).^2./(bin4b4Varint);

% Figure contains: each row: (HH, 4b, Testing), and then the Sample, fit
% and residuals for each respectively.
figure()
subplot(3,3,1)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',binHH4intSample,'FaceColor','flat')
title('HH 4-tag Intermediate Testing Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,2)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',fitHH4intScale,'FaceColor','flat')
title('HH 4-tag Intermediate 1D fits')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,3)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',residualsHH2D,'FaceColor','flat')
title('Residuals HH 1D fits')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')

subplot(3,3,4)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',bin4b4intSample,'FaceColor','flat')
title('4b 4-tag Intermediate Testing Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,5)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',fit4b4intScale,'FaceColor','flat')
title('4b 4-tag Intermediate 1D fits')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,6)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',residuals4b2D,'FaceColor','flat')
title('Residuals 4b 1D fits')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')

subplot(3,3,7)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',binSumBothSample,'FaceColor','flat')
title('Signal+Background 4-tag Intermediate Testing Sample')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,8)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',binFittedValuesSample,'FaceColor','flat')
title('Signal+Background 4-tag Intermediate 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Number of Events')
subplot(3,3,9)
histogram2('XBinEdges',binValue14b4int,'YBinEdges',binValue24b4int,'BinCounts',residualsSum2D,'FaceColor','flat')
title('Residuals Signal+Background 2D fit')
xlabel('m_{H0} / GeV')
ylabel('m_{H1} / GeV')
zlabel('Residuals')

% Chi2: 2D fit values:
chi24b = sum(sum(residuals4b2D));
chi2HH = sum(sum(residualsHH2D));
chi2Sum = sum(sum(residualsSum2D));

chi24btrain = sum(sum(residuals4btrain2D));

fprintf('\nManual Chi2 Values (2D):\n 4b: %s\n HH: %s\n Sum: %s\n',chi24b,chi2HH,chi2Sum)

% Chi2: 1D fit values:
chi24bmH0 = sum(sum(residuals4bmH0));
chi24bmH1 = sum(sum(residuals4bmH1));
chi2HHmH016 = sum(sum(residualsHHmH016));
chi2HHmH116 = sum(sum(residualsHHmH116));
chi2HHmH032 = sum(sum(residualsHHmH032));
chi2HHmH132 = sum(sum(residualsHHmH132));

fprintf('\nManual Chi2 Values (1D):\n 4b mH0 (16 bins): %s\n 4b mH1 (16 bins): %s\n HH mH0 (32 bins): %s\n HH mH1 (32 bins): %s\n HH mH0 (16 bins): %s\n HH mH1 (16 bins): %s\n',chi24bmH0,chi24bmH1,chi2HHmH032,chi2HHmH132,chi2HHmH016,chi2HHmH116)

% Manually calculating chi2 for different alpha (fraction of HH) in a
% minimising chi2 fit:
alpha = 0;
bin4b4intNorm = bin4b4int/sum(sum(bin4b4int));
binHH4intNorm = binHH4int/sum(sum(binHH4int));
for i=1:10000
    alpha = alpha + 0.00001;
    for j=1:(numBins-1)
        for k=1:(numBins-1)
            manualBinsSimple(j,k) = N.*( (1-alpha)*fit4b4int(j,k) + alpha*fitHH4int(j,k) );
            alphaZeroCalc(j,k) = N.*( (1-alpha)*bin4b4intNorm(j,k) + alpha*binHH4intNorm(j,k) );
        end
    end
    residualsmanual2DSimple = (binSumBothSample - manualBinsSimple).^2./(binHH4VarintSample + bin4b4VarintSample);
    residualsAlphaZero = (binSumBothSample - alphaZeroCalc).^2./(binHH4VarintSample + bin4b4VarintSample);
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
title('Chi2 v. Alpha - Alpha Zero')
xlabel('alpha - fraction diHiggs')
ylabel('chi^{2}')

%Cross-sections:
fprintf('%s\n %s\n %s\n %s\n %s\n %s\n',sum(sum(int4b4tagTrainWeight)),sum(sum(intHH4tagTrainWeight)),sum(sum(bin4b4intFitmH0))/4000,sum(sum(bin4b4intFitmH1))/4000,sum(sum(binHH4intFitmH032))/4000,sum(sum(binHH4intFitmH132))/4000)
xsec_HH_test = sum(sum(intHH4tagSampleWeight));
xsec_4b_test = sum(sum(int4b4tagSampleWeight));

%S/root(B)
radiusSignal = [7,14,21,28];
signalCount = zeros(1,length(radiusSignal));
backgroundCount = zeros(1,length(radiusSignal));

for k = 1:length(radiusSignal)
    for i=1:length(int4b4tagSampleWeight)
        if (int4b4tagSamplemH0(i) - optimal_x_HHmH0(3))^2 + (int4b4tagSamplemH1(i) - optimal_x_HHmH1(3))^2 < radiusSignal(k)^2
            backgroundCount(k) = backgroundCount(k) + int4b4tagSampleWeight(i)*luminosityScaling;
        end
    end
    for i=1:length(intHH4tagSampleWeight)
        if (intHH4tagSamplemH0(i) - optimal_x_HHmH0(3))^2 + (intHH4tagSamplemH1(i) - optimal_x_HHmH1(3))^2 < radiusSignal(k)^2
            signalCount(k) = signalCount(k) + intHH4tagSampleWeight(i)*luminosityScaling;
        end
    end
end

signalSig = signalCount./sqrt(backgroundCount);

%S/root(B) full sample (not just testing)
radiusSignalFull = [7,14,21,28];
signalCountFull = zeros(1,length(radiusSignalFull));
backgroundCountFull = zeros(1,length(radiusSignalFull));

for k = 1:length(radiusSignalFull)
    for i=1:length(int4b4tag.weight)
        if (int4b4tag.m_H0(i) - optimal_x_HHmH0(3))^2 + (int4b4tag.m_H1(i) - optimal_x_HHmH1(3))^2 < radiusSignalFull(k)^2
            backgroundCountFull(k) = backgroundCountFull(k) + int4b4tag.weight(i)*luminosityScaling;
        end
    end
    for i=1:length(intHH4tag.weight)
        if (intHH4tag.m_H0(i) - optimal_x_HHmH0(3))^2 + (intHH4tag.m_H1(i) - optimal_x_HHmH1(3))^2 < radiusSignalFull(k)^2
            signalCountFull(k) = signalCountFull(k) + intHH4tag.weight(i)*luminosityScaling;
        end
    end
end

signalSigFull = signalCountFull./sqrt(backgroundCountFull);