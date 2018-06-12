function [binEventSum,binEventSumVar,binValue1,binValue2] = HistoBins(weight,variable1,variable2,luminosityScaling,numBins,Plots,name,nameVar1,nameVar2)
%Plots == 1 (make plots), Plots == 0 (don't)

%Creates bins for a bivariate histogram (Var1,Var2) according to numBins.
%Also outputs variance of each bin
%total bins = (numBins)^2
%Sums up weights*luminosityScaling in each bin region

% ~~~
binValue1 = linspace(85,165,numBins); %Max value in each bin
binValue2 = linspace(85,165,numBins); % Setup for 85 to 165 GeV mass distributions
% ~~~ 

binEventSum = zeros(numBins-1,numBins-1);
binEventSumVar = zeros(numBins-1,numBins-1);
binEventSumErrNorm = zeros(numBins-1,numBins-1);

for i=1:length(weight) %This loop saves weight of each sample into correct bin
    for j=2:numBins
    for k=2:numBins
        if variable1(i)<=binValue1(j) && variable1(i)>=binValue1(j-1) && variable2(i)<=binValue2(k) && variable2(i)>=binValue2(k-1)
            binEventSum(j-1,k-1) = binEventSum(j-1,k-1) + weight(i)*luminosityScaling; %Adds on number of events (weight*scaling) to correct bin
            binEventSumVar(j-1,k-1) = binEventSumVar(j-1,k-1) + (weight(i)*luminosityScaling)^2;
            if binEventSum(j-1,k-1) ~= 0, binEventSumErrNorm(j-1,k-1) = sqrt(binEventSumVar(j-1,k-1))/binEventSum(j-1,k-1); end
        end
    end
    end
end

if Plots == 1
    histogramPlots(binEventSum,binValue1,binValue2,name,nameVar1,nameVar2)
    histogramPlotsErr(binEventSumErrNorm,binValue1,binValue2,name,nameVar1,nameVar2)        
end