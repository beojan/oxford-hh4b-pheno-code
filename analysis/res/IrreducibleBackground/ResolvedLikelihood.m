% Likelihood test for resolved channel

%From the training sample we get the distributions:
% HH: Multiplication of the CB functions: CrystalBall2(optimal_x_HHmH0,mH0)*CrystalBall2(optimal_x_HHmH1,mH1)
% 4b: The 2D function: fitGauss2D4b4res to be normalised

fitGauss2D4b4resSUM = sum(sum(fit4b4res2D));
fitGauss2D4b4res.a = fitGauss2D4b4res.a/fitGauss2D4b4resSUM; %Divide the distribution by cross-section
fitGauss2D4b4res.g = fitGauss2D4b4res.g/fitGauss2D4b4resSUM;

weights = [res4b4tagSampleWeight;resHH4tagSampleWeight]; %testing weights
mH0 = [res4b4tagSamplemH0;resHH4tagSamplemH0]; %testing mH0
mH1 = [res4b4tagSamplemH1;resHH4tagSamplemH1]; %testing mH1

%To save running a loop for the sum, calculating fHH and f4b here:
f1_HH = CrystalBall2(optimal_x_HHmH0,mH0).*CrystalBall2(optimal_x_HHmH1,mH1);
f2_4b = fitGauss2D4b4res(mH0,mH1);

alphaAll = -0.1:0.0001:0.1; %Range of fraction Signal to BG to calculate for
LoglikeliHood = zeros(1,length(alphaAll));
for i=1:length(alphaAll)
    alpha = alphaAll(1,i);
    logProb = weights.*log(alpha.*f1_HH + (1-alpha).*f2_4b);
    LoglikeliHood(1,i) = -sum(logProb);
end
figure()
plot(alphaAll,LoglikeliHood)
title('LogL v. Alpha')
xlabel('alpha')
ylabel('Log(Likelihood)')

[minLogL,indLogL] = min(LoglikeliHood);
AlphaMin = alphaAll(indLogL);

minLogL1sigma = minLogL+0.5;
for i=2:length(LoglikeliHood)
    if LoglikeliHood(i-1) <= minLogL1sigma && LoglikeliHood(i) >= minLogL1sigma, fprintf('1 Sigma Lower bound is %.4f\n',alphaAll(i)); end
    if LoglikeliHood(i-1) >= minLogL1sigma && LoglikeliHood(i) <= minLogL1sigma, fprintf('1 Sigma Upper bound is %.4f\n',alphaAll(i)); end
end