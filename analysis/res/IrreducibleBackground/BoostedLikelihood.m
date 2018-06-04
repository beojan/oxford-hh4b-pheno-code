% Likelihood test for BOOSTED channel

%From the training sample we get the distributions:
% HH: Multiplication of the CB functions: CrystalBall2(optimal_x_HHmH0,mH0)*CrystalBall2(optimal_x_HHmH1,mH1)
% 4b: Multiplcation of the normalised 1D fits:
% fitmH04b4int(mH0)*fitmH14b4int(mH1)

fitmH04b4boostSUM = sum(sum(bin4b4boostFitmH0)); %invLogistic, only normalise 'a'
fitmH04b4boost.a = fitmH04b4boost.a/fitmH04b4boostSUM;
fitmH14b4boostSUM = sum(sum(bin4b4boostFitmH1)); %exp, normalise 'a' and 'd'
fitmH14b4boost.a = fitmH14b4boost.a/fitmH14b4boostSUM;


weights = [boost4b4tagSampleWeight;boostHH4tagSampleWeight]; %testing weights
mH0 = [boost4b4tagSamplemH0;boostHH4tagSamplemH0]; %testing mH0
mH1 = [boost4b4tagSamplemH1;boostHH4tagSamplemH1]; %testing mH1

%To save running a loop for the sum, calculating fHH and f4b here:
f1_HH = CrystalBall2(optimal_x_HHmH0,mH0).*CrystalBall2(optimal_x_HHmH1,mH1);
f2_4b = fitmH04b4boost(mH0).*fitmH14b4boost(mH1);

alphaAll = -0.1:0.0001:0.1; %Range of fraction Signal to BG to calculate for
LoglikeliHood = zeros(1,length(alphaAll));
for i=1:length(alphaAll)
    alpha = alphaAll(1,i);
    logProb = weights.*log(alpha.*f1_HH + (1-alpha).*f2_4b);
    LoglikeliHood(1,i) = -sum(logProb);
end
figure()
plot(alphaAll,LoglikeliHood)
title('LogL v. Alpha Boosted')
xlabel('alpha')
ylabel('Log(Likelihood)')

[minLogL,indLogL] = min(LoglikeliHood);
AlphaMin = alphaAll(indLogL);

minLogL1sigma = minLogL+0.5;
for i=2:length(LoglikeliHood)
    if LoglikeliHood(i-1) <= minLogL1sigma && LoglikeliHood(i) >= minLogL1sigma, fprintf('1 Sigma Lower bound is %.4f\n',alphaAll(i)); end
    if LoglikeliHood(i-1) >= minLogL1sigma && LoglikeliHood(i) <= minLogL1sigma, fprintf('1 Sigma Upper bound is %.4f\n',alphaAll(i)); end
end