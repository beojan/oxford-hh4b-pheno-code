% Likelihood test for intermediate channel

%From the training sample we get the distributions:
% HH: Multiplication of the CB functions: CrystalBall2(optimal_x_HHmH0,mH0)*CrystalBall2(optimal_x_HHmH1,mH1)
% 4b: Multiplcation of the normalised 1D fits:
% fitmH04b4int(mH0)*fitmH14b4int(mH1)

fitmH04b4intSUM = sum(sum(bin4b4intFitmH0)); %exp, only normalise 'a'
fitmH04b4int.a = fitmH04b4int.a/fitmH04b4intSUM;
fitmH14b4intSUM = sum(sum(bin4b4intFitmH1)); %polynomial, normalise 'a','b','c'
fitmH14b4int.a = fitmH14b4int.a/fitmH14b4intSUM;
fitmH14b4int.b = fitmH14b4int.b/fitmH14b4intSUM;
fitmH14b4int.c = fitmH14b4int.c/fitmH14b4intSUM;


weights = [int4b4tagSampleWeight;intHH4tagSampleWeight]; %testing weights
mH0 = [int4b4tagSamplemH0;intHH4tagSamplemH0]; %testing mH0
mH1 = [int4b4tagSamplemH1;intHH4tagSamplemH1]; %testing mH1

%To save running a loop for the sum, calculating fHH and f4b here:
f1_HH = CrystalBall2(optimal_x_HHmH0,mH0).*CrystalBall2(optimal_x_HHmH1,mH1);
f2_4b = fitmH04b4int(mH0).*fitmH14b4int(mH1);

alphaAll = -0.1:0.0001:0.1; %Range of fraction Signal to BG to calculate for
LoglikeliHood = zeros(1,length(alphaAll));
for i=1:length(alphaAll)
    alpha = alphaAll(1,i);
    logProb = weights.*log(alpha.*f1_HH + (1-alpha).*f2_4b);
    LoglikeliHood(1,i) = -sum(logProb);
end
figure()
plot(alphaAll,LoglikeliHood)
title('LogL v. Alpha Intermediate')
xlabel('alpha')
ylabel('Log(Likelihood)')

[minLogL,indLogL] = min(LoglikeliHood);
AlphaMin = alphaAll(indLogL);

minLogL1sigma = minLogL+0.5;
for i=2:length(LoglikeliHood)
    if LoglikeliHood(i-1) <= minLogL1sigma && LoglikeliHood(i) >= minLogL1sigma, fprintf('1 Sigma Lower bound is %.4f\n',alphaAll(i)); end
    if LoglikeliHood(i-1) >= minLogL1sigma && LoglikeliHood(i) <= minLogL1sigma, fprintf('1 Sigma Upper bound is %.4f\n',alphaAll(i)); end
end