function [] = histogramPlotsErr(binEventSumErrNorm,binValue1,binValue2,name,nameVar1,nameVar2)

% Plot histograms to same style each time (ERRORS)

figure() %2d tile plot, error of each bin
histogram2('XBinEdges',binValue1,'YBinEdges',binValue2,'BinCounts',binEventSumErrNorm,'DisplayStyle','tile')
colorbar
title(sprintf('%s fractional errors', name))
xlabel(sprintf('%s', nameVar1))
ylabel(sprintf('%s', nameVar2))
zlabel('Number of Events')
%saveas(gcf,fullfile([pwd '/figures'], sprintf('%s-fracErrors', name)),'epsc')