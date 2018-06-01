function [] = histogramPlots(binEventSum,binValue1,binValue2,name,nameVar1,nameVar2)

% Plot histograms to same style each time (EVENTS)

figure() %2d tile plot, weight value shown as a colour
histogram2('XBinEdges',binValue1,'YBinEdges',binValue2,'BinCounts',binEventSum,'FaceColor','flat')%,'DisplayStyle','tile')
colorbar
title(sprintf('%s', name))
xlabel(sprintf('%s', nameVar1))
ylabel(sprintf('%s', nameVar2))
zlabel('Number of Events')
%saveas(gcf,fullfile([pwd '/figures'], sprintf('%s', name)),'epsc')
