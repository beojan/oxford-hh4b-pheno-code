function [] = histogram1Dboth(binSum,binSumErr,binValue1,binValue2,name,nameVar1,nameVar2)

% Sum in each variable to create 2 seperate 1 dimensional plots
% Also correctly calculates the 1D errors

% histogram('BinEdges',edges,'BinCounts',counts)

% Errors:
binSumErr1 = sqrt(sum(binSumErr.^2,1));
binSumErr2 = sqrt(sum(binSumErr.^2,2));

binValue1mid = zeros(length(binSum)-1,1);
binValue2mid = zeros(length(binSum)-1,1);

% Errorbar Coordinates
for i=2:length(binValue1)
    binValue1mid(i-1) = 0.5*(binValue1(i-1) + binValue1(i));
    binValue2mid(i-1) = 0.5*(binValue2(i-1) + binValue2(i));
end

% var1:
figure()
hold on
histogram('BinEdges',binValue1,'BinCounts',sum(binSum,1),'facealpha',0.5);
errorbar(binValue1mid,sum(binSum,1),binSumErr1,'bo');
hold off
title(sprintf('%s - Distribution in %s', name, nameVar2))
xlabel(sprintf('%s', nameVar2))
ylabel('Number of Events')
%saveas(gcf,fullfile([pwd '/figures'], sprintf('%s-1D-distIn-%s', name, nameVar2)),'epsc')

% var2:
figure()
hold on
histogram('BinEdges',binValue2,'BinCounts',sum(binSum,2),'facealpha',0.5);
errorbar(binValue2mid,sum(binSum,2),binSumErr2,'bo');
hold off
title(sprintf('%s - Distribution in %s', name, nameVar1))
xlabel(sprintf('%s', nameVar1))
ylabel('Number of Events')
%saveas(gcf,fullfile([pwd '/figures'], sprintf('%s-1D-distIn-%s', name, nameVar1)),'epsc')

