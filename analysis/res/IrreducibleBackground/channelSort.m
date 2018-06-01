function [sampleRes,sampleInt,sampleBoost] = channelSort(sample)
%Sort into 3 tables -> res - int - boost

vectorRes = strcmp(string(sample{:,2}),'res'); %analysis_channel column
vectorInt = strcmp(string(sample{:,2}),'inter');
vectorBoost = strcmp(string(sample{:,2}),'boost');

sampleRes = sample;
sampleInt = sample;
sampleBoost = sample;

sampleRes(vectorRes(:)~=1,:) = [];
sampleInt(vectorInt(:)~=1,:) = [];
sampleBoost(vectorBoost(:)~=1,:) = [];

        
        
