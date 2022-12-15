%% Code for Johnston and Painter to plot relevant information. Called by Homing_Script.m.

% Plot average x position
figure(3); hold on; plot(tSave/24,xPosition); title('Average X Position');
xlabel('Time'); ylabel('Average X Position'); box on;

% Plot measure of population clustering
figure(4); hold on; plot(tSave/24,clusterMeasure); title('Measure of Clustering');
xlabel('Time'); ylabel('Measure of Clustering'); box on;

% Plot average distance to target location
figure(5); hold on; plot(tSave/24,distanceToGoal); title('Average Distance to Goal');
xlabel('Time'); ylabel('Average Distance to Goal'); box on;
patch([tSave fliplr(tSave)]/24, [(distanceToGoal+stdDistanceToGoal)' ...
fliplr((distanceToGoal-stdDistanceToGoal)')], [0 0 1], ...
'LineStyle', 'none','FaceAlpha', 0.2); ylim([0 12e5]); xlim([0 31]);

% Plot the fraction of the population that are detectable neighbours
figure(6); hold on; plot(tSave/24,meanNeighbours./nIndividualsRemaining); title('Average Proportion of Neighbours'); xlabel('Time'); ylabel('Average Proportion of Neighbours'); box on;
tCoarse = 0.5:floor(tEnd/24)-0.5;

% Plot the daily average of the number of detectable neighbours
figure(7); hold on; plot(tCoarse,meanCoarseNeighbours); title('Average Neighbours'); xlabel('Time'); ylabel('Average Neighbours'); box on;
patch([tCoarse fliplr(tCoarse)], [(meanCoarseNeighbours+stdCoarseNeighbours)' ...
    fliplr((meanCoarseNeighbours-stdCoarseNeighbours)')], [0 0 1], ...
'LineStyle', 'none','FaceAlpha', 0.2); ylim([0 nIndividualsStart]); xlim([0 floor(tEnd/24)]);

% Plot the fraction of individuals that have yet to reach the target
figure(8); hold on; plot(tSave/24,nIndividualsRemaining/nIndividualsStart);
title('Proportion of Individuals Remaining'); xlabel('Time'); ylabel('Proportion of Individuals Remaining'); box on;

% Plot the number of individuals that have yet to reach the target
figure(9); hold on; plot(tSave/24,nIndividualsRemaining); title('Number of Individuals Remaining'); xlabel('Time'); ylabel('Number of Individuals Remaining'); box on;
patch([tSave fliplr(tSave)]/24, [(nIndividualsRemaining+stdNIndividualsRemaining)' ...
    fliplr((nIndividualsRemaining-stdNIndividualsRemaining)')], [0 0 1], ...
'LineStyle', 'none','FaceAlpha', 0.2); ylim([0 nIndividualsStart]); xlim([0 31]);

%% Setup world map
f10 = figure(10);
set(gcf,'position',[0,0,1000,1000])
worldmap(latLimit,lonLimit);

%% Noise background plot 
geoshow(repmat(shippingLat',numel(shippingLon),1),repmat(shippingLon,1,numel(shippingLat)), ...
    shippingNoise(:,:,shippingSnapshot),'DisplayType','surface','CData',shippingNoise(:,:,shippingSnapshot), ...
    'ZData',zeros(size(shippingNoise(:,:,shippingSnapshot))));
hold on;

%% Ocean current direction map
if ~exist('endLatFlow','var')
    endLatFlow = latFlowVelocity(:,:,flowSnapShot);
    endLonFlow = lonFlowVelocity(:,:,flowSnapShot);
end
s = 12; quiverm(repmat(latFlowGrid(s/2:s:end),1,numel(lonFlowGrid(s/2:s:end)))',repmat(lonFlowGrid(s/2:s:end)',numel(latFlowGrid(s/2:s:end)),1)', ...
    endLatFlow(s/2:s:end,(s/2:s:end)),endLonFlow(s/2:s:end,(s/2:s:end)),'r')


%% Median and spread of trajectories
if ~exist('frequentLat','var')
    generateLocationKernel;                         % Generate the spread in the trajectories
end

% Plot the region of the frequently observed locations
geoshow(frequentLat(frequentBoundary),frequentLon(frequentBoundary), ...
    'DisplayType','polygon','FaceColor','w','facealpha',.3,'EdgeColor','w')

if ~exist('medianPath','var')
    generateMedianPath;                             % Generate the median trajectory
end

% Plot the median trajectory
geoshow(medianPath(:,1),medianPath(:,2),'Color','w','linewidth',2)

%% Plot land map
colormap(cbar)                                                              
geoshow(landLat,landLon,'Color',[114 191 68]/255,'LineStyle','none','Marker','.','MarkerSize',5);   % Overlay land locations.
geoshow(newWorldLat,newWorldLon,'Color',[123 247 65]/255,'linewidth',2);                            % Overlay land locations.
geoshow(latGoal,lonGoal,'Color','r','LineStyle','none','Marker','.','MarkerSize',125)               % Draw goal location.
geoshow(latGoal,lonGoal,'Color','w','LineStyle','none','Marker','.','MarkerSize',100)               % Draw goal location.
geoshow(latGoal,lonGoal,'Color','r','LineStyle','none','Marker','.','MarkerSize',75)                % Draw goal location.
geoshow(latGoal,lonGoal,'Color','w','LineStyle','none','Marker','.','MarkerSize',50)                % Draw goal location.
geoshow(latGoal,lonGoal,'Color','r','LineStyle','none','Marker','.','MarkerSize',25)                % Draw goal location.

% Clear all variables except those needed to generate figures
clearvars -except latGoal lonGoal frequentLat frequentLon frequentBoundary cbar ...
    latFlowGrid lonFlowGrid endLatFlow endLonFlow tSave medianPath ...
    xPosition clusterMeasure distanceToGoal meanNeighbours nIndividualsRemaining ...
    nIndividualsStart latLimit lonLimit savePositionLat savePositionLon flowSnapShot ...
    stdNIndividualsRemaining stdNeighbours stdDistanceToGoal stdCoarseNeighbours ...
    meanCoarseNeighbours latStart lonStart nSavePoints tEnd
