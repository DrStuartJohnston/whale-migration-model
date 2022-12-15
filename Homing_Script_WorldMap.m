%% Code to simulate the whale migration model described in "Avoidance, confusion
%% or isolation? Modelling the impact of noise pollution on whale migration" by
%% Johnston and Painter. This is the highest level script for the mode. Note this file
%% requires the look-up table for the concentration parameter, which can found
%% at https://melbourne.figshare.com/articles/dataset/kappaCDFLookupTable_mat/14551614.

clear all
close all

addpath('./ParameterSets/');
addpath('./SyntheticShippingRoutes/');
addpath('./DataToLoad/');

% Choose appropriate terrain map.
terrainMap = 'NorthSea';                                        % Choose terrain map, here always 'NorthSea'.
parameterSet = 'Pristine';                                      % Choose a parameter set

% Choices are 'Pristine' 'LargeNoise' 'NoNoiseAvoidance' 'HighSensitivityAvoidance'
% 'InterSensitivityAvoidance' 'LowSensitivityAvoidance' 'CurrentNoise' 'CurrentShipping'
% 'IncreasedShipping' 'IncreasedShippingSlowdown' 'ReducedInformationHigh'
% 'ReducedInformationLow' 'LandAvoidance'

load('HighDefNorthSea_RelevantOnly.mat');                       % Load land boundaries.
latLimit = [48 62];                                             % North sea map latitude limits.
lonLimit = [-16 9];                                             % North sea longitude limits.

nRepeats = 1;                                                   % Number of realisations of the model.
nSavePoints = 501;                                              % Number of time points to save model output.

load('kappaCDFLookupTable.mat');                                % Load the lookup table for estimating the vM concentration parameter.
load('shippingNoise.mat');                                      % Load noise map.

finalTime = zeros(nRepeats,1);                                  % Time for all individuals to arrive at the goal.
xPosition = zeros(nSavePoints,nRepeats);                        % Mean position (x) of the population.
yPosition = zeros(nSavePoints,nRepeats);                        % Mean position (y) of the population.
clusterMeasure = zeros(nSavePoints,nRepeats);                   % Measure of clustering of the population.
meanNeighbours = zeros(nSavePoints,nRepeats);                   % Mean number of neighbours within perceptual range.
distanceToGoal = zeros(nSavePoints,nRepeats);                   % Mean distance to goal of the population.
meanDifferenceDirection = zeros(nSavePoints,nRepeats);          % Mean error in heading relative to target.
nIndividualsRemaining = zeros(nSavePoints,nRepeats);            % Number of individuals remaining in the simulation (i.e. yet to arrive at goal).
majorityGone = zeros(nRepeats,1);                               % Time for 90% of the individuals to arrive at the goal.

backgroundFieldType = 'Fixed';                                  % Type of background field (always 'Fixed');
noiseInfluence = 'Information';                                 % Choose type of noise influence either 'Information' or 'Range' (always 'Information')
signalType = 'Intended';                                        % 'Intended signal direction of target, 'Actual' observe directions of other individuals (always 'Intended')
signalIfArrived = 'yes';                                        % Individuals indicate direction of target once they have arrived (always 'Yes')

% Load appropriate parameter set 
if strcmpi(parameterSet,'Pristine')                             % Results in Figure 1/2.
    parametersPristine;
elseif strcmpi(parameterSet,'LargeNoise')                       % Results in SI Figure 2.
    parametersLargeNoise;
elseif strcmpi(parameterSet,'NoNoiseAvoidance')                 % Results in SI Figure 1.
    parametersNoNoiseAvoidance;
elseif strcmpi(parameterSet,'HighSensitivityAvoidance')         % Results in Figure 3.
    parametersHighSensitivityAvoidance;
elseif strcmpi(parameterSet,'InterSensitivityAvoidance')        % Results in Figure 3.
    parametersInterSensitivityAvoidance;
elseif strcmpi(parameterSet,'LowSensitivityAvoidance')          % Results in Figure 3.
    parametersLowSensitivityAvoidance;
elseif strcmpi(parameterSet,'CurrentNoise')                     % Results in Figure 2.
    parametersCurrentNoise;
elseif strcmpi(parameterSet,'CurrentShipping')                  % Results in Figure 5.
    parametersCurrentShipping;
elseif strcmpi(parameterSet,'IncreasedShipping')                % Results in Figure 5.
    parametersIncreasedShipping;
elseif strcmpi(parameterSet,'IncreasedShippingSlowdown')        % Results in Figure 5.
    parametersIncreasedShippingSlowdown;
elseif strcmpi(parameterSet,'ReducedInformationHigh')           % Results in Figure 4.
    parametersReducedInformation;
elseif strcmpi(parameterSet,'ReducedInformationLow')            % Results in Figure 4.
    parametersReducedInformation;
elseif strcmpi(parameterSet,'LandAvoidance')                    % Results in SI Figure 3.
    parametersLandAvoidance;
else
    error('Choose an appropriate parameter set')
end

if strcmpi(includeFlowField,'no')                               % Set flow field multipler to zero if flow not included.
    flowField = 0;                                              % Flow field (unused).
    flowDirection = 0;                                          % Flow direction (unused).
    flowVelocity = 0;                                           % Flow velocity (unused).
else
    flowField = 1;
end

if strcmpi(includeWind,'yes')                                   % Combine wind noise with shipping noise
    load('WindNoise.mat');                                      % Load wind noise.
    meanWindNoise(meanWindNoise<70|isnan(meanWindNoise)) = 70;  % Remove data artefacts.
    windEnergy = 10.^(meanWindNoise/10)*1e-12;                  % Convert to energy.
    shippingEnergy = 10.^(shippingNoise/10)*1e-12;              % Convert to energy.
    totalEnergy = windEnergy+shippingEnergy;                    
    shippingNoise = 10*log10(totalEnergy/1e-12);                % Convert back to dB.
    clear windNoise                                             % Clear unused variables
    clear windEnergy
    clear totalEnergy
    clear shippingEnergy
end

depthAvoidanceWeight = @(x) 0.5-0.5*tanh(0.5*(x-depthAvoidanceLevel));      % Weight movement away from shallow water.
noiseAvoidanceWeight = @(x) 0.5+0.5*tanh(0.5*(x-noiseAvoidanceLevel));      % Weight movement away from noise signal.

totalStepCountLoop = 0;                                                     % Number of reorientation events.
nHistDirection = 60;                                                        % Number of points in histograms.
directionHist = zeros(nHistDirection-1,1);                                  % Predefine direction histogram.
cbar = [linspace(40,115,20)',linspace(36,213,20)',linspace(108,236,20)']/255; %Define colormap.
projection = projcrs(54030,'Authority','ESRI');                             % Appropriate lat-lon x-y projection.

% Calculate the relevant latitude and longitude ranges based on terrain
% map.
relevantLand = (coarseLat>(min(latLimit)-5)&coarseLat<(max(latLimit)+5)&coarseLon>(min(lonLimit)-5)&coarseLon<(max(lonLimit)+5))|isnan(coarseLat);

% Calculate x,y co-ordinates based on projection of latitude and
% longitude.
[coarseLandX,coarseLandY] = projfwd(projection,coarseLat(relevantLand),coarseLon(relevantLand));

domainWidth = lonLimit(2)-lonLimit(1);              % Width of the domain.
domainHeight = latLimit(2)-latLimit(1);             % Height of the domain.                       

plotEvery = 2000;                                   % Plot results every X hours.
nHours = size(shippingNoise,3);                     % Number of hours in the data

%% Load relevant information for the appropriate scenario.
if strcmpi(terrainMap,'NorthSea')
    latGoal = 61;                                   % Lat location of target.
    lonGoal = -5;                                   % Lon location of target.
    [goalLocationX,goalLocationY] = projfwd(projection,latGoal,lonGoal);  % Calculate x,y position of target location.
    startLat = 53;                                  % Starting latitude for north sea map.
    startLon = 4;                                   % Starting longitude for north sea map.
    load('NorthSeaFlowJuly1July10.mat');            % Load current data from HYCOM.
    NorthSeaFlow = hawaii_soest_7e38_7a7b_afxhffc2; 
    lonFlowGrid = NorthSeaFlow.longitude;           % Grid points for latitude for currents.
    latFlowGrid = NorthSeaFlow.latitude;            % Grid points for longitude for currents.
    lonFlowVelocity = permute(NorthSeaFlow.water_u,[4,3,1,2])*3600; % Convert longitude velocity to lat/lon/time.
    latFlowVelocity = permute(NorthSeaFlow.water_v,[4,3,1,2])*3600; % Convert latitude velocity to lat/lon/time.
    lonFlowVelocity(isnan(lonFlowVelocity)) = 0;    % Remove NaNs.
    latFlowVelocity(isnan(latFlowVelocity)) = 0;    % Remove NaNs.
    nDays = size(latFlowGrid,3);                    % Number of days in the current data.
    
    % Bathymetry functionality
    load('BathymetryGrid.mat');                     % Load bathymetry data.
    depthGrid = abs(depthGrid);                     % Convert to depth below water.
    [depthXGrid,depthYGrid] = projfwd(projection,repmat(depthLatGrid',1,1001),repmat(depthLonGrid,1001,1)); % Calculate x and y values of the grid for depth.
    calculateDepthSlope;                            % Function to calculate slope in bathymetry
    
else
    error('Choose an appropriate terrain map')
end

navigationField = @(x,y) atan2(goalLocationY-y,goalLocationX-x) ;         % Direction of target.

calculateNoiseSlope;                                                      % Function to calculate slope in noise field.

savePositionLon = zeros(nSavePoints,nIndividualsStart,nRepeats);          % Store longitude position.
savePositionLat = zeros(nSavePoints,nIndividualsStart,nRepeats);          % Store latitude position.
flowTimeStep = 24;                                                        % Number of hours per timestep in the current data.

%% Perform nRepeats realisations of the simulation
for iRepeat = 1:nRepeats
    checkPlot = 0;                          % Counting variable to determine when to plot results.
    majorityCheck = 0;                      % Check if 90% of population has arrived at the target.
    iRepeat                                 % Print the realisation.
    
    nIndividuals = nIndividualsStart;       % Number of indivuduals in the simulation.
    
    t = 0;                                  % Initialise time counter.
    
    defineBackgroundFields;                 % Define noise and background fields.
   
    removalStore = [];                                          % Keep track of individuals to remove.
    
    initialPosition = zeros(nIndividuals,2);
    initialPosition(:,2) = startLat+rand(nIndividuals,1)-0.5;   % Initial latitude of individuals.
    initialPosition(:,1) = startLon+rand(nIndividuals,1)-0.5;   % Initial longitude of individuals.
    individualList = 1:nIndividuals;                            % List of individuals.
    lonPosition = initialPosition(:,1);                         % Longitude position of individuals.
    latPosition = initialPosition(:,2);                         % Latitude position of individuals.
    pairDistances = zeros(nIndividuals);
    [xPositionStart,yPositionStart] = projfwd(projection,initialPosition(:,2),initialPosition(:,1)); %x,y locations of individuals.
    position = [xPositionStart,yPositionStart];                         % Position of individuals.
    pairDistanceVec = pdist(position);                                  % Calculate distances between all pairs of individuals.
    pairDistances(triu(ones(nIndividuals)==1,1)) = pairDistanceVec;     % Set pair distances for i =/= j.
    pairDistances(tril(ones(nIndividuals)==1,-1)) = pairDistanceVec;    % Set pair distances for i =/= j.
    
    turningTime = exprnd(runTime,nIndividuals,1);                       % Calculate durations of run events.
    timeToUpdate = turningTime;                                         % Calculate time until reorientation events.
    
    heading = zeros(nIndividuals,1);                                    % Headings of individuals.
    
    % Sample individual headings based on inherent information.
    for i = 1:nIndividuals
        heading(i) = circ_vmrnd(navigationField(position(i,1),position(i,2)), ...
            navigationStrengthField(lonPosition(i),latPosition(i)),1);
    end
    
    communicatedDirection = heading;                                    % Direction that individuals believe the target location is.
    meanPosition = zeros(tEnd,2);                                       % Calculate mean position of the population.
    tSave = linspace(0,tEnd,nSavePoints);                               % Time points where the data will be saved.
    tSaveCount = 1;                                                     % Count of time points saved.
    totalStepCount = 0;                                                 % Number of steps taken.
    
    % Main loop of the individual simulations, run until end of simulation
    % or all individuals have arrived at the target.
    while t < tEnd && nIndividuals > 0
        shippingSnapshot = mod(ceil(t/timeStep)-1,nHours)+1;
        flowSnapShot = mod(ceil(t/flowTimeStep)-1,nDays)+1;
        % Uncomment for visualisation of simulation while running (can slow
        % down simulation).
        if t > checkPlot
            %plotWhileRunning;                                           
        end
        
        totalStepCount = totalStepCount + 1;                            % Keep track of total steps taken.
        
        [nextUpdate,nextAgent] = min(timeToUpdate);                     % Calculate next reorientation event time and agent.
        timeToUpdate = timeToUpdate - nextUpdate;                       % Update time to update for all individuals.
        timeElapsed = nextUpdate;                                       % Calculate time step length.
        t = t+nextUpdate;                                               % Update time.

        oldPosition = position;                                         % Retain previous positions.
        
        calculateFlowAtLocation;                                        % Calculate flow field at current location.
        
        % Update the position of all individuals.
        position = position + velocity*timeElapsed*[cos(heading),sin(heading)] + flowField*timeElapsed*[xVelocityAtLocation,yVelocityAtLocation];
        
        % Check that no individuals crossed onto land.
        checkBoundary;
        [latPosition,lonPosition] = projinv(projection,position(:,1),position(:,2)); %Calculate latitude and longitude from x,y positions
        pairDistanceUpdate;                                             % Update pair distances for all pairs of individuals.
        pairDistances(1:nIndividuals+1:end) = 1e10;                     % Avoid influence of pairs of identical individuals.
        closestShipLatIndex = zeros(nIndividuals,1);                    % Calculate closest index on shipping grid
        closestShipLonIndex = zeros(nIndividuals,1);
        % Calculate the closest index on the shipping noise map to the
        % location of the agent.
        for iIndividuals = 1:nIndividuals
            [~,closestShipLatIndex(iIndividuals)] = min(abs(shippingLat-latPosition(iIndividuals)));
            [~,closestShipLonIndex(iIndividuals)] = min(abs(shippingLon-lonPosition(iIndividuals)));
        end
        % Find individuals within the perceptual range of the individual
        % undergoing reorientation.
        
        calculateNoiseAtLocation;                                       % Calculate noise from shipping traffic
        avoidShippingNoise;                                             % Calculate direction away from shipping traffic    
        calculateDepthAtLocation;                                       % Calculate depth at current location
        avoidShallowWater;                                              % Calculate direction away from shallow water
        noiseAtNextAgent = backgroundNoiseAtLocation;                   % Noise at location
        depthAtNextAgent = depthAtLocation;                             % Depth at location
        neighbours = find((signalDB(pairDistances(nextAgent,:)) - noiseAtNextAgent) ...
              > -minimumNoiseOverlap & pairDistances(nextAgent,:) < 1e9 & ...
              signalDB(pairDistances(nextAgent,:)) > minimumHearing); % Check which signals can be detected over background noise
        nNeighbours = numel(neighbours);                                % Number of individuals within perceptual range.
        [minDistance,closestAgent] = min(pairDistances(nextAgent,:));   % Find closest agent.
        
        % Calculate level of inherent information available at the location
        % in the absence of noise.
        currentInherentInformation = navigationStrengthField(lonPosition(nextAgent),latPosition(nextAgent));
        
        % Reduce information due to confusion if appropriate.
        if strcmpi(reduceInformation,'yes')
           currentInherentInformation = minimumInherentInformation + (backgroundStrength-minimumInherentInformation) * ...
               (0.5 - 0.5*tanh((noiseAtNextAgent-informationMidpoint)/informationDecay));
        end
        
        % Calculate sample heading based on inherent information only.
        potentialHeading = circ_vmrnd(navigationField(position(nextAgent,1),position(nextAgent,2)),...
             currentInherentInformation,1);

            % Update heading based on other observed individuals if number of
            % neighbours exceeds zero.
            if nNeighbours > 0
                % Repulsion mechanism unused in this work.
                if minDistance < repulsionDistance
                    heading(nextAgent) = atan2(-position(closestAgent,2)+position(nextAgent,2), ...
                        -position(closestAgent,1)+position(nextAgent,1));
                % Alignment mechanism.
                elseif minDistance < alignDistance
                    successSignal = [];
                    
                    % Increase neighbours if cooperative navigation present (i.e.
                    % if individuals signal their success).
                    if strcmpi(signalIfArrived,'yes')                       % Use information if individuals signal their success.
                        distanceFromGoal = sqrt((position(nextAgent,1)-goalLocationX)^2+(position(nextAgent,2)-goalLocationY)^2);   % Distance to target.
                        if (signalDB(distanceFromGoal) - noiseAtNextAgent > -minimumNoiseOverlap) && ...
                                (signalDB(distanceFromGoal) > minimumHearing)                                                       % Check if success signal can be detected.
                            successSignal = ones(nIndividualsStart-nIndividuals,1);
                            successDirection = navigationField(position(nextAgent,1),position(nextAgent,2));
                            successSignal = successDirection*successSignal;
                            nNeighbours = nNeighbours + numel(successSignal);                                                       % Update neighbours to include successfully arrived.
                        end
                    end
                    
                    % Synthesise information from others based on whether
                    % signal is intent or actual behaviour.
                    if strcmpi(signalType,'Intended')
                       	bestGuessHeading = circ_mean([circ_mean([communicatedDirection(neighbours);successSignal]);potentialHeading],[1-alpha;alpha]); % MLE of heading.
                        alphaLookup = [communicatedDirection(neighbours);successSignal;potentialHeading];                                 % Set of observed headings.
                    elseif strcmpi(signalType,'Actual')
                        bestGuessHeading = circ_mean([circ_mean([heading(neighbours);successSignal]);potentialHeading],[1-alpha;alpha]);    % MLE of heading.
                        alphaLookup = [heading(neighbours);successSignal;potentialHeading];                                               % Set of observed headings.
                    end
                    
                    w = [(1-beta)*ones(nNeighbours,1);beta*nNeighbours];                                % Weighting of observed headings.
                    circ_kappa_script;                                                                  % Calculate estimate of concentration parameter.
                    bestGuessStrength = kappa;                                                          % Estimate of concentration parameter.
                    avoidNoiseWeight = noiseAvoidanceWeight(noiseAtNextAgent);                          % Calculate weighting of noise avoidance response.
                    avoidShallowWaterWeight = depthAvoidanceWeight(depthAtNextAgent);                   % Calculate weighting of land avoidance response.
                    
                    % Renormalise weights if necessary - if both noise and
                    % land avoidance response required.
                    if avoidNoiseWeight+avoidShallowWaterWeight > 1
                        totalWeight = avoidNoiseWeight+avoidShallowWaterWeight;
                        avoidNoiseWeight = avoidNoiseWeight/totalWeight;
                        avoidShallowWaterWeight = avoidShallowWaterWeight/totalWeight;
                    end
                    
                    % Combine avoidance and navigation responses to select
                    % new heading.
                    tmpHeading = circ_vmrnd(bestGuessHeading,bestGuessStrength,1);
                    communicatedDirection(nextAgent) = tmpHeading;
                    heading(nextAgent) = circ_mean([tmpHeading;directionAwayFromNoise;directionOfDeepestWater], ...
                        [1-avoidNoiseWeight-avoidShallowWaterWeight;...
                        avoidNoiseWeight;avoidShallowWaterWeight]);           % Set new heading.
                    
                % Attraction mechanism unused.
                elseif minDistance < attractDistance
                    heading(nextAgent) = atan2(mean(position(neighbours,2))-position(nextAgent,2),...
                        mean(position(neighbours,1))-position(nextAgent,1));
                end
            else % Otherwise individual navigation behaviour
                
                avoidNoiseWeight = noiseAvoidanceWeight(noiseAtNextAgent);              % Calculate weighting of noise avoidance response.
                avoidShallowWaterWeight = depthAvoidanceWeight(depthAtNextAgent);       % Calculate weighting of land avoidance response.
                
                % Renormalise weights if necessary
                if avoidNoiseWeight+avoidShallowWaterWeight > 1
                    totalWeight = avoidNoiseWeight+avoidShallowWaterWeight;
                    avoidNoiseWeight = avoidNoiseWeight/totalWeight;
                    avoidShallowWaterWeight = avoidShallowWaterWeight/totalWeight;
                end
                
                % Combine avoidance and navigation behaviour.
                communicatedDirection(nextAgent) = potentialHeading;
                heading(nextAgent) = circ_mean([potentialHeading;directionAwayFromNoise;directionOfDeepestWater], ...
                        [1-avoidNoiseWeight-avoidShallowWaterWeight;...
                        avoidNoiseWeight;avoidShallowWaterWeight]); 
            end

        timeToUpdate(nextAgent) = exprnd(runTime,1);                        % New duration of run.
        pairDistances(1:nIndividuals+1:end) = 0;                            % Set pair distances to zeros for identical individuals.
        
        % Storage of data at specific time points
        if t > tSave(tSaveCount)
            saveDataWorldMap;
        end
        
        % Determine which individuals have arrived at the target and remove
        % from simulation.
        removal = [];
        for i = 1:nIndividuals
            if sqrt((position(i,1)-goalLocationX)^2+(position(i,2)-goalLocationY)^2) < goalDistance
                removal = [removal;i];
            end
        end
        removalStore = [removalStore;individualList(removal)'];             % Store list of removed individuals.
        individualList(removal) = [];                                       % Remove individuals from list.
        position(removal,:) = [];                                           % Remove individuals from position.
        heading(removal) = [];                                              % Remove individuals from heading.
        communicatedDirection(removal) = [];                                % Remove individuals from communicatedDirection.
        timeToUpdate(removal) = [];                                         % Remove individuals from reorientation.
        nIndividuals = nIndividuals - numel(removal);                       % Number of individuals remaining.
    end
    
    finalTime(iRepeat) = t;                                                 % Final time in the simulation.
    
end

xPositionMean = mean(xPosition,2);                                          % Mean of average x position across realisation loop.
clusterMeasure = mean(clusterMeasure,2);                                    % Mean of clustering across realisation loop.
stdDistanceToGoal = std(distanceToGoal,0,2);                                % Standard deviation of distance to goal across realisation group.
distanceToGoal = mean(distanceToGoal,2);                                    % Mean of average distance to goal across realisation loop.
coarseNeighbours = zeros(floor(tEnd/24),nRepeats);

% Calculate daily average of detected neighbours
for i = 1:floor(tEnd/24)
   coarseNeighbours(i,:) = mean(meanNeighbours((i-1)*floor(nSavePoints*24/tEnd)+1:i*floor(nSavePoints*24/tEnd),:),1);
end

stdNeighbours = std(meanNeighbours,0,2);                                    % Standard deviation of detected neighbours across realisation loop.
stdCoarseNeighbours = std(coarseNeighbours,0,2);                            % Standard deviation of daily average of detected neighbours across realisation loop.
meanNeighbours = mean(meanNeighbours,2);                                    % Mean of average number of neighbours across realisation loop.
meanCoarseNeighbours = mean(coarseNeighbours,2);                            % Mean of daily average of detected neighbouts across realisation loop.
meanDifferenceDirection = mean(meanDifferenceDirection,2);                  % Mean of difference between heading and target across realisation loop.
stdNIndividualsRemaining = std(nIndividualsRemaining,0,2);                  % Standard deviation of number individuals remaining across realisation loop.
nIndividualsRemaining = mean(nIndividualsRemaining,2);                      % Mean of number individuals remaining across realisation loop.

% Reshape lat/lon saved locations into a 2D matrix
savePositionLon = reshape(savePositionLon,nSavePoints,nIndividualsStart*nRepeats);
savePositionLat = reshape(savePositionLat,nSavePoints,nIndividualsStart*nRepeats);

clear kappaCDF                                                              % Clear CDF to avoid saving over and over.

% Plot relevant results
plotResults;