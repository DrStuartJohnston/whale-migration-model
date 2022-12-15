%% Code for Johnston and Painter to save relevant data at specific time points.
%% Called by Homing_Script_WorldMap.m

xPosition(tSaveCount,iRepeat) = mean(position(:,1));                        % Mean position (x) of the population.
yPosition(tSaveCount,iRepeat) = mean(position(:,2));                        % Mean position (y) of the population.

distanceToGoal(tSaveCount,iRepeat) = sum(sqrt((position(:,1)-goalLocationX).^2 + ...
    (position(:,2)-goalLocationY).^2)/nIndividualsStart);                   % Distance to the target location

clusterMeasure(tSaveCount,iRepeat) = sum(pairDistanceVec)/max((nIndividuals*(nIndividuals-1)),1);  % Measure of population clustering.                                                                      
nNeighbours = zeros(nIndividuals,1);                                        % Current number of observed neighbours.
diffDirection = zeros(nIndividuals,1);                                      % Difference in direction between heading and target.

% Loop over individuals
for i = 1:nIndividuals
    % Calculate the detectable neighbours
    neighbours = find((signalDB(pairDistances(i,setdiff(1:nIndividuals,i))) - ...
            shippingNoise(closestShipLonIndex(i),closestShipLatIndex(i),shippingSnapshot)) > -minimumNoiseOverlap & ...
            signalDB(pairDistances(i,setdiff(1:nIndividuals,i))) > minimumHearing);
    nNeighbours(i) = numel(neighbours);                                     % Number of observed neighbours.
    diffDirection(i) = abs(heading(i)+pi - mod(navigationField(lonPosition(i),latPosition(i))+pi,2*pi));
    % Direction between heading and target.
    if diffDirection(i) > pi
        diffDirection(i) = pi - mod(diffDirection(i),pi);
    end
end

meanNeighbours(tSaveCount,iRepeat) = mean(nNeighbours);                                                 % Average number of observed neighbours.
meanDifferenceDirection(tSaveCount,iRepeat) = mean(diffDirection);                                      % Average difference between heading and target.
nIndividualsRemaining(tSaveCount,iRepeat) = nIndividuals;                                               % Number of individuals yet to arrive at target.
savePositionLat(tSaveCount,individualList,iRepeat) = latPosition';                                      % Store latitude.         
savePositionLon(tSaveCount,individualList,iRepeat) = lonPosition';                                      % Store longtitude.
savePositionLat(tSaveCount,removalStore,iRepeat) = savePositionLat(max(1,tSaveCount-1),removalStore);   % Remove latitude of arrived individuals.
savePositionLon(tSaveCount,removalStore,iRepeat) = savePositionLon(max(1,tSaveCount-1),removalStore);   % Remove longtitude of arrived individuals.     
tSaveCount = tSaveCount + 1;                                                                            % Increase counter of number of saved points.
directionHist = directionHist + histcounts(mod(heading,2*pi),linspace(0,2*pi,nHistDirection))';         % Generate histogram of headings.

% Check if 90% of individals have arrived at target and if so, store time.
if nIndividuals/nIndividualsStart <= 0.1 && majorityCheck == 0              
    majorityGone(iRepeat) = t;
    majorityCheck = 1;
end