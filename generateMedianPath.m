%% Code to generate a median trajectory from the location observations

% Note that large numbers of trajectories can be infeasible to analyse.

coarseness = 1;                                                 % Coarseness parameter (i.e. how many skipped observations)
pathsUsed = 10;                                                 % Total number of trajectories used
pathIndices = randsample(1:size(savePositionLat,2),pathsUsed);  % Randomly select trajectories

% Calculate and sort the directions of the initial trajectories
startGradients = (savePositionLat(1,pathIndices)-startLat)./(savePositionLon(1,pathIndices)-startLon);
[sortGradients,order] = sort(startGradients);
chosenStart = sortGradients(ceil(pathsUsed/2));                 % Select a trajectory in the middle of the spread
chosenIndex = order(ceil(pathsUsed/2));                         % Currently used trajectory

nPoints = (nSavePoints-1)/coarseness;                           % Number of points on the trajectory

medianPath = [savePositionLat(1,chosenIndex),savePositionLon(1,chosenIndex)];       % Initialise points on the median trajectory

% Ensure all lat/lon points are non-zero (i.e. if simulation finished
% early apply final location to all save points).
for i = 1:pathsUsed
    nOb = nnz(savePositionLat(:,pathIndices(i)));
    savePositionLat(nOb+1:end,i) = savePositionLat(nOb,i);
    savePositionLon(nOb+1:end,i) = savePositionLon(nOb,i);
end

lineSegments_latStart = savePositionLat(1:coarseness:end-coarseness,pathIndices);   % Lat start point of each trajectory segment
lineSegments_latEnd = savePositionLat(coarseness+1:coarseness:end,pathIndices);     % Lat end point of each trajectory segment
lineSegments_lonStart = savePositionLon(1:coarseness:end-coarseness,pathIndices);   % Lon start point of each trajectory segment
lineSegments_lonEnd = savePositionLon(coarseness+1:coarseness:end,pathIndices);     % Lon end point of each trajectory segment

allSegments = [lineSegments_latStart(:), lineSegments_lonStart(:), ...
    lineSegments_latEnd(:), lineSegments_lonEnd(:)];                                % Matrix of all start and end points
totalIndex = pathsUsed*nPoints;                                                     % Total number of data pointss

% Code here develped by U. Murat Erdem. See SI document for the full reference.
allIntersects = lineSegmentIntersect(allSegments,allSegments);                      % Calculate all line segment intersections

%% Generate median trajectory
i = 1;
while i < nPoints
    storeChosenIndex = chosenIndex;                                                 % Current trajectory index
    ignoreIndex = (chosenIndex-1)*nPoints+1:chosenIndex*nPoints;                    % Ignore all points on that trajectory (no self intersection)
    diffIndex = setdiff(1:totalIndex,ignoreIndex);                                  % Ignore all corresponding indices
    tmp = lineSegmentIntersect([lineSegments_latStart(i,chosenIndex), ...           % Find intersections with non-ignored line segment components
                                lineSegments_lonStart(i,chosenIndex), ...
                                lineSegments_latEnd(i,chosenIndex), ...
                                lineSegments_lonEnd(i,chosenIndex)], ...
                                allSegments(diffIndex,:));
    intersectIndices = find(tmp.intAdjacencyMatrix==1);                             % Find all intersections with current segment
    possibleIntersectLat = tmp.intMatrixX(intersectIndices);                        % Lat locations of intersections
    possibleIntersectLon = tmp.intMatrixY(intersectIndices);                        % Lon locations of intersections
    [~,closestIntersect] = min((lineSegments_latStart(i,chosenIndex)-possibleIntersectLat).^2 + ...
      (lineSegments_lonStart(i,chosenIndex)-possibleIntersectLon).^2);              % Find closest intersection
    % Check that intersections exist
    if ~isempty(closestIntersect)
        newPathLatStart = possibleIntersectLat(closestIntersect);                   % Lat location of new trajectory segment
        allSegments(diffIndex(intersectIndices(closestIntersect)),1) = possibleIntersectLat(closestIntersect);
        newPathLonStart = possibleIntersectLon(closestIntersect);                   % Lon location of new trajectory segment
        allSegments(diffIndex(intersectIndices(closestIntersect)),2) = possibleIntersectLon(closestIntersect);
        
        % Perturb start points of trajectory to avoid self-intersection
        lineSegments_latStart(i,chosenIndex) = possibleIntersectLat(closestIntersect) + ...
            0.00001*(lineSegments_latEnd(i,chosenIndex)-possibleIntersectLat(closestIntersect));
        lineSegments_lonStart(i,chosenIndex) = possibleIntersectLon(closestIntersect) + ...
            0.00001*(lineSegments_lonEnd(i,chosenIndex)-possibleIntersectLon(closestIntersect));
        
        % Included perturbed points in the matrix of all segments
        allSegments((chosenIndex-1)*nPoints+i,1) = lineSegments_latStart(i,chosenIndex);
        allSegments((chosenIndex-1)*nPoints+i,2) = lineSegments_lonStart(i,chosenIndex);
        
        % Find index of the newly chosen trajectory
        chosenIndex = ceil(diffIndex(intersectIndices(closestIntersect))/nPoints);
        
        % Update line segment component to start at the intersection
        lineSegments_latStart(mod(diffIndex(intersectIndices(closestIntersect))-1,nPoints)+1,chosenIndex) = possibleIntersectLat(closestIntersect);
        lineSegments_lonStart(mod(diffIndex(intersectIndices(closestIntersect))-1,nPoints)+1,chosenIndex) = possibleIntersectLon(closestIntersect);
        
        % Update median trajectory
        medianPath = [medianPath;newPathLatStart,newPathLonStart];
        
        % Update number of points in the trajectory
        i = mod(diffIndex(intersectIndices(closestIntersect))-1,nPoints)+1;
    
    % If no intersection on this segment, update median trajectory and move
    % to the next segment
    elseif isempty(closestIntersect)
        chosenIndex = storeChosenIndex;
        medianPath = [medianPath; lineSegments_latEnd(i,chosenIndex), lineSegments_lonEnd(i,chosenIndex)];
        i = i+1;
    end
end

% Calculate the point on the median trajectory that is closest to the
% target location and consider as final point on the trajectory.
[~,closestToGoal] = min((medianPath(:,1)-latGoal).^2 + (medianPath(:,2)-lonGoal).^2);
medianPath = medianPath(1:closestToGoal,:);

% Clear unnecessary variables
clear tmp
clear allIntersects