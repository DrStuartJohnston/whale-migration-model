%% Calculate depth at location based on interpolation from grid points

diffFlowLatGrid = latFlowGrid(2)-latFlowGrid(1);                                                    % Space between grid points in latitude.
diffFlowLonGrid = lonFlowGrid(2)-lonFlowGrid(1);                                                    % Space between grid points in longitude.

[~, tmpLatIndex] = min(abs(latFlowGrid-latPosition(nextAgent)));                                    % Closest latitude grid point.
latDifference = latFlowGrid(tmpLatIndex)-latPosition(nextAgent);                                    % Check if greater or less than grid point.
latWeights = [1-abs(latDifference/diffFlowLatGrid),abs(latDifference/diffFlowLatGrid)];             % Calculate interpolation weights.
if latDifference > 0                                                                                % Select grid indices.
    latIndices = [tmpLatIndex, min(tmpLatIndex+1,numel(latFlowGrid))];
else
    latIndices = [tmpLatIndex, max(tmpLatIndex-1,1)];
end

[~, tmpLonIndex] = min(abs(lonFlowGrid-lonPosition(nextAgent)));                                    % Closest longitude grid point.
lonDifference = lonFlowGrid(tmpLonIndex)-lonPosition(nextAgent);                                    % Check if greater or less than grid point.
lonWeights = [1-abs(lonDifference/diffFlowLonGrid),abs(lonDifference/diffFlowLonGrid)];             % Calculate interpolation weights.
if lonDifference > 0                                                                                % Select grid indices.
    lonIndices = [tmpLonIndex, min(tmpLonIndex+1,numel(lonFlowGrid))];
else
    lonIndices = [tmpLonIndex, max(tmpLonIndex-1,1)];
end

xVelocityAtLocation = latWeights(1)*lonWeights(1)*lonFlowVelocity(lonIndices(1),latIndices(1),flowSnapShot) + ...  % Flow based on interpolation weights
    latWeights(2)*lonWeights(1)*lonFlowVelocity(lonIndices(2),latIndices(1),flowSnapShot) + ...                % and location.
    latWeights(1)*lonWeights(2)*lonFlowVelocity(lonIndices(1),latIndices(2),flowSnapShot) + ...
    latWeights(2)*lonWeights(2)*lonFlowVelocity(lonIndices(2),latIndices(2),flowSnapShot);

yVelocityAtLocation = latWeights(1)*lonWeights(1)*latFlowVelocity(lonIndices(1),latIndices(1),flowSnapShot) + ...  % Flow based on interpolation weights
    latWeights(2)*lonWeights(1)*latFlowVelocity(lonIndices(2),latIndices(1),flowSnapShot) + ...                % and location.
    latWeights(1)*lonWeights(2)*latFlowVelocity(lonIndices(1),latIndices(2),flowSnapShot) + ...
    latWeights(2)*lonWeights(2)*latFlowVelocity(lonIndices(2),latIndices(2),flowSnapShot);