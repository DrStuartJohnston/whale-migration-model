%% Calculate depth at location based on interpolation from grid points

[~, tmpLatIndex] = min(abs(depthLatGrid-latPosition(nextAgent)));                                   % Closest latitude grid point.
latDifference = depthLatGrid(tmpLatIndex)-latPosition(nextAgent);                                   % Check if greater or less than grid point.
latWeights = [1-abs(latDifference/diffDepthLatGrid),abs(latDifference/diffDepthLatGrid)];           % Calculate interpolation weights.
if latDifference > 0                                                                                % Select grid indices.
    latIndices = [tmpLatIndex, min(tmpLatIndex+1,numel(depthLatGrid))];
else
    latIndices = [tmpLatIndex, max(tmpLatIndex-1,1)];
end

[~, tmpLonIndex] = min(abs(depthLonGrid-lonPosition(nextAgent)));                                   % Closest longitude grid point.
lonDifference = depthLonGrid(tmpLonIndex)-lonPosition(nextAgent);                                   % Check if greater or less than grid point.
lonWeights = [1-abs(lonDifference/diffDepthLonGrid),abs(lonDifference/diffDepthLonGrid)];           % Calculate interpolation weights.
if lonDifference > 0                                                                                % Select grid indices.
    lonIndices = [tmpLonIndex, min(tmpLonIndex+1,numel(depthLonGrid))];
else
    lonIndices = [tmpLonIndex, max(tmpLonIndex-1,1)];
end

depthAtLocation = latWeights(1)*lonWeights(1)*depthGrid(latIndices(1),lonIndices(1)) + ...  % Depth based on interpolation weights
    latWeights(2)*lonWeights(1)*depthGrid(latIndices(2),lonIndices(1)) + ...                % and location.
    latWeights(1)*lonWeights(2)*depthGrid(latIndices(1),lonIndices(2)) + ...
    latWeights(2)*lonWeights(2)*depthGrid(latIndices(2),lonIndices(2));
