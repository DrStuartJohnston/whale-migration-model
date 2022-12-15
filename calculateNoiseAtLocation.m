%% Calculate noise at location based on interpolation from grid points

[~, tmpLatIndex] = min(abs(shippingLat-latPosition(nextAgent)));                                    % Closest latitude grid point.
latDifference = shippingLat(tmpLatIndex)-latPosition(nextAgent);                                    % Check if greater or less than grid point.
latWeights = [1-abs(latDifference/diffShippingLatGrid),abs(latDifference/diffShippingLatGrid)];     % Calculate interpolation weights.
if latDifference > 0                                                                                % Select grid indices.
    latIndices = [tmpLatIndex, min(tmpLatIndex+1,numel(shippingLat))];
else
    latIndices = [tmpLatIndex, max(tmpLatIndex-1,1)];
end

[~, tmpLonIndex] = min(abs(shippingLon-lonPosition(nextAgent)));                                    % Closest longitude grid point.
lonDifference = shippingLon(tmpLonIndex)-lonPosition(nextAgent);                                    % Check if greater or less than grid point.
lonWeights = [1-abs(lonDifference/diffShippingLonGrid),abs(lonDifference/diffShippingLonGrid)];     % Calculate interpolation weights.
if lonDifference > 0                                                                                % Select grid indices.
    lonIndices = [tmpLonIndex, min(tmpLonIndex+1,numel(shippingLon))];
else
    lonIndices = [tmpLonIndex, max(tmpLonIndex-1,1)];
end

backgroundNoiseAtLocation = latWeights(1)*lonWeights(1)*backgroundNoise(lonIndices(1),latIndices(1),shippingSnapshot) + ...  % Noise based on interpolation weights
    latWeights(2)*lonWeights(1)*backgroundNoise(lonIndices(2),latIndices(1),shippingSnapshot) + ...                % and location.
    latWeights(1)*lonWeights(2)*backgroundNoise(lonIndices(1),latIndices(2),shippingSnapshot) + ...
    latWeights(2)*lonWeights(2)*backgroundNoise(lonIndices(2),latIndices(2),shippingSnapshot);