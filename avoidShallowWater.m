%% If it is too shallow, avoid shallow water by heading in the direction of deepest water

latSlopeAtLocation = latWeights(1)*lonWeights(1)*depthSlopeLat(latIndices(1),lonIndices(1)) + ...   % Slope in latitude direction 
    latWeights(2)*lonWeights(1)*depthSlopeLat(latIndices(2),lonIndices(1)) + ...                    % based on interpolation weights
    latWeights(1)*lonWeights(2)*depthSlopeLat(latIndices(1),lonIndices(2)) + ...                    % and location.
    latWeights(2)*lonWeights(2)*depthSlopeLat(latIndices(2),lonIndices(2));

lonSlopeAtLocation = latWeights(1)*lonWeights(1)*depthSlopeLon(latIndices(1),lonIndices(1)) + ...   % Slope in longitude direction 
    latWeights(2)*lonWeights(1)*depthSlopeLon(latIndices(2),lonIndices(1)) + ...                    % based on interpolation weights
    latWeights(1)*lonWeights(2)*depthSlopeLon(latIndices(1),lonIndices(2)) + ...                    % and location.
    latWeights(2)*lonWeights(2)*depthSlopeLon(latIndices(2),lonIndices(2));

directionOfDeepestWater = mod(atan2(latSlopeAtLocation,lonSlopeAtLocation)+2*pi,2*pi);              % Direction of deepest water.