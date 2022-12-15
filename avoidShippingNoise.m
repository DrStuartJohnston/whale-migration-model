%% If it is too noisy, avoid loud noise by heading in the quietest direction

latSlopeAtLocation = latWeights(1)*lonWeights(1)*shippingSlopeLat(lonIndices(1),latIndices(1)) + ...   % Slope in latitude direction 
    latWeights(2)*lonWeights(1)*shippingSlopeLat(lonIndices(2),latIndices(1)) + ...                    % based on interpolation weights
    latWeights(1)*lonWeights(2)*shippingSlopeLat(lonIndices(1),latIndices(2)) + ...                    % and location.
    latWeights(2)*lonWeights(2)*shippingSlopeLat(lonIndices(2),latIndices(2));

lonSlopeAtLocation = latWeights(1)*lonWeights(1)*shippingSlopeLon(lonIndices(1),latIndices(1)) + ...   % Slope in longitude direction 
    latWeights(2)*lonWeights(1)*shippingSlopeLon(lonIndices(2),latIndices(1)) + ...                    % based on interpolation weights
    latWeights(1)*lonWeights(2)*shippingSlopeLon(lonIndices(1),latIndices(2)) + ...                    % and location.
    latWeights(2)*lonWeights(2)*shippingSlopeLon(lonIndices(2),latIndices(2));

directionAwayFromNoise = mod(atan2(-latSlopeAtLocation,-lonSlopeAtLocation),2*pi);                      % Quietest direction.