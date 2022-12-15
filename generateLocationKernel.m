%% Code to generate the spread in the individual trajectories.

kernelGrid = zeros(size(shippingXGrid));                                        % Density of individual trajectory locations
[shipLatGrid,shipLonGrid] = projinv(projection,shippingXGrid,shippingYGrid);    % Superimpose onto shipping noise grid
sd = 0.01;                                                                      % Spread in the density kernel
kernelThreshold = 5e-2;                                                         % Density threshold to be included in the spread
nObservations = nnz(savePositionLat);                                           % Number of observations

% Loop over each non-boundary location
for i = 2:numel(shippingLat)-1
    midLatL = (shippingLat(i)+shippingLat(i-1))/2;                              % Calculate lat midpoint of next grid square
    midLatR = (shippingLat(i)+shippingLat(i+1))/2;                              % Calculate lat midpoint of next grid square   
    for j = 2:numel(shippingLon)-1
        midLonL = (shippingLon(j)+shippingLon(j-1))/2;                          % Calculate lon midpoint of next grid square
        midLonR = (shippingLon(j)+shippingLon(j+1))/2;                          % Calculate lon midpoint of next grid square
        for k = 1:nIndividualsStart*nRepeats
            for l = 1:nSavePoints
                % If an individual is located between the four midpoints,
                % apply a Gaussian kernel at that location and evaluate at
                % all grid points.
                if savePositionLat(l,k) > midLatL && savePositionLat(l,k) < midLatR && ...
                        savePositionLon(l,k) > midLonL && savePositionLon(l,k) < midLonR
                    kernelGrid = kernelGrid + exp(-((shipLatGrid-shippingLat(i)).^2+(shipLonGrid-shippingLon(j)).^2)/sd)/sd;
                end                    
            end
        end
    end
end


kernelGrid = kernelGrid/nObservations;                              % Normalise by the number of observations
frequentLat = double(shipLatGrid(kernelGrid>kernelThreshold));      % Calculate the frequently observed grid points
frequentLon = double(shipLonGrid(kernelGrid>kernelThreshold));      % Calculate the frequently observed grid points
frequentBoundary = boundary([frequentLat,frequentLon],0.6);         % Generate a boundary from the frequently observed grid points
