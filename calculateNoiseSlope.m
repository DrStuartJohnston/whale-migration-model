%% Calculate the slope in the noise in latitude and longitude basis vectors.

diffShippingLatGrid = shippingLat(2)-shippingLat(1);                         % Spacing of latitude grid for noise.
diffShippingLonGrid = shippingLon(2)-shippingLon(1);                         % Spacing of longitude grid for noise.

% Calculate x and y values of the grid for noise.
[shippingXGrid,shippingYGrid] = projfwd(projection,repmat(shippingLat,1,numel(shippingLon)), ...
    repmat(shippingLon',numel(shippingLat),1));

% Calculate noise slope in longitude basis vector.
shippingSlopeLon = zeros(size(backgroundNoise)); 
shippingSlopeLon(2:end-1,:,:) = (backgroundNoise(3:end,:,:)-backgroundNoise(1:end-2,:,:))/(2*diffShippingLatGrid);
shippingSlopeLon(1,:,:) = (backgroundNoise(2,:,:)-backgroundNoise(1,:,:))/diffShippingLatGrid;
shippingSlopeLon(end,:,:) = (backgroundNoise(end,:,:)-backgroundNoise(end-1,:,:))/diffShippingLatGrid;

% Calculate noise slope in latitude basis vector.
shippingSlopeLat = zeros(size(backgroundNoise));
shippingSlopeLat(:,2:end-1,:) = (backgroundNoise(:,3:end,:)-backgroundNoise(:,1:end-2,:))/(2*diffShippingLonGrid);
shippingSlopeLat(:,1,:) = (backgroundNoise(:,2,:)-backgroundNoise(:,1,:))/diffShippingLonGrid;
shippingSlopeLat(:,end,:) = (backgroundNoise(:,end,:)-backgroundNoise(:,end-1,:))/diffShippingLonGrid;