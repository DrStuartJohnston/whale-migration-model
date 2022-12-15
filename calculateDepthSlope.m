%% Calculate the slope in the bathymetry in latitude and longitude basis vectors.

diffDepthLatGrid = depthLatGrid(2)-depthLatGrid(1);                         % Spacing of latitude grid for depth.
diffDepthLonGrid = depthLonGrid(2)-depthLonGrid(1);                         % Spacing of longitude grid for depth.

% Calculate x and y values of the grid for depth.
[depthXGrid,depthYGrid] = projfwd(projection,repmat(depthLatGrid',1,1001),repmat(depthLonGrid,1001,1));

% Calculate depth slope in latitude basis vector.
depthSlopeLat = zeros(size(depthGrid)); 
depthSlopeLat(2:end-1,:) = (depthGrid(3:end,:)-depthGrid(1:end-2,:))/(2*diffDepthLatGrid);
depthSlopeLat(1,:) = (depthGrid(2,:)-depthGrid(1,:))/diffDepthLatGrid;
depthSlopeLat(end,:) = (depthGrid(end,:)-depthGrid(end-1,:))/diffDepthLatGrid;

% Calculate depth slope in longitude basis vector.
depthSlopeLon = zeros(size(depthGrid));
depthSlopeLon(:,2:end-1) = (depthGrid(:,3:end)-depthGrid(:,1:end-2))/(2*diffDepthLonGrid);
depthSlopeLon(:,1) = (depthGrid(:,2)-depthGrid(:,1))/diffDepthLonGrid;
depthSlopeLon(:,end) = (depthGrid(:,end)-depthGrid(:,end-1))/diffDepthLonGrid;