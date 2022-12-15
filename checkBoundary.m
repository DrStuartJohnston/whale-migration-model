% Check if individual is located with coarse polygon approximation of the
% coastline.
crossedBoundary = inpolygon(position(:,1),position(:,2),coarseLandX,coarseLandY);

% If so, abort the jump.
position(crossedBoundary,:) = oldPosition(crossedBoundary,:);