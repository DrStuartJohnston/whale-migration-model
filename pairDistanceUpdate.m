% Update distances between pairs of agents in the simulation.

pairDistances = zeros(nIndividuals);
pairDistanceVec = pdist(position);                                  % Calculate distances between all pairs of individuals.
pairDistances(triu(ones(nIndividuals)==1,1)) = pairDistanceVec;     % Set pair distances for i =/= j.
pairDistances(tril(ones(nIndividuals)==1,-1)) = pairDistanceVec;    % Set pair distances for i =/= j.     
