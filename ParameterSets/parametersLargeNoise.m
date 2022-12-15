shippingNoise = 70*ones(301,281,1440);              % Flat noise field
tmp = repmat(1:51,51,1);
shippingNoise(200:250,160:210,:) = repmat(300 - sqrt((26-tmp).^2+(26-tmp').^2),1,1,1440);
% Set region of extremely high noise

includeWind = 'no';                                 % Include Farcas wind data
includeFlowField = 'yes';                           % Include data about ocean currents
reduceInformation = 'no';                           % Reduce inherent information due to background noise.

nIndividualsStart = 10;                             % Number of individuals in the simulation at start.
velocity = 6000;                                    % Speed of individuals (m/h).
runTime = 1;                                        % Mean reorientation time (h).
tEnd = 744;                                         % End of simulation (h).
alpha = 10/20;                                      % Weighting of observations for heading calculation.
beta = 10/20;                                       % Weighting of observations for concentration calculation.
sensingRange = 1e10;                                % Perceptual range of individuals (m) [unused].
depthSensingRange = sensingRange;                   % Perceptual range of individuals to detect depth (m) [unused].
callDB = 178;                                       % Source level of signal (180)
signalDB = @(distance) callDB - 17.8*log10(distance);  % Signal level as a function of distance.
backgroundNoise = shippingNoise;                    % Ambient noise (DB)
minimumHearing = 88;                                % Minimum received level of signal to be detected.
noiseAvoidanceLevel = 120;                          % Noise avoidance threshold parameter.
depthAvoidanceLevel = 30;                           % Depth of water (m) that is deemed to be avoided.
minimumNoiseOverlap = 5;                            % Ability of individuals to detect noise (i.e. can detect signal at X DB below background)
backgroundStrength = 2;                             % Background information level.
repulsionDistance = 0;                              % Repulsion mechanism (unused).
alignDistance = sensingRange;                       % Alignment distance (always = sensing range).
attractDistance = sensingRange;                     % Attraction mechanism (unused).
goalDistance = 50000;                               % Distance from goal (m) to be counted as "arrived".