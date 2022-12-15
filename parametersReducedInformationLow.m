includeWind = 'yes';                                % Include Farcas wind data
includeFlowField = 'no';                            % Include data about ocean currents
reduceInformation = 'yes';                          % Reduce inherent information due to background noise.

shippingNoise(shippingNoise<70|isnan(shippingNoise)) = 70; % Eliminates artefacts in data;
nIndividualsStart = 100;                            % Number of individuals in the simulation at start.
velocity = 6000;                                    % Speed of individuals (m/h).
runTime = 1;                                        % Mean reorientation time (h).
tEnd = 744;                                         % End of simulation (h).
alpha = 10/20;                                      % Weighting of observations for heading calculation.
beta = 10/20;                                       % Weighting of observations for concentration calculation.
sensingRange = 1e10;                                % Perceptual range of individuals (m) [unused].
depthSensingRange = sensingRange;                   % Perceptual range of individuals to detect depth (m) [unused].
callDB = 178;                                       % Source level of signal (dB)
signalDB = @(distance) callDB - 17.8*log10(distance);  % Signal level as a function of distance.
backgroundNoise = shippingNoise;                    % Ambient noise (dB)
minimumHearing = 90;                                % Minimum received level of signal to be detected.
noiseAvoidanceLevel = 120;                          % Noise avoidance threshold parameter.
depthAvoidanceLevel = 30;                           % Depth of water (m) that is deemed to be avoided.
minimumNoiseOverlap = 5;                            % Ability of individuals to detect noise (i.e. can detect signal at X DB below background)
backgroundStrength = 1;                             % Background information level.
repulsionDistance = 0;                              % Repulsion mechanism (unused).
alignDistance = sensingRange;                       % Alignment distance (always = sensing range).
attractDistance = sensingRange;                     % Attraction mechanism (unused).
goalDistance = 50000;                               % Distance from goal (m) to be counted as "arrived".
minimumInherentInformation = 0;                     % Minimum level of inherent information.
informationMidpoint = 110;                          % Midpoint in tanh function for inherent information.
informationDecay = 10;                              % Scaling in tanh function for inherent information.