function [peak, tw] = runPureFTCs(param) % param = [EIgain, sigmaEI, EEgain, sigmaEE]
%%% Runs AC_Model to calculate and plot channel FTCs.
% Inputs:
% |----- param: array of parameters, in the format [EIgain, sigmaEI, EEgain, sigmaEE]
% Outputs:
% |----- peak: FTC peak (Maximum firing rate value)
% |----- tw: FTC tw in octaves

% Load center frequencies of units to use as starting point
load('../F.mat')

% Channels to inspect
channels = [55]; % 1:98 for full

% Data points on tuning curve
n_runs = 300; %300 recommended minimum

responses = zeros(98,n_runs);
responsesAxis = zeros(98,n_runs);
peaks = zeros(numel(channels),1); tws = peaks; 

% Obtain channel responses
for j = 1:numel(channels)

    car_freq = F(channels(j)+1); % Current center frequency (CF), Hz   

    % Prepare range of frequencies around CF to use as input
    left_end = (2*car_freq)/(sqrt(2)+1.5);
    right_end = (sqrt(2)+0.5)*left_end;
    freq_range = logspace(log10(left_end),log10(right_end),n_runs);
    responsesAxis(channels(j),:) = freq_range;

    % Record unit activity around CF
    for i = 1:numel(freq_range)
        [a1OutputEx_Temp, ~] = AC_Model_orig(freq_range(i)', param);
        %[a1OutputEx_Temp, a1OutputIn_Temp] = AC_Model_orig(0.05*sin(2*pi*freq_range(i)*t));
        close all
        activityDistributionE = mean(a1OutputEx_Temp');
        responses(channels(j),i) = activityDistributionE(1,channels(j));
    end
end

% Plot FTCs
figure(1);
for i = 1:numel(channels)
    plot(responsesAxis(channels(i),:),responses(channels(i),:))
    hold on
end

title('Frequency Tuning Curves'); xlabel('Frequency (Hz)'); ylabel('A1 Excitatory Activity'); 

% Calculate Tuning Widths
for i = 1:numel(channels)
   % Find half-peak index+value
   [max_v, max_i] = max(responses(channels(i),:));
   peaks(i) = max_v;
   half_peak = max_v/2;
   
   % Find indices of half-peak value
   left_set = responses(channels(i),1:max_i);
   right_set = responses(channels(i),max_i+1:end);
   [~, half_left_i] = min(abs(half_peak*ones(size(left_set))-left_set));
   [~, half_right_i] = min(abs(half_peak*ones(size(right_set))-right_set));
   
   % Convert to Octaves
   tws(i) = log2(responsesAxis(channels(i),half_right_i+numel(left_set))/responsesAxis(channels(i),half_left_i));
   disp(strcat('Done with number=',num2str(i)))  
end
end
