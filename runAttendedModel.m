%% runAttendedModel.m
% Loads ripples from RippleBin.mat, plots output on those ripples using AC_Model_attend.m

load('RippleBin.mat')

%load ripple sounds
freq = [300 4000];
rate = [3 10];
ripples = [];

for ifreq = 1:numel(freq)
    for irate = 1:numel(rate)
        [ripples(irate + (ifreq-1)*numel(freq),:),Fs] = audioread(['Ripple_',...
            num2str(freq(ifreq)),'Hz_',num2str(rate(irate)),'Hz_energy1_ramp10.wav']);
    end
end

ripple_resp = cell(4,3);

for i = 1:size(ripples,1)
    for param = 1:3
        EE1 = AC_Model_attend(ripples(i,:)',param-1);
        ripple_resp_init_4000{i,param} = mean(EE1'); % Change 300/4000 depending on attended frequency in AC_Model_attend
        disp(strcat('Done with i=',num2str(i),' and param=',num2str(param)))
    end
end

%% plot results original, unaveraged

close all
title_1 = 'Case:';
title_2 = ["Unattended", "Parameter Set 1", "Parameter Set 2"];
title_3 = [", Attend Low", ", Attend High"];

%load('ripple_resp_init_300.mat'); load('ripple_resp_init_4000.mat')
figure;

for param = 1:3
    subplot(2,3,param)
    for i = 1:4
        plot(ripple_resp_init_300{i,param},'Linewidth',2)
        hold on
    end
    xlabel('Channels 1 to 98'); ylabel('Mean firing rate over time');
    title(strcat(title_1,title_2(param),title_3(1)));
    legend('300Hz/3Hz','300Hz/10Hz','4000Hz/3Hz','4000Hz/10Hz','location','Northwest');
    ylim([0 8])
end

for param = 1:3
    subplot(2,3,param+3)
    for i = 1:4
        plot(ripple_resp_init_4000{i,param},'Linewidth',2)
        hold on
    end
    xlabel('Channels 1 to 98'); ylabel('Mean firing rate over time');
    title(strcat(title_1,title_2(param),title_3(2)));
    legend('300Hz/3Hz','300Hz/10Hz','4000Hz/3Hz','4000Hz/10Hz','location','Northwest')
    ylim([0 8])
end

%% plot results, averaged

rippleRespAvgd300 = cell(2,3); rippleRespAvgd4000 = rippleRespAvgd300;

for i = 1:3
    rippleRespTemp = [ripple_resp_init_300{1,i}; ripple_resp_init_300{2,i}; ripple_resp_init_300{3,i}; ripple_resp_init_300{4,i}];
    rippleRespAvgd300{1,i} = mean(rippleRespTemp(1:2,:));
    rippleRespAvgd300{2,i} = mean(rippleRespTemp(3:4,:));
    
    rippleRespTemp = [ripple_resp_init_4000{1,i}; ripple_resp_init_4000{2,i}; ripple_resp_init_4000{3,i}; ripple_resp_init_4000{4,i}];
    rippleRespAvgd4000{1,i} = mean(rippleRespTemp(1:2,:));
    rippleRespAvgd4000{2,i} = mean(rippleRespTemp(3:4,:));
end

for param = 2:3
    figure;
    plot(rippleRespAvgd300{1,param}-rippleRespAvgd4000{1,param},'Linewidth',2)
    hold on
    plot(rippleRespAvgd300{2,param}-rippleRespAvgd4000{2,param},'Linewidth',2)
    
    xlabel('Channels 1 to 98'); ylabel('Attend Low-Attend High Firing Rate');
    title(strcat('Paramater set:',num2str(param-1)));
    legend('300Hz','4000Hz','location','Southwest')
end



