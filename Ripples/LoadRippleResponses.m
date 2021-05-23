
clear all;
close all

dirdata = '..';


%%load ripple sounds
freq = [300 4000];
rate = [3 10];
ripples = [];

for ifreq = 1:numel(freq)
    for irate = 1:numel(rate)
        [ripples(irate + (ifreq-1)*numel(freq),:),Fs] = audioread([dirdata,'\Ripples\Ripple_',num2str(freq(ifreq)),'Hz_',num2str(rate(irate)),'Hz_energy1_ramp10.wav']);
        
        %%switch on to play ripple
%         sound(ripples(irate + (ifreq-1)*numel(freq),:),Fs);
%         pause;
    end
end


%%load and plot fMRI responses to ripple sounds
load(['Ripplebin.mat']);
[nsubj,nMask,nMap,nFreqBin] = size(rippleBin); %8 subjects, 5 masks (see nameMSK), 5 sound maps (i.e., response to a sound throughout the cortex; see nameMSK), 15 BF frequency bins (see freqnew)

%these data can be plotted in many ways; below is one example
imsk = 4; %nameMSK{imsk} will tell you which mask (i.e., area of the brain) you will plot

figure(1)
ste_temp = std(squeeze(rippleBin(:,imsk,2,:)-rippleBin(:,imsk,4,:)))./sqrt(numel(subjects)-1); %Response to 300 Hz in Attend low (i.e., attend to 300 Hz) - Attend High (i.e., attend to 4 kHz)
errorbar(freqnew, mean(squeeze(rippleBin(:,imsk,2,:)-rippleBin(:,imsk,4,:))),ste_temp,'Linewidth',2);
hold all
ste_temp = std(squeeze(rippleBin(:,imsk,3,:)-rippleBin(:,imsk,5,:)))./sqrt(numel(subjects)-1);%Response to 4 kHz in Attend low - Attend High 
errorbar(freqnew, mean(squeeze(rippleBin(:,imsk,3,:)-rippleBin(:,imsk,5,:))),ste_temp,'Linewidth',2);

ylim = get(gca,'YLim'); plot([300 300],ylim,'k'); plot([4000 4000],ylim,'k'); set(gca,'fontsize', 12); %add vertical lines at 300 Hz and 4 kHz
legend({'Ripple 300 Hz', 'Ripple 4 kHz'});
set(gca,'xscale','log');
set(gca,'XTick',(2.^(linspace(log2(min(freqnew)), log2(max(freqnew)), min(numel(freqnew),5)))));
set(gca,'XTickLabel',round(2.^(linspace(log2(min(freqnew)), log2(max(freqnew)), min(numel(freqnew),5)))/100)/10);
ylabel('PSC ALow - AHigh');
xlabel('BF (kHz)');
title('fMRI ripple responses') %title(nameMSK{imsk});

hold on
plot(xlim,[0 0],'Handlevisibility','off')


