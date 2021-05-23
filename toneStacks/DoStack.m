function [ripple,widthOct] = DoStack(freq,nbins,CF)
%%%called by toneStackFunc.m
%%%adapted from DoRipple.m by Michelle Moerel originally for creating ripples

fs = 16000;
t = 0:1/fs:1;
NoTime = length(t)-1;  t=t(1:NoTime);

%set the starting freq, and number of freq bins (based on 'freq' and 'width')
[~, CFthis] = min(abs(freq-CF));
%CFthis = (CFthis-(ceil((nbins-1)/2)):1:CFthis+(ceil((nbins-1)/2)));
CFthis = 60*(-(ceil((nbins-1)/2)):1:(ceil((nbins-1)/2))) + CFthis*ones(1,nbins+1);
widthOct = log2(CF(CFthis(end))/CF(CFthis(1))); %width of ripple in oct
NoFreq = numel(CFthis);
%f0 = CF(CFthis(1)); %starting freq

ripple = zeros(1,length(t));
phase = pi/2;
for n = 1:NoFreq
    ripple = ripple + sin(2*pi*CF(CFthis(n))*t+phase);
end
