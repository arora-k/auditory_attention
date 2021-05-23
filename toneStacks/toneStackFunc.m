function s = toneStackFunc(cen_freq)
%%%create tone stacks based on ERB spaced axis
%%%adapted from rippleFunc.m by Michelle Moerel originally for creating ripples
% Inputs:
% |----- cen_freq: frequency of central tone in the tone stack
% Outputs:
% |----- s: 1x16000 tone stack stimulus
close all;

%%fixed variables
fs = 16000;             %sampling freq; needs to be 16000
DM = 0.6;               %modulation depth
CF_old = 440*2 .^((-31:97)/24); %create frequency axis; this is the one used for the fMRI experiment sounds
CF = MakeErbCFs(CF_old(1),CF_old(end),10000);  %ERB frequency axis; stick to original axis as close as possible
%load('F.mat'); CF = F(2:end-1); %model channel freqs

%%variables to be set
rippleParam = struct;
rippleParam.freq = cen_freq;  %center frequencies of ripples (in Hz)
%rippleParam.rate = 3;      %temporal modulation frequency (in Hz)
%rippleParam.width = 1;          %width of ripple frequency band (in oct)
rippleParam.width = 10;          %width of ripple frequency band (in number of model channels)

%%initialize output
nrsounds = numel(rippleParam.width)*numel(rippleParam.freq)*numel(rippleParam.rate);
sounds = zeros(nrsounds,fs);
LookupVect = zeros(nrsounds,5);

%create dynamic-rippled-noise
isnd = 1;
for iwidth = 1:numel(rippleParam.width)
    width = rippleParam.width(iwidth);
    
    for ifreq = 1:numel(rippleParam.freq)
            
            %%create stack
            [ripple,widthOct] = DoStack(rippleParam.freq(ifreq),width,CF);
            
            %%normalize power
            en_s = sum(ripple.^2);
            ripple = ripple*sqrt(100/en_s);
            ripple = ripple- mean(ripple);
            sounds(isnd,:) = ripple';
    end
end

s = sounds;

function y=HzToErbRate(x)
% Convert Hz to ERB rate
%
%   y=HzToErbRate(x)
%
%   y = HzToErbRate(x) converts the frequency X (in Hz) to
%   the eqivalent ERB number.
%
%   See also ERBRATETOHZ, MAKEERBCFS.

% !---
% ==========================================================
% Last changed:     $Date: 2012-10-28 13:02:39 +0000 (Sun, 28 Oct 2012) $
% Last committed:   $Revision: 210 $
% Last changed by:  $Author: ch0022 $
% ==========================================================
% !---

y=(21.4*log10(4.37e-3*x+1));

% [EOF]
end
function cfs = MakeErbCFs(mincf,maxcf,numchans)
% Make a series of center frequencies equally spaced in ERB-rate.
%
%   cfs = MakeErbCFs(mincf,maxcf,numchans)
%
%   This function makes a vector of center frequenies
%   equally spaced on the ERB-rate scale.
%
%   cfs = MakeErbCFs(mincf,maxcf,numchans) creates numchans
%   centre frequencies between mincf and maxcf.
%
%   Adapted from code written by: Guy Brown, University of
%   Sheffield and Martin Cooke.
%
%   See also ERBRATETOHZ, HZTOERBRATE.

% !---
% ==========================================================
% Last changed:     $Date: 2012-10-28 13:02:39 +0000 (Sun, 28 Oct 2012) $
% Last committed:   $Revision: 210 $
% Last changed by:  $Author: ch0022 $
% ==========================================================
% !---

cfs = ErbRateToHz(linspace(HzToErbRate(mincf),HzToErbRate(maxcf),numchans));

% [EOF]
end
function y=ErbRateToHz(x)
% Convert ERB rate to Hz.
%
%   y = ErbRateToHz(x)
%
%   y = ErbRateToHz(x) converts the ERB number x to the
%   eqivalent frequency y (in Hz).
%
% See also HZTOERBRATE.

% !---
% ==========================================================
% Last changed:     $Date: 2012-10-28 13:02:39 +0000 (Sun, 28 Oct 2012) $
% Last committed:   $Revision: 210 $
% Last changed by:  $Author: ch0022 $
% ==========================================================
% !---

y=(10.^(x/21.4)-1)/4.37e-3;

% [EOF]
end
function [x] = RampSound(x,fs,ramptime)
%
%ramp sound onset and offset according to ramptime (in ms)
%

nrsamples = numel(x);
rampsamples = round(ramptime/1000*fs);
x_ramp = x;

rampON = linspace(0,1,rampsamples);
rampOFF = linspace(1,0,rampsamples);

x_ramp(1:rampsamples) = x(1:rampsamples).*rampON';
x_ramp(nrsamples-(rampsamples-1):nrsamples) = x(nrsamples-(rampsamples-1):nrsamples).*rampOFF';
x = x_ramp;
end

end
