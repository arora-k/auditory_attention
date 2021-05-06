function [a1Output_E] = AC_Model_orig(carrier_freq, param, channel_span)
%%% Modified version of Auditory Cortex Model (Description below)
%%% Last Edited: 06.05.21 by Kabir Arora
% Changes from original AC_Model.m:
% |----- Only includes A1
% |----- Added section to modify post-peripheral WCCM input
% Inputs:
% |----- carrier_freq: frequency (Hz) of pure sine wave used as model input
% |----- param: array of parameters, in the format [EIgain, sigmaEI, EEgain, sigmaEE]
% |----- channel_span: In case of multi-channel input, array of channels (1-98) to receive input
% Outputs:
% |----- a1Output_E: 98 x 16000 array of excitatory A1 output


%% Spectro-temporal processing in a two-stream computational model of the Auditory Cortex

% This script embedded functions from :
%       Wilson Cowan cortical model 
%       Wilson HR. Computation by excitatory and inhibitory networks. In: Spikes, Decisions & Actions: Dynamical Foundations of Neuroscience. Oxford UK: Oxford University Press; 1999. p. 88–115.
%
%       4th order gammatone filterbank, Ma N, Green P, Barker J, Coy A. Exploiting correlogram structure for robust speech recognition with multiple speech sources. Speech Commun. 2007;
%       Available at: http://www.dcs.shef.ac.uk/~ning/resources/gammatone/

%       LIN implementation based on NSL toolbox
%       Powen Ru, Multiscale Multirate Spectro-Temporal Auditory Model\
%       Available at: https://isr.umd.edu/Labs/NSL/Software.htm

% The script generates an amplitude modulated tone signal, which is then processed by the 4 simulated areas:
%                         A1 and R (core)
%                         Slow and Fast (belt)
% It computes strongest osciallation frequency across output channel closest to carrier frequency for all 4 areas, the Vector Strength of these oscillations and the mean firing rates.

% Author: Isma Zulfiqar

close all; tic
%disp(['Modeling spectro-temporal processing in 2 core (A1,R) and 2 belt (Slow, Fast) areas of the AC:' newline])

%% Initializing Frequency axis

Size = 100;  % Spatial (Frequency) size of array - number of units
lim1_frqaxis = 50;
lim2_frqaxis = 8000;
F = MakeErbCFs(lim1_frqaxis,lim2_frqaxis,Size);  % ERB frequency axis

%disp(['Model spans ' num2str(lim1_frqaxis) ' Hz to ' num2str(lim2_frqaxis) ' Hz on tonotopic axis...'])

%% Generating input sound 
duration = 1; % second
Fs = 16000;
mod_rate = 6; % Hz
mod_depth = 1; % from 0 to 1
%carrier_freq = 1000; % Hz

%disp(['Input signal is ' num2str(carrier_freq) ' Hz tone, amplitude modulated at ' num2str(mod_rate) ' Hz...'])

t = [0:(duration*Fs)-1]'/Fs;
s = sin(2*pi*carrier_freq*t);
power = sum(s.^2);
%figure(1); plot(t,s); xlabel('Time'); ylabel('Amplitude'); title('Sound waveform');

% =============================================================================================
% Peripheral Processing Stage
% =============================================================================================
  
%% Applying gammatone filterbank
%disp('Passing the sound through gammatone filterbank ... ')

bm_upsampled = gammatoneFast(s,F,Fs); 
%color='k'; figureSpec(1:length(bm_upsampled),F,squeeze(bm_upsampled'),2,'Gammatone filter output',color);

%% Lateral Inhibitory Network
%disp('Running tonotopic derivative and half-wave rectification ...')

y2_h = bm_upsampled(:,Size);
y3_h = 0; y4 = []; y5 =[]; v5 = [];
for ch = Size-1:-1:1
    y2 = bm_upsampled(:,ch);
    
    % ===== Lateral inhibition (masked by higher
    % (frequency) spatial response)
    
    y3 = y2 - y2_h;
    y2_h = y2;
    
    % ==== half wave rectification --- > y4
    y4(:,ch) = max(y3,0);
    
    
    % temporal integration window ---> y5
    y5 = filter(1, [1 -0.9845],y4(:,ch));
    v5(:,ch) = decimate(y5,64);
    
    
end
F_new = F(2:Size-1);
bm = y4(:,2:Size-1);
%color='k'; figureSpec(1:64:length(bm_upsampled),F_new,squeeze(v5(:,2:end)'),3,'Temporally integrated Spectrogram (LIN output)',color);

% converting LIN output to high reslution model input
Last = floor(1000 * (length(s)/Fs));  %last time in computation in ms
DelT = 1000/Fs;  % dt - time resolution of the simulations in ms
time = 0:DelT:Last;
step = DelT/1000; % ms
timei = 0:step:Last/1000;  % time axis with DelT for simulations (in s.)
FS = 1/step;
bm = bm(ceil(linspace(1,length(s),length(timei))),:);
bm = 1E3*bm;

Stim = bm';

%[~,channel] = min(abs(F-carrier_freq(1))); % computing channel with CF closes to carrier freq
    
% =============================================================================================
% Cortical Processing Stage
% =============================================================================================
    
%% inizialize the network and define the inputs
%disp('Initializing network...')

EE1 = zeros(Size-2, length(time)); IN1 = EE1; EEresp1 = EE1; INresp1 = EE1;
EE2 = zeros(Size-2, length(time)); IN2 = EE2; EEresp2 = EE2; INresp2 = EE2;
EE3 = zeros(Size-2, length(time)); IN3 = EE3; EEresp3 = EE3; INresp3 = EE3;
EE4 = zeros(Size-2, length(time)); IN4 = EE4; EEresp4 = EE4; INresp4 = EE4;
input1 = zeros(1,length(time)); input2 = EE1; input3 = EE1;

%% Stimulus
%disp('Initial input is only for excitatory units, inhibitory input is set to zero...')

P = Stim(:,1:length(time));  % input to excitatory network
%color='k'; figureSpec(timei,F_new,Stim,4,'Model Input', color); 
Q=zeros(Size-2,length(time));  % input to the IN population is set to zero

%% Constant input over time, single/multi channel

% P_new = zeros(98,16001);
% 
% mean_input = mean(P(55, :),'all');
% 
% for chan = 1:98
%     if (ismember(chan, channel_span))
%         P_new(chan, :) = mean_input*ones(size(P_new(chan, :)));
%     end
% end
% 
% P = P_new;%./(sqrt(numel(channel_span)));

%color='k'; figureSpec(timei,F_new,P,4,'Model Input', color);
%% Model Parameters
%disp('Setting model parameters and checking for stability...')

Xsyn1 = 20*(-5:5); % filter width
Xsyn2 = 20*(-5:5);
Xsyn3 = 20*(-2:2);
Xsyn4 = 20*(-2:2);

% naka rushton constants
% AI
thExc1 = 60;
thIn1 = 80;
max1 = 100; %100

% R
thExc2 = 60;
thIn2 = 80;
max2 = 100;

% slow
thExc3 = 60;
thIn3 = 80;
max3 = 100;

% fast
thExc4 = 60;%20;
thIn4 = 80;%40;
max4 = 100;%12.5;

% maximum synaptic strength for active transient response
EEgain = param(3); %1.5; 
EIgain = param(1); %1.3;
IIgain = 1.5; %1.5;
IEgain = EIgain; %1.3

% tau - time constant of the network in ms; 
DT1 = 10; % A1
DT2 = 20; % R
DT3 = linspace(300,200,Size-2)'; % Slow, Time constants change across tonotopic axis
DT4 = linspace(3,1,Size-2)'; % Fast, Time constants change across tonotopic axis

% spatial spread, EE = excitatory to excitatory, EI = excitatory to
% inhibitory, IE = inhibitory to excitatory, II = inhibitory to inhibitory

% A1
sigmaEE1 = param(4); %40; 
sigmaEI1 = param(2); %160;
sigmaIE1 = sigmaEI1;
sigmaII1 = 10; %10;
%disp('For A1:'); checkSigmas(sigmaEE1, sigmaEI1, sigmaIE1, sigmaII1, thExc1, thIn1, EEgain, EIgain, IEgain, IIgain , max1);

% R
sigmaEE2 = 40;
sigmaEI2 = 160;
sigmaIE2 = sigmaEI2;
sigmaII2 = 10;
%disp('For R:'); checkSigmas(sigmaEE2, sigmaEI2, sigmaIE2, sigmaII2, thExc2, thIn2, EEgain, EIgain, IEgain, IIgain , max2);

% Slow
sigmaEE3 = 20;
sigmaEI3 = 80;
sigmaIE3 = sigmaEI3;
sigmaII3 = 10;
%disp('For Slow:'); checkSigmas(sigmaEE3, sigmaEI3, sigmaIE3, sigmaII3, thExc3, thIn3, EEgain, EIgain, IEgain, IIgain , max3);

% Fast
sigmaEE4 = 200;
sigmaEI4 = 300;
sigmaIE4 = sigmaEI4;
sigmaII4 = 30; 
%disp('For Fast:'); checkSigmas(sigmaEE4, sigmaEI4, sigmaIE4, sigmaII4, thExc4, thIn4, EEgain, EIgain, IEgain, IIgain , max4);

% synaptic weights
synEE1 = EEgain*exp(-abs(Xsyn1)./sigmaEE1);
synEI1 = EIgain*exp(-abs(Xsyn1)./sigmaEI1);
synII1 = IIgain*exp(-abs(Xsyn1)./sigmaII1); 
%synIE1 = IEgain*exp(-abs(Xsyn1)./sigmaIE1);

synEE2 = EEgain*exp(-abs(Xsyn2)./sigmaEE2); 
synEI2 = EIgain*exp(-abs(Xsyn2)./sigmaEI2);
synII2 = IIgain*exp(-abs(Xsyn2)./sigmaII2); 

synEE3 = EEgain*exp(-abs(Xsyn3)./sigmaEE3); 
synEI3 = EIgain*exp(-abs(Xsyn3)./sigmaEI3);
synII3 = IIgain*exp(-abs(Xsyn3)./sigmaII3); 

synEE4 = EEgain*exp(-abs(Xsyn4)./sigmaEE4); 
synEI4 = EIgain*exp(-abs(Xsyn4)./sigmaEI4);
synII4 = IIgain*exp(-abs(Xsyn4)./sigmaII4); 

%% Simulating response

%disp('Simulating responses... ');
%% AI 


smoothing_kernel = [0.5 1 0.5]; % weight of connections from input to single unit
for T = 2:length(time)  %Loop in ms, Euler solution method
    
    input1 = conv(P(:,T),smoothing_kernel,'same')./sum(smoothing_kernel);
    
    EEresp1(:,T) = NeuralConv(synEE1, EE1(:,T-1)) - NeuralConv(synEI1, IN1(:,T-1)) + input1;
    EEresp1(:,T) = (EEresp1(:,T).*(EEresp1(:,T) > 0)).^2;
    INresp1(:,T) = NeuralConv(synEI1, EE1(:,T-1)) - NeuralConv(synII1, IN1(:,T-1)) + Q(:,T);
    INresp1(:,T) = (INresp1(:,T).*(INresp1(:,T) > 0)).^2;

    EE1(:,T) = EE1(:,T-1) + (DelT/DT1)*(-EE1(:,T-1) + (max1)*EEresp1(:,T)./(thExc1^2 + EEresp1(:,T)));
    IN1(:,T) = IN1(:,T-1) + (DelT/DT1)*(-IN1(:,T-1) + (max1)*INresp1(:,T)./(thIn1^2 + INresp1(:,T)));
end

a1Output_E = squeeze(EE1);

%color='b'; figureSpec(timei,F_new,squeeze(EE1),5,'A1 - Excitatory Response',color);
% color='r'; figureSpec(timei,F_new,squeeze(IN1),6,'A1 - Inhibitory Response',color);


%% Stability conditions for the model
function checkSigmas(sigmaEE, sigmaEI, sigmaIE, sigmaII, thExc, thIn, EEgain, EIgain, IEgain, IIgain, n)

for i=1:length(sigmaEE)
% condition 1: inhibition spread is larger than excitation
if ~((sigmaEI(i) == sigmaIE(i)) && (sigmaIE(i) > sigmaEE(i)))
    error('The sigmaIE/sigmaEI should be greater than sigmaEE.')
end

% condition 2: remove unstable states

if ~((n *( (2 * EEgain * sigmaEE(i)) - (2 * IEgain * sigmaIE(i)))) < thExc && ...
        (n *( (2 * EIgain * sigmaEI(i)) - (2 * IIgain * sigmaII))) > thIn)
    error(['Try changing values of parameter ' num2str(i)]);
end

end
    disp('Condition for activity spread satisfied by sigma values.')
    disp('Condition for stability of response satisfied by set parameters.')
end

%% Neural Convolution
function Result = NeuralConv(fltr, inputs)
% Convolves fltr with inputs and removes extraneous values
%fltr is assumed to have an odd number of elements and be centered
%Replicates inputs for periodic boundary conditions

% spread functions
Sz = length(inputs);
Xx = conv(fltr, inputs);
extra = fix(length(fltr)/2);
Xx = Xx(1 + extra: length(Xx) - extra);
Result = Xx;
end

%% Draw figures
function figureSpec(t,f,cdata1, fig, text, color)
%CREATEFIGURE(CDATA1)
%  CDATA1:  image cdata

%  Auto-generated by MATLAB on 22-Apr-2015 23:18:29
% Create figure

figure1 = figure(fig);
clf;

% Create axes
axes1 = axes('Parent',figure1,...
    'Position',[0.13 0.287878787878788 0.693788546255507 0.687121212121215],...
    'Layer','top');
%% Uncomment the following line to preserve the X-limits of the axes
xlim(axes1,[t(1) t(end)]);
%% Uncomment the following line to preserve the Y-limits of the axes
ylim(axes1,[f(1) f(end)]);

box(axes1,'on');
hold(axes1,'all');
set(axes1,'YTick',floor(linspace(f(1),f(end),5)));
ticks = floor(linspace(1,length(f),5));
set(axes1,'YTickLabel',{num2str(floor(f(ticks(1)))) num2str(floor(f(ticks(2)))) num2str(floor(f(ticks(3)))) num2str(floor(f(ticks(4)))) num2str(floor(f(ticks(5))))});
set(axes1,'XTickLabel',{'' '' '' '' '' ''});
set(axes1, 'FontSize', 12);
% Create image
imagesc(t,f,cdata1,'Parent',axes1,'CDataMapping','scaled');
ylabel('Frequency (Hz)');
title(text);
% colormap jet;

% Create axes
axes2 = axes('Parent',figure1,...
    'Position',[0.127063142437592 0.0867003367003374 0.695256975036711 0.138047138047137]);

box(axes2,'on');
p=plot(t,mean(cdata1,1),color); set(p, 'LineWidth', 2);
xlim(axes2,[t(1) t(end)]);
xlabel('Time (sec)'); set(axes2, 'FontSize', 12);
ylabel('Mean firing rate');

% Create axes
axes3 = axes('Parent',figure1,...
    'Position',[0.86568281938326 0.29040404040404 0.0814537444933923 0.684343434343435]);
box(axes3,'on');

p=plot(mean(cdata1,2),1:length(f),color); set(p, 'LineWidth', 2);
ylim(axes3,[1 length(f)]);
set(axes3,'YTick',length(f)/4:length(f)/4:length(f));
ticks = floor(length(f)/4:length(f)/4:length(f));
set(axes3,'YTickLabel',{'' '' '' ''});
set(axes3, 'FontSize', 12);
xlabel('Mean firing rate');


end

%% Gammatone filter functions
function [bm,env,delay] = gammatoneFast(x,cfs,fs,align)
% Produce an array of responses from a Gammatone filter via FFT
%   
%   bm = gammatoneFast(x,cfs,fs)
%   bm = gammatoneFast(...,align)
%   [bm,env] = gammatoneFast(...)
%   [bm,env,delay] = gammatoneFast(...)
%
%   This function takes an input vector and passes it
%   through a bank of fourth-order gammatone filters, with
%   centre frequencies specified by cfs. The function
%   returns a matrix, with each row/column corresponding to
%   a filter output with a centre frequency determined by
%   the corresponding element in cfs. The orientation of the
%   output is determined by the orientation of the input: if
%   x is a row vector then the output will contain one row
%   for each filter output, and vice versa.
% 
%   Centre frequencies may be any value below the Nyquist
%   rate (determined by the sampling frequency fs).
%   Typically centre frequencies are equally spaced on the
%   ERB-rate scale and may be calculated thus:
%
%   cfs = MakeErbCFs(low_cf,high_cf,numchans)
%   
%   where low_cf is the lowest frequency in the bank,
%   high_cf is the highest, and numchans is the numbers of
%   filters in the bank.
% 
%   bm = gammatoneFast(...,align) allows phase alignment to
%   be applied. With align=false, no alignment is applied
%   (default). With align=true, fine structure and envelope
%   alignment is applied so that the impulse response peaks
%   occurs at t=0.
% 
%   [bm,env] = gammatoneFast(...) returns the instantaneous
%   envelopes env for each filter.
% 
%   [bm,env,delay] = gammatoneFast(...) returns the delay
%   (in samples) removed by the phase alignment of each
%   gammatone filter, i.e. delays(n)=0 if align=false. delay
%   is a vector the same size as cfs.
%
%   Based on code written by ZZ Jin, adapted by DLW in
%   Jan'07 and JF Woodruff in Nov'08
% 
%   See also MakeErbCFs.

% !---
% ==========================================================
% Last changed:     $Date: 2012-10-28 13:02:39 +0000 (Sun, 28 Oct 2012) $
% Last committed:   $Revision: 210 $
% Last changed by:  $Author: ch0022 $
% ==========================================================
% !---

if nargin < 3
    fs = 16000; % default sampling frequency
end
if nargin < 4
    align = false; % default phase alignment
end

% check inputs
assert(isvector(x) & isnumeric(x),'x must be a vector')
assert(isvector(cfs) & isnumeric(cfs),'cfs must be a vector')
assert(isscalar(fs),'fs must be a scalar')
assert(islogical(align) & numel(align)==1,'align must be logical')

% number of frequency channels
numchans = length(cfs);

filterOrder = 4; % filter order
gL = 2^nextpow2(0.128*fs); % gammatone filter length at least 128 ms
b = 1.019.*24.7.*(4.37.*cfs./1000+1); % rate of decay or bandwidth

gt = zeros(gL,numchans);  % Initialise IR
tc = zeros(size(cfs));  % Initialise time lead
phase = 0;

tpt=(2*pi)/fs;
gain=((1.019.*b.*tpt).^filterOrder)./6; % based on integral of impulse

tmp_t = (0:gL-1)/fs;

% calculate impulse response
for i = 1:numchans
    if align
        tc(i) = (filterOrder-1)./(2*pi*b(i));
        phase = -2*pi*cfs(i)*tc(i);
    end
    gt(:,i) = gain(i)*fs^3*tmp_t.^(filterOrder-1).*exp(-2*pi*b(i)*tmp_t).*cos(2*pi*cfs(i)*tmp_t+phase);
end

% if input is row vector, transpose to column vector
rot = false;
if size(x,1)==1
    x = x';
    rot = true;
end

% gammatone filtering using FFTFILT
bm = fftfilt(gt,repmat(x,1,numchans));

% Hilbert envelope
env = abs(hilbert(bm));

% delay due to time lead
delay = round(tc.*fs);

% remove time lead
for i = 1:numchans
    bm(:,i) = [bm(delay(i)+1:end,i); zeros(delay(i),1)];
    env(:,i) = [env(delay(i)+1:end,i); zeros(delay(i),1)];
end

% transpose output if necessary
if rot
    bm = bm';
    env = env';
end

% [EOF]
end
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
end
