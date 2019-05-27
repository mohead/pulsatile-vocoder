function [stCI] = createCI()
%% function [stCI] = createCI()
% This function is a library of function to generate electrodogramms and auralization of these for cochlear implants.
% It should allow for auralization of arbitrary electrodogramms of cochlear implant systems.
%
% The auralization stage contains a function to simulate the spatial spread
% of excitation, followed by a Gammatone-analysis and synthesis filterbank,
% which simulates the frequency to place mismatch due to mismatched
% physiological best frequencies and the electrode-array place in the
% cochlea.
%
% Input parameters:
%       None, this function will return function handels, which have their
%       own input parameters-description
% Optional input parameters:
%	All optional input parameters can be evoked with 'key'-'value' pairs.
%   See function setParameters for the details
%
% Output parameters:
%	stCI = struct containing handles to the following functions:
%     stCI.Simulation = takes a mono acoustic signal and will generate a electrodogramm.
%     stCI.Auralisation = takes an electrodogramm as input and will
%     auralize it. You will have to know the bandwise rms and overall rms afterwards.
%     stCI.setParameter = Define Parameters for the CI Simulation.
%     stCI.testGammatonedelayline = function to test, wheter the delay of
%     the synthesis filterbank is ok or not (when the delay is too short,
%     interference might result, when adding the signals together. It is
%     optimal, when each first maxima of each channel occurs at the same
%     time).
%
% Author : Ben Williges <ben.williges@uni-oldenburg.de>
%%=======================================================================%%

% add necessary folders to MATLABs PATH
addpath([fileparts(which(mfilename)),filesep,'Gammatonefilterbank']);

stCI.Simulation = @CI_Simulation;
stCI.Auralisation = @CI_auralisation;
stCI.setParameter = @setParameter;
stCI.testGammatonedelayline = @test_gammatonefilterbankdelay;

end
%% Begin of main functions

function [electrodogramm, CISIM_parameter] = CI_Simulation(signal,fs, coding_strat, CISIM_parameter)
% function [electrodogramm, CISIM_parameter] = CI_Simulation(signal,fs, coding_strat, CISIM_parameter)
% This function will create an electrodogramm, given an monaural input signal sampled at fs. The coding strategy used
% is defined by the string coding_strat. The struct CISIM_parameter contains additional parameters, which are set by
% the inputParser-function setParameter().
%
% Input parameters:
%       signal = one channel vector containing the signal in samples
%       fs = sampling frequency of signal in Hz
%       coding_strat = coding strategy as a string. See subfunction coding_strategy() for a full list.
%       Currently only 'CI' is supported, which does a CIS processing
%       CISIM_parameter = struct containing additional parameters, like center frequencies of the filter bank, MCL, TCL, Volume, pps and so on. For the complete list and
%       for creating the parameter structure use the inputParser-function setParameter.
%
% Output parameters:
%       electrodogramm = electrodogramm of the signal in the dimension channels x time.
%       CISIM_parameter = Struct, which contains the parameters. It will also contains some additional parameters, like
%       the signals after the filterbank, in case you want to analyse them further.

% Preemphasis filter, see Phd-Thesis Laneau, 2006, When the deaf listen to music - Pitch perception in Cochlear implants, p. 51
w     = 2*1200/fs;
[b,a] = butter(1,w,'high');
signal     = filter(b,a,signal);
% Resampling to vocoder-frequency
if fs ~= CISIM_parameter.voc_sampling_frequency_hz;
    signal = resample(signal,CISIM_parameter.voc_sampling_frequency_hz,fs);
end

% Create analyse -Filterbank
analyzer_stim = Gfb_Analyzer_new(CISIM_parameter.voc_sampling_frequency_hz,...
    CISIM_parameter.gamma_order_stimulation, ...
    CISIM_parameter.center_frequencies_hz_stimulation,...
    CISIM_parameter.bandwidth_factor);

% First analysis-filterbank
%1st analyse:
[y1, analyzer_stim] = Gfb_Analyzer_process(analyzer_stim,signal');

%CI simulation:

% Calculate envelope after first filterbank
env_signal_new = abs(y1); %Envelope
weights = repmat(CISIM_parameter.weights,1,size(env_signal_new,2));
env_signal_new = sqrt(weights.*(real(y1).^2+imag(y1).^2));
fine_signal_new = real(y1); %Fine structure information

CISIM_parameter.env = env_signal_new;

% calculate rms per analysis-channel, which is needed later for auralization:
CISIM_parameter.rms_per_channel = rms2(y1,2);

% calculate pps (not all strategies might need that)
CISIM_parameter.block_delay = calculate_delay(CISIM_parameter.voc_sampling_frequency_hz,CISIM_parameter.pps);
%calculate pulse-length in samples:
CISIM_parameter.len_pulse = ceil((2*CISIM_parameter.pulselength+CISIM_parameter.ipg)*CISIM_parameter.voc_sampling_frequency_hz); %Round towards infinity to make sure, that exact timing is always possible

% lowpass filter envelope:
env_signal_new = lowpass_filter_env(env_signal_new, CISIM_parameter.center_frequencies_hz_stimulation, CISIM_parameter.voc_sampling_frequency_hz);
CISIM_parameter.env_lp = env_signal_new;

% sample envelope or fine structure depending on coding strategy
[electrodogramm, CISIM_parameter] = coding_strategy(env_signal_new, fine_signal_new, coding_strat,CISIM_parameter);

% Compression of pulses and mapping to Clinical Units
electrodogramm = process_compression_ci(electrodogramm,CISIM_parameter.B,CISIM_parameter.M,CISIM_parameter.alpha_c);
electrodogramm = converttoCU(electrodogramm, CISIM_parameter.TCL, CISIM_parameter.MCL, CISIM_parameter.Volume);

%debug: plot electrodogramm:
if CISIM_parameter.debug
    electrodogramm_plot = electrodogramm;
    % fine structure
    plotChannels(fine_signal_new, CISIM_parameter.voc_sampling_frequency_hz, 1,'Finestructure');
    xlabel('Time [s]');
    ylabel('Channel');
    ylim([0.8 12.5]);
    % original envelope signal
    plotChannels(env_signal_new, CISIM_parameter.voc_sampling_frequency_hz, 1, 'Envelope');
    xlabel('Time [s]');
    ylabel('Channel');
    % sampled data
    plotChannels(electrodogramm_plot, CISIM_parameter.voc_sampling_frequency_hz, 1/1200, strcat('Coding strategy: ', coding_strat));
    xlabel('Time [s]');
    ylabel('Elektrode');
    ylim([0.8 12.5]);
end

end
%%
function [outsig] = CI_auralisation(electrodogramm,fs,auralisation_parameter)
% function [outsig] = CI_auralisation(electrodogramm,fs,auralisation_parameter)
% This function will auralise a electrodogramm (channel x time) to a  monaural audio signal at sampling rate fs. The struct auralisation_parameter contains the parameters
% used for this function, like frequency to place mapping, spatial spread and smallband acoustic rms-values for each channel.
% Input parameters:
%       electrodogramm = electrodogramm of the signal in the dimension channels x time.
%       fs = sampling frequency of signal in Hz
%       CISIM_parameter = struct containing additional parameters, like frequency to place mapping, spatial spread and smallband acoustic rms-values for each channel.
%
% Output parameters:
%       outsig = Output signal, resampled to fs.

% convert CUs back to normal amplitudes (inverts loudness growth function/ CI compression):
electrodogramm = inverseCUConversion(electrodogramm,auralisation_parameter.TCL, ...
    auralisation_parameter.MCL, auralisation_parameter.Volume);
electrodogramm = inverseCICompression(electrodogramm,auralisation_parameter.B, ...
    auralisation_parameter.M, ...
    auralisation_parameter.alpha_c);
% Simulate spatial spread
% Interaction between electrodes
[electrodogramm,mspatial_weights] = CIinteract(electrodogramm, auralisation_parameter.max_elecs_affected, auralisation_parameter.distance_electrodes, auralisation_parameter.lambda);
% Create analyse -Filterbank
analyzer_aural = Gfb_Analyzer_new(auralisation_parameter.voc_sampling_frequency_hz,...
    auralisation_parameter.gamma_order_auralisation,...
    auralisation_parameter.center_frequencies_hz_auralisation,...
    auralisation_parameter.bandwidth_factor);
% Next Analyse-Filterbank (modulates puls-pattern to desired frequency region on basiliar
% membrane)
electrodogramm = Gfb_Analyzer_process(analyzer_aural, electrodogramm);
% Scale to correct channel-rms
[electrodogramm, auralisation_parameter.gain_per_channel] = setSignaltoRMS(electrodogramm,auralisation_parameter.rms_per_channel,'linear');
% Signal-Synthesis:
synthesizer_aural = Gfb_Synthesizer_new (analyzer_aural, 1/100); %1/100 = Delay
if auralisation_parameter.debug
    delay_output = Gfb_Delay_process(synthesizer_aural.delay, electrodogramm);
    plotChannels(delay_output,auralisation_parameter.voc_sampling_frequency_hz,1, 'r');
    title('After synthesis delay');
    % calculate and plot frequency response of auralisation filters:
    impulse = zeros(1*auralisation_parameter.voc_sampling_frequency_hz,1); impulse(1) = 1; %creates direc impulse
    % Create analyse -Filterbank for auralisation
    analyzer_aural = Gfb_Analyzer_new(auralisation_parameter.voc_sampling_frequency_hz,...
        auralisation_parameter.gamma_order_auralisation,...
        auralisation_parameter.center_frequencies_hz_auralisation,...
        auralisation_parameter.bandwidth_factor);
    % Next Analyse-Filterbank (modulates puls-pattern to desired frequency region on basiliar
    % membrane)
    impulse_out = Gfb_Analyzer_process(analyzer_aural, impulse);
    frequency_response = fft(real(impulse_out)');
    frequency = [0:length(impulse)-1] * auralisation_parameter.voc_sampling_frequency_hz / length(impulse);
    Gfb_plot(160, [0, auralisation_parameter.voc_sampling_frequency_hz/2, -40, 0], ...
        strcat('frequency response of the auralisation filterbank with spatial spread constant lambda [m] = ', num2str(auralisation_parameter.lambda)), ...
        'frequency / Hz', 'filter response / dB', ...
        frequency, 20 * log10(abs(frequency_response)));
    hold on; plot((1:length(frequency)),-3.*ones(size(frequency)),'k-');
    disp(' ');
    disp('Figure 160 shows the frequency response of the auralisation filterbank.');
    % Also add the spatial spread to the figure:
    plot(repmat(auralisation_parameter.center_frequencies_hz_auralisation,size(electrodogramm,1),1)',20*log10(mspatial_weights)','x--','Markersize',15);
end
outsig = Gfb_Synthesizer_process(synthesizer_aural, electrodogramm);
outsig = outsig'; %Make mono channel
% Resample to original sampling frequency
outsig = resample(outsig,fs,auralisation_parameter.voc_sampling_frequency_hz);
end
%%
function parameter = setParameter(signal,fs,coding_strat,varargin)
% function parameter = setParameter(signal,fs,coding_strat,varargin)
% Input parameters:
%       signal = one channel vector of the signal you want
%       to put through the CI simulation and auralization
%       fs = sampling frequency of the input signal in hz.
%       varargin = Key-Value pairs of optional parameters, for all
%       optional parameters, see below.
%
% Output parameters:
%	parameter = struct containing all the parameters needed for the
%	auralization (e.g. number of channels of the electrodogramm,
%   values for compression

% Input parser (make all parameters except the first three optional)
p = inputParser; %creates Input parser structure
p.addRequired('signal', @(x) isnumeric(x) && size(signal,2) == 1);
p.addRequired('fs', @isnumeric);
p.addRequired('coding_strat', @ischar); % Can be any acceptable switch-case from coding_strategy()
% PPS-Rate over all channels. To obtain maximum PPS per Channel divide this number by
% number of center_frequencies_hz_stimulation
p.addParamValue('voc_sampling_frequency_hz', 48000, @isnumeric);
% Bandwidth of Filters
p.addParamValue('bandwidth_factor', 3.*[1 1 1 1 1 1 1 1 1 1 1 1], @isnumeric); %ERB-Bandwidth of the analysis/synthesis filters. If one electrode needs to be disabled, change the bandwidth accordingly
% Order of Gammatonefilterbank (stimulation and auralisation analysis
% filterbank)
p.addParamValue('gamma_order_stimulation', 3, @isnumeric);
p.addParamValue('gamma_order_auralisation', 3, @isnumeric);
% Middle frequency of analysis-Filterbank
p.addParamValue('center_frequencies_hz_stimulation', [120 235 384 579 836 1175 1624 2222 3019 4084 5507 7410], @isnumeric);
p.addParamValue('pulselength', 40e-6, @isnumeric); % Default pulse-length in seconds.
p.addParamValue('ipg', 2.1e-6, @isnumeric); % Default inter-pulse-gap in seconds.
% Higher version with 12 Electrodes using greenwoods equation and an
% electrode insertion depth of 24 mm in a 32 mm cochlear duct (Flex24 Electrode)
% old auralization frequencies *with* included dominance region
% between 500 Hz and 1.2 kHz [390 550 759 1031 1384 1843 2440 3216 4225 5537 7243 9460]
p.addParamValue('center_frequencies_hz_auralisation', [120 235 384 1794 2181 2646 3203 3872 4675 5638 6793 8179], @isnumeric);
p.addParamValue('B', 0.0156, @isnumeric); %Base level of compression (see Diss Tamas Harzcos, p. 18 (4/256)
p.addParamValue('M', 0.5859, @isnumeric); %Saturation level of compression (see Diss Tamas Harzcos, p. 18 (150/256)
p.addParamValue('Volume', 1, @isnumeric); %Volume control of Med-El-CIs [0 1], to be applied in the CU-conversion step, not directly needed for Vocoder-auralisation
p.addParamValue('alpha_c', 416.2, @isnumeric); %controls the steepness of the compression function (see Diss Tamas Harzcos, p. 18
p.addParamValue('TCL', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]', @isnumeric);
p.addParamValue('MCL', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]', @isnumeric);
p.addParamValue('max_elecs_affected', 5, @isnumeric); %After how many electrodes the spatial spread should be neglected (1% of original amplitude remaining at this electrode)
p.addParamValue('lambda', 0.0036, @isnumeric); % exponential decay constant for spatial spread
p.addParamValue('distance_electrodes', 0.0028, @isnumeric); %Distance between adjunct electrodes [m] for calculation of spatial spread
p.addParamValue('pps', 800, @isnumeric);
p.addParamValue('electrodeselmethod', 'random', @ischar); % or 'sequential' (fixed order)
p.addParamValue('debug',0, @isnumeric);
p.addParamValue('weights',[1 1 1 1 1 1 1 1 1 1 1 1]', @isnumeric)
% Validation
p.parse(signal, fs, coding_strat, varargin{:});
parameter = p.Results;
% One additional validation:
assert(isequal(length(parameter.bandwidth_factor),length(parameter.center_frequencies_hz_stimulation)),'You must specify an Gammatone ERB-bandwidth for each analyzing channel');
end
%%
function [electrodogramm, parameter] = coding_strategy(envelope, finestructure, coding_strat,parameter)
switch coding_strat
    case 'CI'
        electrodogramm = pulsatile_sampling_CIS(envelope,parameter.len_pulse, parameter.block_delay, parameter.electrodeselmethod);
        electrodogramm = assertsequentialStimulation(electrodogramm);
    otherwise
        error('This Vocoder configuration is not known');
end
end
%%
function [sampled_puls_pattern] = pulsatile_sampling_CIS(signal_to_sample,pulselength_samples, pps_blockshift, sortmethod)
% function [sampled_puls_pattern] = pulsatile_sampling_CIS(signal_to_sample,pulselength_samples, pps_blockshift, sortmethod)
% Function to perform the CIS coding strategy on a given multichannel
% signal_to_sample. The pps per channel is coded in pps_blockshift,
% The pulselength of each biphasic pulse + inter-pulse-gap is coded in
% pulselength_samples. Sortmethod determines, whether the channel timing within
% one pps-blockshift is 'sequential' as is needed for
% CI stimulation or 'random', which is needed for later auralization.
%
% Input parameters:
%       signal_to_sample = multi-channel signal (channel x samples) you want to sample.
%       Usually, this is the envelope.
%       pulselength_samples = pulse-duration (anodic+cathodic+inter-pulse-gap) in samples
%       pps_blockshift = time-shift after each stimulation round across
%       all electrodes in samples. Corresponds to 1/pps*nrchannels,
%       e.g. total stimulation rate
%       sortmethod = string, either 'random' (for later auralization),
%       or 'sequential' (for CI stimulation), see
% Output parameters:
%	parameter = struct containing all the parameters needed for the
%	auralization (e.g. number of channels of the electrodogramm,
%   values for compression

% Create pulsatile pattern with correct sampling pulse and randomized sequence between electrodes
sampled_puls_pattern = zeros(size(signal_to_sample));
pulse = ones(1,pulselength_samples); % Generate pulse. Note that this pulselength contains negative Phase, inter-phase-gap and positive phase!
n_samples = 1;
min_value = 0;
nrchans = size(signal_to_sample,1);
% now estimate how much samples needs to be shifted within one
% pps-shift, e.g. going from channel 1 to channel 2 at the total
% stimulation rate (taking into acount the total pulselength in samples)
if length(pps_blockshift) > 1
    for ii = 1:length(pps_blockshift)
        channelshift(ii) = round((pps_blockshift(ii)-pulselength_samples*nrchans)/nrchans);
    end
else
    channelshift = round((pps_blockshift-pulselength_samples*nrchans)/nrchans);
end
chanshift_idx = 1;
pps_idx = 1;
if nrchans*pulselength_samples > min(pps_blockshift)
    error('Using this pps, number of channels, pulselength, and fs would result in parallel stimulation! Please adjust these three factors, so that the followind equation is satisfied, with x = [1,2,3,4,...,100]: ((x * M) + pulselength_samples)*pps*M/M = fs');
end

size_signal = size(signal_to_sample,2);
while n_samples <= size_signal-1;
    rand_channel = sortElectrodesperBlock([1:size(signal_to_sample,1)], sortmethod);
    for n_kanaele = 1:nrchans;
        %Break before last sample
        if n_samples > length(signal_to_sample)-pulselength_samples;
            break
        end
        % Pulsatile sampling if amplitude of finestructure > Threshold (currently 0)
        if signal_to_sample(rand_channel(n_kanaele), n_samples) > min_value;
            sampled_puls_pattern(rand_channel(n_kanaele),[n_samples : n_samples+pulselength_samples-1]) = ...
                signal_to_sample(rand_channel(n_kanaele),n_samples).*pulse; %sampling with puls
        end
        if signal_to_sample(rand_channel(n_kanaele), n_samples) > min_value && n_kanaele < size(signal_to_sample,1)
            n_samples = n_samples+pulselength_samples+channelshift(chanshift_idx); % Shift to next channel (sequential stimulation)
            if chanshift_idx == length(chanshift_idx)
                chanshift_idx = 0;
            end
            chanshift_idx = chanshift_idx + 1;
        end
        if n_samples == length(signal_to_sample)-1;
            n_samples = n_samples+1;
            break
        end
    end
    n_samples = n_samples+pps_blockshift(pps_idx)-(nrchans-1)*(pulselength_samples+channelshift(chanshift_idx)); %Forwards to next block (according to pps)
    if pps_idx == length(pps_blockshift)
        pps_idx = 0;
    end
    pps_idx = pps_idx +1;
end

end
%%
function [melectrodenumbers] = sortElectrodesperBlock(electrodenrs,choosemethod)
switch choosemethod
    case 'random'
        rand_idx = randperm(length(electrodenrs));
        melectrodenumbers = electrodenrs(rand_idx);
    case 'sequential'
        melectrodenumbers = sort(electrodenrs,'descend');
    otherwise
        error('You must select either ''random'' or ''sequential'' for ordering the electrodes per block!');
end
end
%%
function electrodogramm = assertsequentialStimulation(electrodogramm)
% create temporal matrix with just zeros (no stimulus) or one (stimulus)
stimpattern = zeros(size(electrodogramm));
stimpattern(electrodogramm>0) = 1;
% Calculate sum:
nr_chans_stimulated = sum(stimpattern,1);
if any(nr_chans_stimulated>1)
    warning('Parallel electrode stimulation detected! Will only keep the electrode with maximum amplitude!');
end
[~,tbins_stimulation] = find(nr_chans_stimulated>1);
for dd = 1:length(electrodogramm(nr_chans_stimulated>1))
    slice = electrodogramm(:,tbins_stimulation(dd));
    [max_val, max_idx] = max(slice);
    slice(:) = 0; %Set all channels which are stimulated at the same time to zero
    slice(max_idx) = max_val; %keep only the maximum one
    electrodogramm(:,tbins_stimulation(dd)) = slice;
end
end
%%
function env_signal_filtered = lowpass_filter_env(env_signal_new, center_frequencies,fs)
env_signal_filtered = nan(size(env_signal_new));
for ii = 1:size(env_signal_filtered,1)
    % design 1rd order low-pass butterworth IIR-Filter
    % Envelope filter coeffs with 200 Hz cutoff-frequency
    [b_env, a_env] = butter(1, (200./(fs/2)));
    env_signal_filtered(ii,:) = filtfilt(b_env, a_env, env_signal_new(ii,:));
end
end
%%
function [delay_block] = calculate_delay(fs,pps)
% this function calculates the delay in samples needed for
% "blockprocessing" at a certain total stimulation rate, which is m*pps.
% If the total stimulation rate (pps*m) is not a even divisor of fs, a
% rounding procedure will apply: The delay in samples will be a vector,
% where each entry either slightly either over or undersamples the
% correct delay. When taking the mean, the correct delay in samples is
% guranteed to be in a range of 1e5.
correct_delay = fs/(pps);
if correct_delay < 1
    error('Too high pps! This value is not supported!');
end
delay_block = ceil(fs/(pps));
while round(mean(delay_block)*1e5) ~= round(correct_delay*1e5)
    difference = correct_delay - mean(delay_block);
    if difference > 0
        delay_block = [delay_block ceil(fs/pps)];
    else
        delay_block = [delay_block floor(fs/pps)];
    end
end
if length(delay_block) > 1
    warning('pps * M is not a divisor of fs. Make sure that the following equation with x = [1,2,3,4,...,100 is satisfied: x*pps*M = fs');
end
end
%%
function [delayed_signal, signal,output_signal] = test_gammatonefilterbankdelay(signal,CISIM_parameter,auralisation_parameter)
% Create analyse -Filterbank
analyzer_stim = Gfb_Analyzer_new(CISIM_parameter.voc_sampling_frequency_hz,...
    CISIM_parameter.gamma_order_stimulation,...
    CISIM_parameter.center_frequencies_hz_stimulation,...
    CISIM_parameter.bandwidth_factor);
% First analysis-filterbank
%1st analyse:
[signal, analyzer_stim] = Gfb_Analyzer_process(analyzer_stim,signal);
frequency_response = fft(real(signal)');
frequency = [0:length(signal)-1] * auralisation_parameter.voc_sampling_frequency_hz / length(signal);
Gfb_plot(80, [0, auralisation_parameter.voc_sampling_frequency_hz/2, -40, 0], ...
    'frequency response of the analysis filterbank', ...
    'frequency / Hz', 'filter response / dB', ...
    frequency, 20 * log10(abs(frequency_response)));
hold on; plot((1:length(frequency)),-3.*ones(size(frequency)),'k-');
disp(' ');
disp('Figure 80 shows the frequency response of the analysis filterbank.');
% Create analyse -Filterbank for auralisation
analyzer_aural = Gfb_Analyzer_new(auralisation_parameter.voc_sampling_frequency_hz,...
    auralisation_parameter.gamma_order_auralisation,...
    auralisation_parameter.center_frequencies_hz_auralisation,...
    auralisation_parameter.bandwidth_factor);
% Next Analyse-Filterbank (modulates puls-pattern to desired frequency region on basiliar
% membrane)
signal = Gfb_Analyzer_process(analyzer_aural, signal);

% Signal-Synthesis:
synthesizer_aural = Gfb_Synthesizer_new (analyzer_aural, 1/100); %1/100 = Delay
delayed_signal = Gfb_Delay_process(synthesizer_aural.delay, signal);
% Process signal completely
output_signal = Gfb_Synthesizer_process(synthesizer_aural, signal);
end

%EOF
