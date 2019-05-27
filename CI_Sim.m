function [ output_sig ] = CI_Sim( input_signal,fs,CIType)
%function [ output_sig ] = CI_Sim( input_signal,fs,CIType )
% This function simulates a Cochlear CI signal processing and will then
% auralize it with a vocoder. The vocoder mimics the electrode-to-place
% mismatch of the contour advance electrode or the med el flex 24 electrode.
%
% Input parameters:
%       input_signal = input signal in samples, one channel only
%       fs = sampling frequency in hz
%       CIType = String, for the CI types. Currently only 'Cochlear' and 'Med-El' are
%       supported
%
% Output parameters:
%       output_sig = Vocoded (e.g. CI simulating signal processing and
%       auralization afterwards) output signal. It has the same sampling
%       frequency and the same (broadband) RMS as the input signal.
%
% Calibration of the CI simulation:
% CIs have a very narrow input dynamic range of about 40 dB, see parameters
% par.B and par.M or params_CI.B and params_CI.M. By default, every signal,
% which gets passed to this function will be scaled to fit into this narrow
% dynamic range. If you experiment requires level differences between
% signals to be preserved (e.g. loudness scaling or ILDs), then set the
% variable gain_to_scale to the value you get for your reference signal
% (e.g. middle loud or 0 degree reference for spatial acoustic). All other input signals will then be
% scaled by this signal, thus preserving level differences between input signals
%
% % Usage example:
% addpath('Vocoder');
% addpath('CI_CODING_STRAT');
% [signal,fs] = audioread(['Vocoder' filesep 'OLSA.wav']);
% 
% % ACE coding strategy (used in Cochlear)
% vocoded_signal = CI_Sim(signal(:,1),fs,'Cochlear');
% sound(vocoded_signal,fs);
% % CIS coding strategy (used in Med-El)
% voc_signal = CI_Sim(signal(:,1),fs,'Med-El');
% sound(voc_signal,fs)
%
% Authors: Ben Williges <ben.williges AT uni-oldenburg.de>

	% load CI function library and get common CI parameters
           CI = createCI();
           CI_Params; %get parameters
        switch CIType
        case 'Cochlear'
           %Ci Cochlear processing left ear
           par=ACEParameters;
           %% Set input signal to correct CI calibration level:
           CalibLevel=10; % Base level (par.B) + x dB [dB] (dynamic range = 40 dB)
           DesiredLevel=20*log10(par.B)+CalibLevel;
           [input_signal_CIsim, gain_to_scale] = setSignaltoRMS(input_signal,DesiredLevel,'dB');
%          gain_to_scale = 3.8625; % CAUTION: This is the calibration
%          variable. Only use this if you know, what you are doing and need
%          to e.g. represent ILDs faitfully or want to do loudness scaling
%          with CI
           input_signal_CIsim = scaleSignalbyRMS(input_signal,gain_to_scale,'dB');
           %% Run CI soding strategy
           [electrodogram, vDistance, StimOrder,vTime, mEnv]=CIStrat(input_signal_CIsim,par,fs);
           %% Auralization of the electrodogramm
           params_CI = CI.setParameter(input_signal,fs,'CI',params_cochlear{:});
           % Set some auralization parameters derived from the CI coding strategy
           params_CI.rms_per_channel = rms2(mEnv,2);
           params_CI.TCL = par.T;
           params_CI.MCL = par.C;
           params_CI.Volume = par.vol;
           TSR = unique(par.pps)*par.n; %Calculate Total stimulation rate = sampling frequency of electrodogramm signal
           electrodogram_upsampled = resample_elec_zeros(electrodogram,TSR,params_CI.voc_sampling_frequency_hz); %Upsampling is needed, because some of the electrodes are located at frequencies > TSR 
           % Perform the auralization:
           [output_sig] = CI.Auralisation(electrodogram_upsampled,fs,params_CI);
        case 'Med-El'
            params_CI = CI.setParameter(input_signal,fs,'CI',params_med_el{:});
            CalibLevel=35; % Base level (par.B) + x dB [dB] (dynamic range = 40 dB). Note that the fitting into the dynamic range only happens after the electrodogramm. There it corresponds to a level which is 10 dB above params_CI.B.
           DesiredLevel=20*log10(params_CI.B)+CalibLevel;
           [input_signal_CIsim, gain_to_scale] = setSignaltoRMS(input_signal,DesiredLevel,'dB');
%          gain_to_scale = 28.8625; % CAUTION: This is the calibration
%          variable. Only use this if you know, what you are doing and need
%          to e.g. represent ILDs faitfully or want to do loudness scaling
%          with CI
           input_signal_CIsim = scaleSignalbyRMS(input_signal,gain_to_scale,'dB');
            [electrodogramm, params_CI] = CI.Simulation(input_signal_CIsim,fs,'CI',params_CI);
            params_CI.bandwidth_factor = ones(size(params_CI.bandwidth_factor)); %Set bandwidth for auralization gammatonefilters to 1.
            [output_sig] = CI.Auralisation(electrodogramm,fs,params_CI);
            %Account for delay
        otherwise
        error(['This CI Type is not know! Allowed ''Cochlear'' or "Med-El". You entered: ' CIType]);
        end
        % set output signal to same rms as input signal
        output_sig = setSignaltoRMS(output_sig,rms2(input_signal,1),'linear');
end

