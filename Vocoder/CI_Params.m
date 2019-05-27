% script for setting vocoder-parameters, which were used to get the
% simulated CI data with the "Cochlear Vocoder" and the "Med-EL Vocoder"

% Auralisation params for Cochlear Contour Advance
params_cochlear = {...
'voc_sampling_frequency_hz', 5*7200,...% This is the sampling frequency of the ACE-Output (n*maxima * 900 pps). The electrodogramm will be upsampled to 5*7200 to have no influence in the highest auralisation filter of fs/2.
'bandwidth_factor', ones(22,1).*1,... %ERB-Bandwidth of the analysis/synthesis filters. If one electrode needs to be disabled, change the bandwidth accordingly
'gamma_order_auralisation', 3,...
'center_frequencies_hz_stimulation', [722 797 939 1075 1248 1428 1668 1943 2282 2659 3149 3655 4209 4847 5458 6057 7018 8322 9647 11057 12388 13879],... %Dummy data, not used, only needed to pass parameter check of createCI().
'center_frequencies_hz_auralisation', [722 797 939 1075 1248 1428 1668 1943 2282 2659 3149 3655 4209 4847 5458 6057 7018 8322 9647 11057 12388 13879],... % Data from Landsberger (2015) for Nucleus Contour Advance Electrode
'max_elecs_affected', 10,... %After how many electrodes the spatial spread should be neglected (1% of original amplitude remaining at this electrode)
'lambda', 0.0036,... % exponential decay constant for spatial spread in dB/m
'distance_electrodes', 0.00099,... %Distance between adjunct electrodes [m] for calculation of spatial spread. Estimated from Landsberger (2015) Data.
'debug',0}; %Switches debug mode on (1) or off (0)

%Params for Med-el Flex24 (also contains stimulation parameters)
params_med_el = {...
'voc_sampling_frequency_hz', 48000,...% Sampling frequency of the electrodogramm
'bandwidth_factor', [1 1 1 1 1 1 1 1 1 1 1 1].*3,... %ERB-Bandwidth of the analysis/synthesis filters. If one electrode needs to be disabled, change the bandwidth accordingly
'gamma_order_stimulation', 3,...
'gamma_order_auralisation', 3,...
'center_frequencies_hz_stimulation', [120 235 384 579 836 1175 1624 2222 3019 4084 5507 7410],...% Middle frequency of analysis-Filterbank
'pulselength', 32e-6,... % Default pulse-length in seconds.
'ipg', 2.1e-6,... % Default inter-pulse-gap in seconds.
'center_frequencies_hz_auralisation', [357 548 689 968 1483 2228 3319 4670 6630 9758 12530 15374],... %Based on Data from Landsberger (Insertion angle converted to physical frequency)
'B', 0.0156,... %Base level of compression (-34 dB FS) (Laneau Ph.D-Thesis, 2005, p. 124),(Harzcos, Ph.D-Thesis, p. 18).
'M', 1.5859,... %Saturation level of compression (+4 dB FS) (Reference: See above, but for 40 dB input dynamic range (Cochlear Clinical Guidance Document, p. 16)
'Volume', 1,... %Volume control of Med-El-CIs [0 1], to be applied in the CU-conversion step, not directly needed for Vocoder-auralisation
'alpha_c', 340.83,...; %controls the steepness of the compression function. Here we use a value from Stefan Fredelakes ACE implementation. There are other values to be found in the literature, (see Diss Tamas Harzcos, p. 18. for a value of 416.2)
'TCL', [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]',...
'MCL', [800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800]',...
'max_elecs_affected', 5,... %After how many electrodes the spatial spread should be neglected (1% of original amplitude remaining at this electrode)
'lambda', 0.0036,... % exponential decay constant for spatial spread (dB/m)
'distance_electrodes', 0.0024,... %Distance between adjunct electrodes [m] for calculation of spatial spread. Estimated from Landberger (2015) Data.
'pps', 800,... % 800 pps
'electrodeselmethod', 'random',...$'random',... % or 'sequential' (fixed order)
'debug',0,... %Switches debug mode on (1) or off (0)
'weights',[0.98 0.98 0.98 0.68 0.68 0.45 0.45 0.2 0.2 0.15 0.15 0.15]'}; %weights of analysis-Filterbanks 
