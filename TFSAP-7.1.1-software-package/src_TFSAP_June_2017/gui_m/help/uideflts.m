% uideflts	Default	values for all uicontrols

% For items selected by	popup-menus, default choice is selected	by giving
% the position of that choice as it appears in the menu.  For example,
% number 1 will	give a default corresponding to	the first menu item.  Note
% that these numbers should NOT	appear in quotes, while	editable text
% defaults MUST	appear in quotes.
%

%TFSA 7.0
% Copyright
% Signal Processing Research Centre
% Queensland University	of Technology
% GPO Box 2434
% Brisbane 4001
% Australia
%
% email: tfsa@qut.edu.au

MAIN_WIN_NAME = 'Welcome to TFSA';
TFSA_VER = { 'TFSAP 7.0' };
TFSA_STR = { 'Time Frequency Signal Analysis' };
AUTHOR = { 'Developed by Professor Boualem Boashash' };
MAINTAINER = {' '};%{'Maintained by Dr. Samir Ouelha'}; 
COPYRIGHT = { 'Copyright 1987-2016' };
WWW_SPR = { 'http://www.bee.qut.edu.au/projects/spr/' };
%WWW_SPR = { '' };
WWW_SPR = MAINTAINER;



%% Set the path to find help files..
ppath = fileparts( which( mfilename ) );

HELP_FILES_PATH = fullfile( ppath, 'help' );
HTML_DOCS = fullfile( ppath, 'htmldocs', 'index.html' );
COPY_HELP = fullfile( HELP_FILES_PATH, 'copyright_tfsa.m' );
WELCOME_HELP = fullfile( HELP_FILES_PATH, 'welcome_tfsa.m' );

BL_TRANS_TITLE = 'Quadratic Transforms';% Previous BILINEAR
ML_TRANS_TITLE = 'Multilinear Transforms';

WINDOW_LIST = ['Rectangular|Hann|Hamming|Bartlett'];
WINDOW_LISTTXT = 'Smoothing Window:';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BILINEAR TRANSFORMATIONS...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BTFD_IN1_TXT = 'INPUT Time Domain Signal:';
BTFD_IN2_TXT = 'INPUT Time Domain Signal 2:';
BTFD_OUT_TXT = 'OUTPUT time-frequency distribution:';
BTFD_DIST_TXT = 'Distribution:';
BTFD_TIME_TXT = 'Time Resolution:';
BTFD_LAG_TXT = 'Lag Resolution:';
BTFD_FFT_TXT = 'FFT Length (zero pad to):';
BTFD_DATAWIN_TXT = 'Lag Window Length (odd):';

BTFD_DIST = [ 'Wigner-Ville|Wigner Distribution|Spectrogram|'...
              'Ambiguity Function|'...
              'Smoothed Wigner-Ville|S-Method|' ...
              'Rihaczek-Margenau-Hill|' ...
              'Exponential|Born-Jordan|Zhao-Atlas-Marks|' ...
              'Cross Wigner-Ville|' ...
              'B-distribution|Modified B-distribution|' ...
              'Ext Modified B-distribution|'...
              'CKD|'...
              'MDD'];

BTFD_IN1_DEF = 'time1';					% input	1 var name
BTFD_IN2_DEF = 'time2';					% input	2 var name
BTFD_OUT1_DEF =	'tfd1';  				% output 1 var name
BTFD_DTYPE_DEF = 1;					% distribution type: 1 = wvd, 2	= smooth, etc
BTFD_WSIZE_DEF = '127';					% window size (odd)
BTFD_TRES_DEF =	'1';					% time resolution

BTFD_WPARAM_DEF = '21';					% smoothing window length (odd)
BTFD_FFT_DEF = '128';				        % fft length
BTFD_WTYPE_DEF = 2;					% window type for smoothed wvd:	1 = rect, etc
BTFD_BPARAM_TXT = 'Parameter Beta:';
BTFD_MBPARAM_TXT = 'Parameter Alpha:';

%%%%%

BTFD_B_PARAM = '0.1';                                   % B-dist. parameter
BTFD_MB_PARAM = '0.1';                                  % modified B-dis. param.d


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add by Samir Ouelha %%%%%%%%%%%%%%%%%
% STFT with overlapping --- Parameter Default Values
BTFD_STFT_OVERLAP_PARAM = '60';                                  % New STFT overlap parameter 
BTFD_STFT_OVERLAP_PARAM_TXT = 'Window overlap:';

BTFD_STFT_WINDOW_LENGTH_PARAM = '63';                                  % New STFT overlap parameter 
BTFD_WINDOW_LENGTH_STFT_PARAMTXT = 'Window length:';

BTFD_PARAM_S_METHOD = '3';
BTFD_PARAM_TXT_S_METHOD = 'L';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMBD --- Parameter Default Values

BTFD_PARAM_TXT_EMB_Alpha = 'Parameter Alpha:';
BTFD_PARAM_EMB_Alpha = '0.1';                                  % Extended M B-dis. param. Alpha 

BTFD_PARAM_TXT_EMB_Beta = 'Parameter Beta:';
BTFD_PARAM_EMB_Beta = '0.1';                                  % Extended M B-dis. param. Beta

%%%%%%%%%%

%%%CKD

BTFD_PARAM_TXT_ECSK_C = 'Parameter C:';
BTFD_PARAM_ECSK_C = '0.1'; 

BTFD_PARAM_TXT_ECSK_D = 'Parameter D:';
BTFD_PARAM_ECSK_D = '0.1'; 

BTFD_PARAM_TXT_ECSK_E = 'Parameter E:';
BTFD_PARAM_ECSK_E = '0.1'; 


%%%MDD

BTFD_PARAM_TXT_MFK_theta='Theta';
BTFD_PARAM_MFK_theta='30';

%%%%%

% 
% -------------------------------------
% 
% 




BTFD_ZAMPARAM_TXT = 'ZAM Parameter ''a'':';
BTFD_CWPARAM_TXT = 'Smoothing Parameter:';
BTFD_ZAMPARAM = '11';					
BTFD_CWPARAM = '11';					
BTFD_SPECPARAM_TXT = 'Smoothing Window Length (odd):';					
BTFD_SPECPARAM = '21';					
BTFD_HELP = fullfile( HELP_FILES_PATH, 'bilinear_tfsa.m' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MULTILINEAR TRANSFORMATIONS...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SAME VALUES FOR BTFD ...
MTFD_IN_TXT = BTFD_IN1_TXT;
MTFD_OUT_TXT = BTFD_OUT_TXT;
MTFD_TIME_TXT = BTFD_TIME_TXT;
MTFD_FFT_TXT = BTFD_FFT_TXT;
MTFD_DIST_TXT = BTFD_DIST_TXT;
MTFD_TRES_DEF = BTFD_TRES_DEF;
MTFD_FFT_DEF = BTFD_FFT_DEF;
MTFD_IN_DEF = BTFD_IN1_DEF;
MTFD_OUT_DEF = BTFD_OUT1_DEF;

MTFD_DATAWIN_TXT = 'Data Window Length';
MTFD_DIST = ['Polynomial-WVD (Order 6)|'...
             'Polynomial-WVD (Order 4)'];
MTFD_DIST_TYPE = 1;
MTFD_PARAM_TXT = 'Interpolation Degree';
MTFD_PARAM = '8';
MTFD_WSIZE_DEF = '128';					
MTFD_HELP = fullfile( HELP_FILES_PATH, 'multilinear_tfsa.m' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Direct Method
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DIR_TRANS_TITLE = 'Direct TF Method';
DIRTFD_IN1_TXT = BTFD_IN1_TXT;
DIRTFD_OUT_TXT = BTFD_OUT_TXT;
DIRTFD_DIST_TXT = BTFD_DIST_TXT;
DIRTFD_TIME_TXT = BTFD_TIME_TXT;
DIRTFD_FFT_TXT = BTFD_FFT_TXT; 
DIRTFD_DATAWIN_TXT = BTFD_DATAWIN_TXT;
DIRTFD_SPECPARAM_TXT = BTFD_SPECPARAM_TXT;
DIRTFD_SPECPARAM = BTFD_SPECPARAM;    
DIRTFD_SPECWIN_TXT = 'Smoothing Window Shape';					

DIRTFD_IN1_DEF = BTFD_IN1_DEF;					% input	1 var name
DIRTFD_OUT1_DEF = BTFD_OUT1_DEF;				% output 1 var name
DIRTFD_DTYPE_DEF = 1;					% distribution type: 1 = wvd, 2	= smooth, etc
DIRTFD_WSIZE_DEF = BTFD_WSIZE_DEF;					% window size (odd)
DIRTFD_TRES_DEF = BTFD_TRES_DEF;					% time resolution
DIRTFD_WPARAM_DEF = BTFD_WPARAM_DEF;					% smoothing window length (odd)
DIRTFD_FFT_DEF = BTFD_FFT_DEF;				        % fft length
DIRTFD_WTYPE_DEF = 2;					% window type for smoothed wvd:	1 = rect, etc
DIRTFD_HELP = fullfile( HELP_FILES_PATH, 'direct_method_tfsa.m' );
DIRTFD_DIST = ['Wigner-Ville|Short-Time Fourier Transform|' ...
            'Short Time Fourier Transform (overlap)|Rihaczek|Windowed-Rihaczek'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
IF_TRANS_TITLE = 'Instantaneous Frequency Estimation';
IF_IN1_TXT = BTFD_IN1_TXT;
IF_OUT_TXT = 'OUTPUT Estimated Frequency:';
IF_DIST_TXT = 'Algorithm:';
IF_TIME_TXT = BTFD_TIME_TXT;
IF_FFT_TXT = BTFD_FFT_TXT; 
IF_DATAWIN_TXT = BTFD_DATAWIN_TXT;
IF_SPECPARAM_TXT = BTFD_SPECPARAM_TXT;
IF_SPECPARAM = BTFD_SPECPARAM;    
IF_SPECWIN_TXT = 'Smoothing Window Shape';
IF_RADIO_TXT = 'Select Smoothing Order';
IF_ARADIO_TXT = 'Select Adaptive Algorithm:';
IF_WPD_PARAM_TXT = 'Smoothing Window Length:';
IF_WPD_PARAM_DEF = 16;
IF_APARAM1_TXT = 'Adaptive Rate Constant:';
IF_APARAM2_TXT = 'Forgetting Factor:';
IF_APARAM_DEF = 0.5;
IF_ZPARAM_TXT = 'Zero Crossing Rate Window';
IF_ZPARAM_DEF = 64;
IF_PPARAM_TXT = 'Polynomial Phase Order';
IF_PPARAM_DEF = 2;
IF_PWVD_PARAM_TXT = 'Degree of Interpolation:';
IF_PWVD_PARAM_DEF = 8;



IF_RADIO1 = '1st Order (FFD)';
IF_RADIO2 = '2nd Order (CFD)';
IF_RADIO3 = '4th Order';
IF_RADIO4 = '6th Order';
IF_ARADIO1 = 'LMS';
IF_ARADIO2 = 'RMS';

IF_IN1_DEF = BTFD_IN1_DEF;					% input	1 var name
IF_OUT_DEF = 'if1';				% output 1 var name
IF_DTYPE_DEF = 1;					% distribution type: 1 = wvd, 2	= smooth, etc
%IF_WSIZE_DEF	= '511';				% window size (odd)
IF_WSIZE_DEF = BTFD_WSIZE_DEF;					% window size (odd)
IF_TRES_DEF = BTFD_TRES_DEF;					% time resolution
IF_WPARAM_DEF = BTFD_WPARAM_DEF;					% smoothing window length (odd)
%IF_FFT_DEF =	'512';					% fft length
IF_FFT_DEF = BTFD_FFT_DEF;				        % fft length
IF_WTYPE_DEF = 2;					% window type for smoothed wvd:	1 = rect, etc
IF_HELP = fullfile( HELP_FILES_PATH, 'if_estimation_tfsa.m' );
IF_DIST = ['Finite Phase Difference|Weighted Phase Difference|' ...
           'Zero-Crossing|Adaptive|LS  Polynomial Coefficients|' ...
           'Peak of Spectrogram|Peak of WVD|Peak of Polynomial ' ...
           'WVD'];
IF_DIST_DEF=1;           



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal SYNTHESIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SYNTH_TRANS_TITLE = 'Signal Synthesis';
SYNTH_IN1_DEF =  BTFD_OUT1_DEF;					%
                                                                % input	1 var name
SYNTH_IN1_TXT = 'INPUT Time Frequency Distribution matrix:';
SYNTH_IN2_DEF = BTFD_IN1_DEF;			
SYNTH_IN2_TXT = 'OPTIONAL Original Time Domain signal:';
SYNTH_OUT1_TXT = 'OUTPUT Time Domain signal:';
SYNTH_OUT1_DEF = 'synth1';				% output 1
                                                        % var name

SYNTH_DIST = ['STFT (Inverse-DFT)|STFT (Overlap-Add)|Modified STFT| ' ...
              'Modified Spectrogram|Wigner-Ville Distribution'];
SYNTH_DTYPE_DEF = 5;
SYNTH_DIST_TXT = 'Distribution:';
SYNTH_DATAWIN_TXT = BTFD_DATAWIN_TXT;
SYNTH_WSIZE_DEF = BTFD_WSIZE_DEF;
SYNTH_WIN_TXT = BTFD_SPECPARAM_TXT;
SYNTH_WIN_DEF = BTFD_SPECPARAM;		
SYNTH_WTYPE_DEF = 1;			
SYNTH_TOL_DEF = '1';                    
SYNTH_TOL_TXT = 'Iterative Tolerance Level:';
SYNTH_CHECKBOX_TXT = 'Original Signal to Reconstruct Phase:';
SYNTH_CHECKBOX_DEF = 'Supply Signal';


SYNTH_HELP = fullfile( HELP_FILES_PATH, 'synthesis_tfsa.m' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PLOT_HELP = fullfile( HELP_FILES_PATH, 'graph_tfsa.m' );
PLOT_TYPES = ['Plot|Image|PseudoColor|Waterfall 1|Waterfall 2|' ...
    'Mesh|Surf|Contour|Tfsapl'];
PLOT_COLOURS = 'Hsv|Gray|Hot|Cool|Bone|Copper|Pink|Prism|Jet';
PLOT_SHADING = 'Flat|Interp';

PLOT_INPUT_TF = 'Time-Freq:';
PLOT_INPUT1 = 'Input Data:';
PLOT_INPUT2 = 'Time:';
PLOT_TIMETXT = 'Time Domain plot';
PLOT_FREQTXT = 'Freq. Domain plot';
PLOT_LOG = 'Log Plot';
PLOT_GRAYSCALE = 'Grayscale';
PLOT_FONTSIZE_LABEL = 'Font Size:';
PLOT_FONTSIZES = ['6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|22|24|' ...
                  '26|28|30'];
PLOT_CMAP = [ 'hsv|hot|gray|bone|copper|pink|white|flag|lines', ...
              'colorcube|vga|jet|prism|cool|autumn|spring|winter|summer'];
PLOT_CMAP_DEFAULT = 11;
PLOT_SHADING = ['flat|interp|faceted'];
PLOT_SHADING_DEFAULT = 1;
PLOT_DEFAULT_FONT_SIZE = 9;
PLOT_TFCOLOURMAP_TEXT = 'TF Colour Map:';
PLOT_SHADING_TEXT = 'TF Shading:';
PLOT_TFLINE = [ 'white|black|cyan|yellow|magenta|red|green|blue' ];
PLOT_TFBACK = PLOT_TFLINE;
PLOT_TFLINE_DEFAULT = 3;
PLOT_TFBACK_DEFAULT = 2;
PLOT_TFLINE_TEXT = 'TF Line:';
PLOT_TFBACK_TEXT = 'TF Background:';
PLOT_TLINE = PLOT_TFLINE;
PLOT_FLINE = PLOT_TFLINE;
PLOT_TBACK = PLOT_TFLINE;
PLOT_FBACK = PLOT_TFLINE;
PLOT_TLINE_DEFAULT = 3;
PLOT_TBACK_DEFAULT = 2;
PLOT_FLINE_DEFAULT = 3;
PLOT_FBACK_DEFAULT = 2;
PLOT_TLINE_TEXT = 'Time Line:';
PLOT_FLINE_TEXT = 'Freq. Line:';
PLOT_TBACK_TEXT = 'Time Background:';
PLOT_FBACK_TEXT = 'Freq. Background:';
PLOT_TGRID = 'Time Grid';
PLOT_FGRID = 'Freq. Grid';
PLOT_TGRID_DEFAULT = 0;
PLOT_FGRID_DEFAULT = 0;
PLOT_TIME = 1;
PLOT_FREQ = 1;

PLOT_ZOOM_LABEL1 = 'Time zoom min:';
PLOT_ZOOM_LABEL2 = 'Time zoom max:';
PLOT_ZOOM_LABEL3 = 'Freq. zoom min:';
PLOT_ZOOM_LABEL4 = 'Freq. zoom max:';
PLOT_ZOOM_MIN = '0';
PLOT_ZOOM_MAX = '1';
PLOT_XLABEL = 'X Axis Label:';
PLOT_YLABEL = 'Y Axis Label:';
PLOT_SAMPLING = 'Sampling Freq(Hz):';
PLOT_RES = 'Time Resolution:';
PLOT_SHADING_LABEL = 'Shading:';
PLOT_FIG_LABEL = 'Figure Title:';

PLOT_F1_DEF = 'tfd1';		% time frequency variable
PLOT_F2_DEF = 'time1';			% time variable
PLOT_PLOT_DEF =	9;			% plot type
PLOT_CMAP_DEF =	9;			% colourmap
PLOT_SHAD_DEF =	1;			% shading model
PLOT_TIT_DEF = '';			% title	of plot
PLOT_XLAB_DEF =	'Frequency (Hz)';	% XLabel
PLOT_XLAB2_DEF = 'Scale';       	% XLabel
PLOT_XLAB3_DEF = 'Wavelet Coefficient';
PLOT_YLAB_DEF =	'Time (secs)';		% YLabel
PLOT_YLAB2_DEF = 'Sample';
PLOT_RES_DEF = '1';			% time resolution
PLOT_SF_DEF = '1';			% sampling frequency
PLOT_SHRINK_TEXT = '<< Collapse';
PLOT_EXPAND_TEXT = 'Expand >>';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Signal Generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SIG_TRANS_TITLE = 'Test Signal Generation';
SIG_IP_VAR1_TXT = 'Signal Type:';
SIG_IP_VAR1_DIST = ['Signal Only|Noise Only|Noisy Signal|Analytic Gen.|Demo Signal'];
SIG_IP_VAR1_DEF = 1;

SIG_OUT_TXT = 'OUTPUT Time-Domain Signal:';
SIG_VAR2_TXT = 'INPUT Real Time Domain Signal:';
SIG_ANAL_DEF = 'z1';
SIG_ANAL2_TXT = 'OUTPUT Complex Analytic Associate:';
SIG_START_TXT = 'Start Frequency (Hz):';
SIG_STD = 'Stand. Dev.:';
SIG_START_DEF = 0.1;
SIG_END_TXT = 'End Frequency (Hz):';
SIG_END_DEF = 0.4;
SIG_SAMPLE_TXT = 'Sampling Frequency (Hz):';
SIG_SAMPLE_DEF = 1;
SIG_LENGTH_TXT = 'Signal Length:';
SIG_LENGTH_DEF = 128;
SIG_CENTRE_TXT = 'Centre Frequency (Hz):';
SIG_CENTRE_DEF = 0.25;
SIG_MOD_TXT = 'Modulation Frequency (Hz):';
SIG_MOD_DEF = 0.02;
SIG_FREQDEV_TXT = 'Frequency Deviation (Hz):';
SIG_FREQDEV2_TXT = 'Number of Steps:';
SIG_FREQDEV_DEF = 1;
SIG_FREQDEV2_DEF = 3;
SIG_SNR_TXT = 'SNR in dB:';
SIG_SNR_DEF = 0.15;

SIG_REAL_TXT = 'Real';
SIG_COMP_TXT = 'Complex';
SIG_ANAL_TXT = 'Analytic';

SIG_OUT_DEF = 'time1';				% output 1 var name

SIG_HELP = fullfile( HELP_FILES_PATH, 'sig_gen_tfsa.m' );
SIG_DIST_TXT = 'Signal Characteristics:';
SIG_DIST = ['Linear FM|Quadratic FM|Cubic FM|Stepped FM|Sinusoidal ' ...
            'FM|Hyperbolic FM'];
SIG_DIST_DEMO = ['Test Signal|Bat Data|Whale Data|Whale Data 2|EEG Data|HRV Signal1|HRV Signal2|Bird Signal|Newborn-EEG-Background-SPM2013|Newborn-EEG-Seizure-SPM2013|HRV-SPM2013|HeartSound-Normal|HeartSound-Abnormal|Adult-EEG-PLED-SPM2013|speech Signal'];
SIG_DEMO_TXT = 'Demo Signals:';

SIG_NOISE_TXT = ['Uniform Noise|Gaussian Noise'];
SIG_DIST_DEF=1;           
SIG_TYPE_DEF=1;

SIG_DEMO_BAT_L = 400;
SIG_BAT_FS = 142000;
SIG_DEMO_TEST_L = 1024;
SIG_TEST_FS = 1;
SIG_DEMO_WHALE_L = 7002;
SIG_WHALE_FS = 8000;
SIG_DEMO_WHALE2_L = 25001;
SIG_WHALE2_FS = 8000;
SIG_DEMO_EEG_L = 1024;
SIG_EEG_FS = 50;

SIG_DEMO_HRV_L = 512;

%added by Asim
SIG_DEMO_HRV2_L=256;
SIG_DEMO_Bird=9761;
SIG_HRV_FS = 2;
SIG_HRV2_FS = 2;
SIG_Bird_FS=24417;

%%%%%
SIG_DEMO_EEG_Background_SPM2013_L = 7680;
SIG_EEG_Background_SPM2013_FS = 256;

SIG_DEMO_EEG_Seizure_SPM2013_L = 7680;
SIG_EEG_Seizure_SPM2013_FS = 256;

SIG_DEMO_HRV_SPM2013_L = 720;
SIG_HRV_SPM2013_FS = 4;

SIG_DEMO_HeartSound_Normal_L = 45470;
SIG_HeartSound_Nnormal_FS = 44100;

SIG_DEMO_HeartSound_Abnormal_L = 18482;
SIG_HeartSound_Abnormal_FS = 8000;

SIG_DEMO_Adult_EEG_PLED_SPM_L = 2000;
SIG_Adult_EEG_PLED_SPM_FS = 50;

SIG_DEMO_Q1W_L = 21086;
SIG_DEMO_Q1W_FS = 8000;

SIG_STDU_DEF =	'1/sqrt(12)';	% STD of Uniform Noise
SIG_STDG_DEF =	'1';		% STD of Gaussian Noise




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TIME-SCALE ANALYSIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TS_TITLE = 'Time-Scale Analysis';
TS_IN1_TXT = BTFD_IN1_TXT;
TS_IN1_DEF = BTFD_IN1_DEF;
TS_IN2_TXT = 'INPUT Wavelet Coefficient Vector:';
TS_OUT_TXT = 'OUTPUT Time-Scale Distribution:';
TS_OUT2_TXT = 'OUTPUT Wavelet Coefficient Vector:';
TS_OUT3_TXT = 'OUTPUT Time Domain Signal:';
TS_OUT_DEF = 'ts1';			% output 1 var name

TS_DTYPE_TXT = 'Distribution:';
TS_DTYPE = ['Daubechies (4 coefficient)|Daubechies (12 ' ...
                'coefficient)|Daubechies (20 coefficient)'];
TS_DTYPE_DEF = 3;			
TS_DIR_TXT = 'Transform Direction:';
TS_DIR = ['Forward|Inverse'];
TS_DIR_DEF = 1;
TS_OUTFORM_TXT = 'Output Format:';
TS_OUTFORM_SCALE = 'Scalogram';
TS_OUTFORM_COEFF = 'Wavelet coefficients';
TS_OUTFORM_DEF = 1;

%TS_WSIZE_DEF =	'511';			% window size (must be odd)
TS_WSIZE_DEF = '127';		       % window	size (must be odd)
TS_TRES_DEF = '1';			% time resolution
TS_WTYPE_DEF = 2;			% window type


TS_HELP = fullfile( HELP_FILES_PATH, 'wavelets_tfsa.m' );






