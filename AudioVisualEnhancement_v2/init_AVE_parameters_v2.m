% Init file. Initialize general parameters for the quest and main
% experiment in AVE using modified stimuli function (Alexis Perez Bellido 1/2/2017)
%run('P:\3018029.02\AudioVisualEnhancement_pedestal\pathdef.m')
%path(ans)

% run the experiment demo detect_AVE_pedestal_v2('demo')  
% run the experiment  detect_AVE('task')
%% Monitor/stimulus size calibration (ADJUST PARAMETERS FOR EACH EXPERIMENT)
%init.screenwidhsize = 33; % in cm Mac 13 = 33
%init.screen_width_resolution = 2560; % in Mac 13 2560 pixels
%init.distance2screen = 40; % in cm

init.screenwidhsize = 53; % in cm Mac 13 = 33
init.screen_width_resolution = 1920; % in Mac 13 2560 pixels
init.distance2screen = 56; % in cm

init.pixel_size = init.screenwidhsize/init.screen_width_resolution; %Pixels are squared, so we only need to know one dimension
init.one_degree_length_incm = tan(deg2rad(1/2))*init.distance2screen*2;
init.one_degree_length_in_px = init.one_degree_length_incm/init.pixel_size;
%init.one_degree_length_in_px = 2*init.distance2screen*tan((1/2)*(pi/180))*(init.screen_width_resolution/init.screenwidhsize);

% add convolve2 function to path
% addpath('./convolve2');
% set basic noisy Gabor parameters (define in each trial contrast and phase
% during the loop(

ppd = round(init.one_degree_length_in_px); % pixels per degree of visual angle

init.orientations = [90 0]; %[315 45 0 0]; % degrees of possible gabor pathces (3 and 4 positions are noise)
%detect.pre_cue = [90 90 315 45]; %orientation to be reported
%detect.post_cue = [315 45 315 45]; %orientation to be reported


%% Experiment imings
init.INIT_DELAY = 1;
%THESE TWO VARIABLES BELOW SHOULD BE LARGER THAN THE MINOR AND MAX SOA
%(otherwise u wont hear the sound sometimes)
init.FIXATION = 1.0;
init.PRE_STIM_CUE = 0.75; % from the begining of trial till the presentation of the visual stimuli
init.POST_STIM_CUE = 0.75; % This will be used later to indicate when the subject can start responding (300 ms ensures that the sound has been already presented
init.RESPONSETIME = 1.5; % time to respond and feedback
init.ITI = 1.5; % BLANK SCREEN time until the nextrial starts IN THE fMRI EXPERIMENT THE TIMING BETWEEN TRIALS WILL BE DETERMINED BY OPTSEQ
init.jitter_range = 0.25;
 
% The STIM AND CUE VALUES ARE NOT NEEDED TO DEFINE THE TRIALS TIMINGS
init.STIMDUR = 0.016 * 2; % 1000/60 hz * 5 = 60 ms
init.ASTIMDUR = init.STIMDUR; % 1000/60 hz * 5 = 60 ms
init.CUEDUR = 0.016 * 10; % after the cue presentation, the pedestal is presented. The shorter the cue is presented the longer the pedestal shows..
init.fs = 44100 %sampling fq


%if(mod(init.ASTIMDUR*fs,1) == 0) % the duration should be 
%serror('Program exit. Stimulus duration is not an integer.')    
%end
%%

% set noisy Gabor parameters
init.cfg          = [];
init.cfg.ppd = ppd; % saving ppd in the experiment
init.cfg.patchsiz = ppd*22.0; % patch size (pix) give degrees 1 degree are ppd pixels
init.cfg.patchenv = ppd*1.5; % patch spatial envelope s.d. (pix)
init.cfg.envscale = 10^10; % patch spatial envelope s.d. (pix)
init.cfg.patchlum = 0.5; % patch ba ckground luminance
init.cfg.freq = 0.5; % spatial frequency (above 5 and adding noise the creation texture process seems to be very slow)

init.cfg.gaborper = ppd*(1/init.cfg.freq); % Gabor spatial period (pix) (Transformation from SF)
init.cfg.noisedim = init.cfg.gaborper/6; % noise dimension (pix)
init.cfg.noisecon = 0.3/3; % noise RMS contrast
%init.cfg.diameter = init.cfg.patchsiz; % This way scale the gaussian evelop and I make the border very abrupt
%init.cfg.gaborang = tilt; % Gabor orientation (rad)


