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

init.screenwidhsize = 50; % in cm Mac 13 = 33 cm, in lab pc = 53 cm, in MEG = 80
init.screen_width_resolution = 1920; % retina display 1280 in pc 1920 % 
init.distance2screen = 75.0; % in cm

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
% The STIM AND CUE VALUES ARE NOT NEEDED TO DEFINE THE TRIALS TIMINGS
init.STIMDUR = 0.016 * 2; % 1000/60 hz * 5 = 60 ms
init.ASTIMDUR = init.STIMDUR; % 1000/60 hz * 5 = 60 ms
init.CUEDUR = 0.016 * 2; % flicker
init.fs = 44100 %sampling fq


%if(mod(init.ASTIMDUR*fs,1) == 0) % the duration should be 
%serror('Program exit. Stimulus duration is not an integer.')    
%end
%%

% set noisy Gabor parameters
init.cfg          = [];
init.cfg.ppd = ppd; % saving ppd in the experiment
init.cfg.patchsiz = ppd*22.0; % patch size (pix) give degrees 1 degree are ppd pixels
init.cfg.patchenv = ppd*5; % patch spatial envelope s.d. (pix)
init.cfg.envscale = 10^100; % patch spatial envelope s.d. (pix)
init.cfg.patchlum = 0.5; % patch ba ckground luminance
init.cfg.freq = 0.5; % spatial frequency (above 5 and adding noise the creation texture process seems to be very slow)

init.cfg.gaborper = ppd*(1/init.cfg.freq); % Gabor spatial period (pix) (Transformation from SF)
init.cfg.noisedim = init.cfg.gaborper/6; % noise dimension (pix)
init.cfg.noisecon = 0.3/3; % noise RMS contrast
%init.cfg.diameter = init.cfg.patchsiz; % This way scale the gaussian evelop and I make the border very abrupt
%init.cfg.gaborang = tilt; % Gabor orientation (rad)

 
