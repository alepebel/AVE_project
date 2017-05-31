function detect_AVE(task)

%% Alexis P�rez Bellido 20-5-2017
% this function run 2 possible task under MEG conditions: 1 task is a
% visual detection task, where the participants have to detect visual
% grating. 2 task is a oddball counting task, where the participants have
% to attend to rings which are black.
% Stimulation can be visual, auditory and audiovisual.

%% Triggers codes using BITSI
% X0 Block number X (e.g. 10)
% X1 Fixation point in block number X (e.g. 11, 21, 31 for 1st, 2nd and 3d block respectively)
% 1X1 1X2 1sr Pedestal (e.g. 111 or 112 center or periphery respectively)
% X2 X3 1st Cue presentation (e.g. 12 or 13, Cue and non cue respectively; important for oddball task)
% 1X3 1X4 Target (e.g. 113 or 114 signal or noise respectively) In noise events actually nothing happens
% X4 X5 2nd Cue presentation (e.g. 14 or 15,  Cue and non cue respectively; important for oddball task)
% X6 Response phase starts (e.g. 16,)
% X7 X8 Sound presentation (e.g. 17,18, audivisual and visual trials respectively)
% 1 2 3 Incorrect, Correct, NoResponse


% to calibrate time use this as 1 and it will draw very big patches on the
% screen that can be measure with an oscilloscope
time_calibration = 0;

if strcmp(task,'demo')
    uiwait(msgbox(['Hello Mr/Ms Experimenter.. Remember to adjust the sounds intensity. The speaker wheel mark should be vertical and windows volume = 20 units if you uses the SPEAKERS. WINDOWS VOLUME SHOULD BE 4 IF YOU USE THE HEADPHONES. Adjust and proceed...'],'modal'));
    disp('Remember to adjust the experiment sounds volume. Speaker wheel mark should be vertical and windows volume 5 units')
end
% Main experiment. Use contrast from quest to adjust stimulus to threshold
% level. This function precreates gabor patches with particular
% characteristics and calculate their energy for each orientation.
% Provide feedback Use detect_AVE('demo') for demo or detect_AVE('Exp') for
% experiment

%clear;
tic

%% Set current directory
%path = '/Volumes/ALEX_EXT/Experimentos/DondersExperiments' %
path = mfilename('fullpath'); % Finds current location of file
[pathstr] = fileparts(path); % Removes file information, saves directory
path = pathstr;

%% Collect subject info

% universal file separator across operative systems is filesep

if strcmp(task, 'demo')
    prompt = {'Contrast center (0.007):','Contrast surround (0.007):'};
    input = inputdlg(prompt);
    initials = 'demo'; % Collect subject initials
    %session = input{1};
    %session  = str2num(session);
    session = 1;
    contrast_center = input{1};
    contrast_center  = str2num(contrast_center);
    contrast_surround = input{2};
    contrast_surround  = str2num(contrast_surround);
    
else
    
    prompt = {'Enter subject initials: ','Contrast center (0.007):','Contrast surround (0.007):'};
    input = inputdlg(prompt);
    initials = input{1}; % Collect subject initials REMEMBER ATACCH ALL THE
    contrast_center = input{2};
    contrast_center  = str2num(contrast_center);
    contrast_surround = input{3};
    contrast_surround  = str2num(contrast_surround);
end


outputFolder = [path, filesep 'results' filesep ,initials,filesep];

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

dataFile = [outputFolder initials,'_AVE_detect.mat'];

if exist(dataFile, 'file')
    disp('This subject has been already tested, checking session');
    load(dataFile);
    n_sessions_runned = length(Experimentalsession);
    uiwait(msgbox(['Session tested previously = ', num2str(n_sessions_runned)],'modal'));
end

if  ~strcmp(task, 'demo')
    % SELECTING SESSION AND CHECKING THAT HAS NOT BEEN TESTED ALREADY
    prompt = {'Enter session: '};
    input = inputdlg(prompt);
    session = input{1};
    session = str2num(session);
    choice = [];
    
    if exist(dataFile, 'file') %
        
        if sum([1:n_sessions_runned] == session)
            choice = questdlg('Session tested previously. Do you want to overwrite?', ...
                'Yes', 'No');
        end
        
        if strcmp(choice,'No')
            clear;
            % disp('Avoiding overwritting data. Session stops')
            error('Avoiding overwritting data. Session stops');
        end
    end
end

% folder to save the raw images
imgFile = [outputFolder 'imgfiles' filesep,'/'];

if ~exist(imgFile,'dir')
    mkdir(imgFile);
end

rawimg_fname = ['IMG_run_',num2str(session),'.mat']

if session == 1 & ~strcmp(task, 'demo') % only ask for these data in the first session of the main experiment
    prompt = {'Enter Age: ','Enter Gender(M/F): '};
    input = inputdlg(prompt);
    Age = input{1};
    Gender = input{2};
    Experimentalsession(session).Age = Age;
    Experimentalsession(session).Gender = Gender;
end


Experimentalsession(session).initials = initials;
Experimentalsession(session).session = session;
Experimentalsession(session).date = datetime;


if (mod(session, 2))
    procedure = 'detection'; % or vdetection or oddball task
else
    procedure = 'oddball'; % or vdetection or oddball task
end

% first experiment is only detection
procedure = 'detection';

Experimentalsession(session).procedure = procedure;

MEG = 1; % 0 behavioral, 1 with MEG.


% creating the  subject file;
save(dataFile,'Experimentalsession');

% loading general parameters
init_AVE_parameters_pedestal_meg;

detect = init;

% bitsi stuff and trigger coding
if MEG
    detect.bitsi = Bitsi('com1');
    detect.triggers = def_detect_AVE_MEG_Triggers();
else
    detect.bitsi = Bitsi('');
    detect.triggers = def_detect_AVE_MEG_Triggers();
end




% detect.bitsi.validResponses = 'ae';
% detect.bitsi.clearResponses;
%  [response time] = detect.bitsi.getResponse(3, true)

% If MAC, select the external device
%[keyboardIndices, productNames] = GetKeyboardIndices
%deviceNumber = 2;

% PREPARING STIMULI COMBINATIONS

reps =  12;

if strcmp(task, 'demo')
    reps =  2;
    
end

%detect.orientations = [90]; % degrees of possible gabor pathces
%detect.pre_cue = [90 90 315 45]; %orientation to be reported
%detect.post_cue = [315 45 315 45]; %orientation to be reported

% in the informative trials, the orientation is determined by the
detect.stim_conditions = {'signal','noise'};
%detect.combination_conditions = {'A_irrel_CCW','A_irrel_CW','A_rel_CCW','A_rel_CW'};
detect.modality_conditions = {'V','AV'};
detect.visual_presentation = {'center','surround'};


% Combinations of conditions and randomization
allCombinations = combvec([1:length(detect.stim_conditions)],[1:length(detect.modality_conditions)],[1:length(detect.visual_presentation)]);
allCombinations = allCombinations'; % transposing all combinations matrix
sizeComb        = length(allCombinations);
nTrials         = reps*sizeComb;
allCombinations = repmat(allCombinations,reps,1);
%generic variable to randomize inside the allCombination
rand_trials_vector = 1 : nTrials;
rand_trials_vector = rand_trials_vector(randperm(nTrials));  % This will serve to select the corresponding arrow during the experiment.

%randomize order
stim_combinations = allCombinations(rand_trials_vector,:);

detect.stim_combinations = stim_combinations;
detect.repetitions = reps;
detect.nTrials = length(stim_combinations);


% creating vector with oddball trials for the attention experiment

if strcmp(procedure,'oddball')
    detect.proportion_oddball = 0.1;
else
    detect.proportion_oddball = 0.0; % no oddball in the visual detection task
end

n_oddballs = round(detect.proportion_oddball*detect.nTrials); % 10 % of the trials are oddballs in the attention task
detect.oddballs_vector = [ones(1,n_oddballs),zeros(1,detect.nTrials-n_oddballs)];
detect.oddballs_vector = detect.oddballs_vector(randperm(detect.nTrials));
detect.oddtions = [1 0; 0 1; 1 1]; % possible oddball trials (it presents the possiblities that the 1st, 2nd or both cues change to white)


%% Experiment events timings
for t_idx = 1 : detect.nTrials % consider randomizing a little bit the ITI in each trial
    
    detect.INIT_DELAY = 2 % time before the first trial
    detect.FIXATION(t_idx) = 0.75 + (1.0-0.75)*rand(1,1); % intertrial time (fixation point)
    
    detect.PRE_STIM_TEMPLATE(t_idx) = 1.0 + (1.5-1.0)*rand(1,1); % when this time finishes the target is drawn
    detect.CUE1(t_idx) = 0.5 + (0.8-0.5)*rand(1,1); % time to flicker the first visual distractor
    
    detect.POST_STIM_TEMPLATE(t_idx) = 1.25; % after this, time to responde
    detect.CUE2(t_idx) = 0.5 + (1.0-0.5)*rand(1,1); % time to flicker the second visual distractor
    
    detect.RESPONSETIME(t_idx) = 1.5;
    detect.FEEDBACK(t_idx) = 0.25;
    
    if (t_idx== 1) % previously the function did not create the timings correctly, this is my modification
        detect.TRIAL_TIMES(t_idx) = detect.INIT_DELAY;
    else
        
        detect.TRIAL_TIMES(t_idx) = detect.TRIAL_TIMES(t_idx-1) + ...
            detect.FIXATION(t_idx) + detect.PRE_STIM_TEMPLATE(t_idx)  + ...
            detect.POST_STIM_TEMPLATE(t_idx) + detect.RESPONSETIME(t_idx) + ...
            detect.FEEDBACK(t_idx); %
        
    end
    
end

% detect.INIT_DELAY   + (t_idx-1)*(detect.FIXATION(t_idx) + detect.PRE_STIM_TEMPLATE(t_idx)  +  detect.POST_STIM_TEMPLATE(t_idx) + detect.RESPONSETIME(t_idx)+detect.FEEDBACK(t_idx)); %

%% RECORDING DATA IN THIS TEMPLATE MATRIX
detect.responses_vector = zeros(detect.nTrials,20); % vector to record conditions 1 session, 2 trial, 3stim time, 4 orientation,5postcue,6intensity, 7 resp, 8 RT, 9 correct
detect.responses_vector = num2cell(detect.responses_vector);
% these will be my variables
detect.responses_vector_variables = {'session','trial','Pre_time','Stim_time','grating_tilt','contrast','SF','visual_space','keypressed','correct','RT','subject'};

[nrows,ncols] = size(detect.responses_vector); % this makes very slow the process. Record the data at the end in txt file

%% setup keyboard basic parameters

KbName('UnifyKeyNames')

%will avoid subject pressing incorrect buttons. It can also improve
%computer performance

%% Video set up %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check script is running on Psychtoolbox-3

%PsychDebugWindowConfiguration % transparent window to debug

%Screen('Preference', 'SkipSyncTests', 1); % In lab I dont need to skip syn
%tests
AssertOpenGL;

% Selects most peripheral screen
whichScreen = max(Screen('Screens'));
%rect=[0 0 800 600]; % mini-screen.
rect=[0 0 640 480]; % mini-screen.

HideCursor;

% Open graphics window on peripheral screen at certain resolution
% Our test screen in the lab will have dimensions 1024 x 768, so this is to
% test the look/feel on that size screen
[win winRect] = Screen('OpenWindow', whichScreen,[]) %,,rect wsswwws

% Running on PTB-3? Abort otherwise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Auditory stimuli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ASTIMDUR  = detect.ASTIMDUR; %sec.
fs = detect.fs; %sampling fq
% creating the whitenoise
white_noise = (rand(1,round(ASTIMDUR*fs))*2-1); % LETS create some whitenoise.
[white_noise] = linramp(white_noise,ASTIMDUR,fs,0.015); % apply ramp on/off

% create a warning beep (beeper works very bad so I create my own
% stimulus
freq = 1000; % hz
beep_noise = sin(linspace(0, ASTIMDUR *freq*2*pi, round(ASTIMDUR*fs)));
[beep_noise] = linramp(beep_noise,ASTIMDUR,fs,0.005); % apply ramp on/of
AVdelay = -0.08; % negative values sound before the visual, postive values after
detect.AVdelay = AVdelay;

InitializePsychSound;
% Perform basic initialization of the sound driver:
nrchannels = 1;

devices = PsychPortAudio('GetDevices')
audiodevice = 6;
pahandle   =  PsychPortAudio('Open',audiodevice, [], 1, fs, nrchannels);
beephandle =  PsychPortAudio('Open',audiodevice, [], 1, fs, nrchannels);

% load memory buffers with auditory stimulus. Ready to fire!.
PsychPortAudio('FillBuffer', pahandle, white_noise); % Auditory stimulus
PsychPortAudio('FillBuffer',beephandle, beep_noise); % warning signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initializing some screen basic parameters

[w, h] = RectSize(winRect);
posx = [w/2]; wMid = [w/2];    %coordinates
posy = [h/2]; hMid = [h/2];
poscenter = [w/2, h/2];

% Basic color codes
black = BlackIndex(win);
white = WhiteIndex(win);
gray = (black + white) / 2 + 1;

% Taking the absolute value of the difference between white and gray will
% help keep the grating consistent regardless of whether the CLUT color
% code for white is less or greater than the CLUT color code for black.
absoluteDifferenceBetweenWhiteAndGray = abs(white - gray);

percetange_stimuli = 0;
Creating_stimuli_screen;

% load basic noisy Gabor parameters
cfg = detect.cfg

% creating center and surround masks
patchsiz = cfg.patchsiz;
[rr cc] = meshgrid(1:patchsiz,1:patchsiz);
mask_radius = patchsiz/2;

mask(1).filter = sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)<=mask_radius/2;  %center in pixels
mask(1).filter = mask(1).filter - (sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)<=(0.75*ppd)); % add inner empty area
%imagesc(mask(1).filter)
mask(2).filter = sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)<=mask_radius;  %surround in pixels
mask(2).filter = mask(2).filter .* (sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)>=mask_radius/2);  %surround in pixels

%imagesc(mask(2).filter)
%% Visual cue size (it should be in the middle between one an another stimuli

% the area of stimulation in the center is 4 times smaller than in the
% surround (fovea see 2 mm)
%sum(sum(mask(2).filter) )
%sum(sum(mask(1).filter) )
% precompute signal energy filters
%gaborang = pi/180*(70:1:110); % filter orientations (rad) % calculate only a small part of the spectrum. Then you can calculate the rest from the images.
%gaborang = pi/180*(60:1:120); % filter orientations (rad) % calculate only a small part of the spectrum. Then you can calculate the rest from the images.
%gaborang = pi/180*(1:1:180); % filter orientations (rad) % calculate only a small part of the spectrum. Then you can calculate the rest from the images.

gaborang = pi/180*(90); % we calculate energy only for the vertical orientation
% precomputing energy filter for center and surround trials

% in principle you can save this values and use them for further
% experiments. However, at least using intranet, takes more time to save
% and load the files than create them on the fly.

%  center_filters_file = [outputFolder 'center_filters.mat'];
%  surround_filters_file = [outputFolder 'surround_filters.mat'];

%if(session == 1 & ~strcmp(task,'demo') )
if ~strcmp(task, 'demo')
    for j = 1:length(gaborang)
        % center
        cfg_center_e{j} = get_patch_energy(setfield(cfg,'gaborang',gaborang(j)));
        cfg_center_e{j}.gabor{1,1} = mask(1).filter.* cfg_center_e{j}.gabor{1,1};
        cfg_center_e{j}.gabor{1,2} = mask(1).filter.* cfg_center_e{j}.gabor{1,2};
        % surround
        cfg_surround_e{j} = get_patch_energy(setfield(cfg,'gaborang',gaborang(j)));
        cfg_surround_e{j}.gabor{1,1} = mask(2).filter.*cfg_surround_e{j}.gabor{1,1};
        cfg_surround_e{j}.gabor{1,2} = mask(2).filter.*cfg_surround_e{j}.gabor{1,2};
        % save(center_filters_file,'cfg_center_e');
        % save(surround_filters_file,'cfg_surround_e');
    end
    
end
%load(center_filters_file)
%load(surround_filters_file)

cfg_x = get_patch_contrast(cfg);

nimages = detect.nTrials; % create as many images as trials

% Pre-create test images (define in each trial contrast and phase
% during the loop)
x = nan(nimages,1); % RMS contrast
e = nan(nimages,length(gaborang)); % signal energies
hbar = waitbar(0.0,'making test images...');
set(get(findobj(hbar,'Type','Axes'),'Title'),'FontSize',16);

% Precreate all the visual stimuli
for i = 1:nimages
    
    if mod(i,5) == 0, waitbar(i/nimages,hbar); end
    target_definition = detect.stim_combinations(i,1); % selecting if it signal or noise
    target_space = detect.stim_combinations(i,3); % center surround
    
    %target_space = 1
    tilt = pi/180*detect.orientations(target_definition);
    cfg.gaborang = tilt; % Gabor orientation (rad)
    
    if strcmp(detect.stim_conditions(target_definition),'noise') % 1 is CCW, 2 is CW and 3 is noise
        cfg.gaborcon = 0; % Gabor Michelson contrast
    else
        
        % select contrast for center or surround
        if  strcmp(detect.visual_presentation(target_space),'center')
            cfg.gaborcon =   contrast_center; % Gabor Michelson contrast
        else
            cfg.gaborcon =   contrast_surround; % Gabor Michelson contrast
        end
        
    end
    
    % presenting center or surround
    if  strcmp(detect.visual_presentation(target_space),'center')
        spatial_filter = mask(1).filter;
        indexes_sf = find(spatial_filter == 0);
    else
        spatial_filter = mask(2).filter;
        indexes_sf = find(spatial_filter == 0);
    end
    
    cfg.gaborphi = rand; % Gabor unit phase (randomize in each trial)
    % clip patch luminance to [0,1] range to avoid ugly peaks of contrast in the
    % patch
    
    percetange_stimuli = (i/nimages)*100;
    Creating_stimuli_screen;
    
    while true
        
        [img noiseimg] = make_gabor_and_noise(cfg);
        %[patchimg] = absoluteDifferenceBetweenWhiteAndGray *
        %patchimg;\
        % making not stimulus part like grey background
        img(indexes_sf ) = 0.5;
        noiseimg(indexes_sf) = 0.5;
        if all(img(:) > 0 & img(:) < 1) % not larger than 1 or 0 which are the rgb limits
            break
        end
    end
    % get weighted RMS constrast
    x(i,1) = get_patch_contrast(cfg_x,img);
    % get signal energies
    if ~strcmp(task, 'demo')
        if strcmp(detect.visual_presentation(target_space),'center') %selecting center surround filters for calculating energy for each orientation in each trial
            for j = 1:length(gaborang)
                e(i,j) = get_patch_energy(cfg_center_e{j},img);
            end
        else
            for j = 1:length(gaborang)
                e(i,j) = get_patch_energy(cfg_surround_e{j},img);
            end
        end
    end
    
    img = (gray * img)/0.5;
    noises(i).noiseimg = (gray * noiseimg)/0.5;    % scale to rgb values
    stimuli(i).img = img; % dont save in the same file, it is to heavy
    detect.stimuli(i).cfg = cfg;
end

close(hbar);

% I wont save the stimulus anymore. I simply calculate the sensitivities
% online

%if ~strcmp(task, 'demo') % only save stimuli in the real Exp task
%    save([imgFile,rawimg_fname],'stimuli','-v7.3'); % if we are not
% manipulating the stimuli noise do not save the texture
%end

%saving stimuli energies
if ~strcmp(task, 'demo')
    detect.energies = e;
end
detect.contrast = x;
detect.orient_filters = 180*gaborang/pi;

% Display instructions whereas stimuli are created

Detectinstructions;

% Wait until instructions have been readed
if MEG
    detect.bitsi.clearResponses
    getResponse(detect.bitsi,1000,true); % seems that when it says false, even if you press a button it doesnt exit the function
else
    KbWait(); % press a key to continue with the script
end

if strcmp(task, 'demo')
    examplestimuli;
    
    if MEG
        detect.bitsi.clearResponses
        getResponse(detect.bitsi,1000,true); % seems that when it says false, even if you press a button it doesnt exit the function
    else
        KbWait(); % press a key to continue with the script
    end
    
    Screen('Close', signalexampletex)
    Screen('Close', noisexampletex)
end


%% Initialize quest parameters / only for the visual detection experiment
if strcmp(procedure, 'detection')
    pThreshold=0.72; % (previous was 0.75 and gamma 0.5)
    tGuessSd = 0.035;
    beta=3.5;delta=0.01;
    gamma=0.5; % gamma is the proportion correct of responses that are expected when stimulus is not present. slope, threshold and minimum value (it is a 3AFC, noise, 2 x tilt )
    % the quest is designed only to update in noise and signal trials
    % Si el sujeto no ve nada, juzgara aleatoriamente, asi que puede acertara
    % el 50% de las veces
    see_quest = 0;
    range = 0.1; % range below and above threshold to be tested
    %grain = 0.001;
    
    % lets initialize one quest for the center and another for the surround
    qcenter =                  QuestCreate(contrast_center,tGuessSd,pThreshold,beta,delta,gamma,[],range,see_quest);
    qcenter.normalizePdf =      1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.
    %qcenter_tTest =             QuestQuantile(qcenter); % assign a first tTest value
    
    qsurround =                 QuestCreate(contrast_surround,tGuessSd,pThreshold,beta,delta,gamma,[],range,see_quest);
    qsurround.normalizePdf =    1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.
    %qsurround_tTest =           QuestQuantile(qsurround); % assign a first tTest value
    
    % preparing the quest that runs across different runs
    if (session == 1)
        % If this is the first run, initialize the intersession quest
        q_center_intersession = qcenter;
        q_surround_intersession = qsurround;
    else
        %If this is not the first run, continue using the quest from the
        %previous run
        %q_intersession = q;
        if strcmp(Experimentalsession(session-1).procedure, 'oddball')
            q_center_intersession = Experimentalsession((session-2)).detect.q_center_intersession;
            q_surround_intersession = Experimentalsession((session-2)).detect.q_surround_intersession;
        else
            if(isfield(Experimentalsession((session-1)).detect, 'q_center_intersession'))
                q_center_intersession = Experimentalsession((session-1)).detect.q_center_intersession;
            else
                q_center_intersession =  qcenter;
            end
            
            if(isfield(Experimentalsession((session-1)).detect, 'q_surround_intersession'))
                q_surround_intersession = Experimentalsession((session-1)).detect.q_surround_intersession;
            else
                q_surround_intersession =  qsurround;
            end
        end
    end
    % thus q is the quest that it is applied in this run, and q_intersession is
    % the quest across sessions.
    
    detect.bitsi.clearResponses;
    % POSITION COORDINATES FOR CIRCULAR OBJECTS(FP AND INNER ANNULUS)
    
    cue_color = 150; % color of the visual oddball in grey scale
    
    color_fp = [0 0 0];
    radious_fp = 1*ppd; % size of fp in degrees
    baseRect_fp = [0 0 radious_fp  radious_fp ]; % create a fake box to set the center of ovals
    centeredRect_fp = CenterRectOnPointd(baseRect_fp, w/2, h/2);
    
    radious_center = 5*ppd;
    baseRect_center = [0 0 radious_center  radious_center ]; % create a fake box to set the center of ovals
    centeredRect_center = CenterRectOnPointd(baseRect_center, w/2, h/2);
    
    % Calibration cue. Use to calibrate the timing wuith the oscilloscope
    baseRect_calibration_cue = [0 0 500 500 ]; % create a fake box to set the center of ovals
    centeredRect_calibration_cue = CenterRectOnPointd(baseRect_calibration_cue, w/2, h/2);
    
    
    
    %% taking pictures
    imageArray = Screen('GetImage',win);
    imwrite(imageArray, [path, filesep 'results' filesep 'peripherycue1.png'])
    
    
    % Opening file to save data
    
    dataadaptoutput = [path, filesep 'results',filesep initials,'_session_',num2str(session),'_AVE_pedestal_meg_detect.txt'];
    
    %% Preparing threat... Initialize loop parameters
    iTrial = 0;
    state = 'Initialize';
    exit_loop = 0;
    
    Screen('FillRect', win, gray);    % Gray background
    vbl = Screen('Flip', win);
    ifi = Screen('GetFlipInterval', win);
    flipLag = .9 / Screen('NominalFrameRate', win);
    % use this line to draw and avoid tearing effects ->  vbl = Screen('Flip', win, vbl + 0.5 * ifi);
    
    %tstart = vbl;
    %count_frames = 0;
    
    init_time = Screen('Flip', win);
    
    detect.bitsi.sendTrigger((session*10)) % sends a trigger marking the beginning of each block
    
    while ~exit_loop % exit the loop when trials are finish
        
        % Updating time
        
        absolute_time = GetSecs;
        
        time = absolute_time - init_time; % measures relative elapsed time
        
        switch state
            
            case 'Initialize'
                %   disp('Initialize')
                
                
                if (iTrial  >= detect.nTrials)
                    exit_loop = 1; % ends the loop
                    
                else % buffer stim and get time for next stimulus presentation
                    
                    iTrial = iTrial + 1;
                    % Trial conditions
                    
                    % select center-surround combinations x modality x
                    % signal/noise
                    stimulus_definition = detect.stim_combinations(iTrial,1); % select orientation target
                    stimulus_condition = detect.stim_conditions(stimulus_definition);
                    modality_definition = detect.stim_combinations(iTrial,2); %
                    modality_condition = detect.modality_conditions(modality_definition);
                    target_space = detect.stim_combinations(iTrial,3);
                    target_definition = detect.stim_combinations(iTrial,1); % orientation definition
                    % modality_condition = 'AV';
                    
                    
                    if target_space == 1 % saving contrast
                        contrast_trial = contrast_center;
                    else
                        contrast_trial = contrast_surround;
                    end
                    
                    disp('____________________________')
                    disp(iTrial)
                    disp (['stimulus presented in the ' char(detect.visual_presentation(target_space))])
                    disp (['contrast level ' num2str( contrast_trial)]);
                    
                    % selecting whether the trials if an oddball or not
                    % (in attention experiment there are oddball trials)
                    oddball = detect.oddballs_vector(iTrial);
                    % oddball = 1;
                    
                    if oddball
                        rand_oddball = randi([1 3]);
                        oddball_cfg =   detect.oddtions(rand_oddball,:);
                    else
                        rand_oddball = 4; % if there is no oddball, name 4 the trial
                        oddball_cfg = [0 0];
                    end
                    
                    %% Here I prepare the targets texture
                    % all the target textures were pre-generatedthe at the begining). Save the matrix always
                    % for regression analyses
                    
                    % Creating texture
                    gabortex = Screen('MakeTexture', win,   stimuli(iTrial).img); %Texture generated
                    noisetex = Screen('MakeTexture', win,  noises(iTrial).noiseimg); %Texture noise generated
                    % imagesc(n)
                    texrect  = Screen('Rect', gabortex); %Extract information about texture size
                    scale    = 1; % keep same proportion...
                    
                    rotAngles = [0]; % we create the gratings already rotated, so no rotation parameter should be applied to the texture.
                    dstRects  = CenterRectOnPoint(texrect*scale , posx, posy);
                    % Inner boundary and VC
                    in_VC_Rects  = CenterRectOnPoint(texrect/2*scale , posx, posy); % VISUAL CUE: diameter in center is half than in the surround
                    out_VC_Rects  = dstRects; % VISUAL CUE: diameter in center is half than in the surround
                    
                    VC_thickness = ppd*0.2;%
                    
                    %Initializing flags / We only want the functions to draw
                    %each item one time to avoid delays
                    Vflag =     0;
                    Aflag =     0;
                    Vprecue =   0;
                    Vpostcue =  0;
                    response =  0;
                    keycode =  [];
                    correct =   0;
                    RT =        0;
                    first_cue = 0;
                    second_cue = 0;
                    resp_trigger = 0;
                    keypressed = 'NoResponse';
                    response_given = 'NoResponse';
                    state =      'Fixationpoint';
                    keyIsDown = 0;
                    triggerpressed = 0;
                    % detect.bitsi.clearResponses;  % empty bitsi buffer
                    empty_resp = 0;
                end
                
                
                %% SETTING UP TRIAL REFERENCE TIMINGS
                TRIAL_FP      = detect.TRIAL_TIMES(iTrial); % DRAWS FP
                TRIAL_1PED    = TRIAL_FP + detect.FIXATION(iTrial); % DRAW 1s PEDESTAL
                TRIAL_CUE1    = TRIAL_1PED + detect.CUE1(iTrial); % DRAW 1ST CUE
                
                TRIAL_TARGET_2PED  =  TRIAL_1PED + detect.PRE_STIM_TEMPLATE(iTrial); % DRAWS TARGET & PEDESTAL
                TRIAL_CUE2        = TRIAL_TARGET_2PED + detect.CUE2(iTrial); % DRAW TARGET
                TRIAL_RESPONSE    = TRIAL_TARGET_2PED + detect.POST_STIM_TEMPLATE(iTrial); %
                TRIAL_FEEDBACK    = TRIAL_RESPONSE + detect.RESPONSETIME(iTrial); %
                TRIAL_ENDS        = TRIAL_FEEDBACK + detect.FEEDBACK(iTrial);
                
                TRIAL_DURATION = TRIAL_ENDS - TRIAL_FP;
                
                %% Response options in the experiment
                % You have to save the response options in each trial
                
                if strcmp(procedure, 'oddball')
                    
                    if rand_oddball == 1 | rand_oddball == 2
                        possible_responses =  [ '0','1'; '1','2'];
                        response_options =  possible_responses(randi(2),:)
                    end
                    
                    if rand_oddball == 3
                        possible_responses =  [ '0', '2'; '1','2'];
                        response_options =  possible_responses(randi(2),:)
                    end
                    
                    if rand_oddball == 4
                        possible_responses =  [ '0', '1'; '0','2'];
                        response_options =  possible_responses(randi(2),:)
                    end
                    
                else
                    response_options = ['Y', 'N']; % left and right response option respectively
                    response_options = response_options(randperm(2));
                end
                
                
            case 'Fixationpoint'
                
                if (time >=   TRIAL_FP)
                    % Display fixaxtion point
                    % disp('fixation point')
                    %  disp(time)
                    
                    Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp); % new and better centered fixation point
                    
                    if time_calibration
                        Screen('FillOval', win, [255 255 255],centeredRect_calibration_cue, 400);
                    end
                    
                    detect.bitsi.sendTrigger((session*10)+1) % sends a trigger marking the beginning of fp
                    
                    Screen('Flip', win);
                    state = 'Pedestal_One';
                end
                
                
            case 'Pedestal_One'
                
                if (time >=  TRIAL_1PED) % DRAW VISUAL CUE
                    
                    if (time >=  TRIAL_CUE1 & first_cue == 0 )
                        Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                        Screen('FillOval', win,[0 0 0], centeredRect_fp, radious_fp);
                        
                        if target_space == 1 % for center oddball
                            if oddball_cfg(1) % this will only be 1 if it is an oddball detection block
                                %         Screen('FrameOval', win, [cue_color cue_color cue_color], in_VC_Rects, VC_thickness );
                                %         Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                                detect.bitsi.sendTrigger((session*10)+2) % sends a trigger marking the presence of a CUE
                                
                            else
                                %         Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                                %         Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                                detect.bitsi.sendTrigger((session*10)+3) % sends a trigger marking the presence of a non CUE
                            end
                        else % for surround  oddball
                            if oddball_cfg(1) % this will only be 1 if it is an oddball detection block
                                %         Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                                %         Screen('FrameOval', win, [cue_color cue_color cue_color], out_VC_Rects, VC_thickness );
                                detect.bitsi.sendTrigger((session*10)+2) % sends a trigger marking the beginning of a CUE
                                
                            else
                                %        Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                                %        Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                                detect.bitsi.sendTrigger((session*10)+3) % sends a trigger marking the presence of a non CUE
                            end
                            
                        end
                        
                        Screen('DrawingFinished', win);
                        Screen('Flip', win);
                        
                        first_cue = 1;
                        Vprecue = 1; % avoids drawing the non cued stimulus
                    end
                    
                    if (time >=  (TRIAL_CUE1 + detect.CUEDUR )) % remove the cue
                        Vprecue = 2 ; % open again the function to draw for a second time the normal pedestal
                    end
                    
                    if (Vprecue == 0 | Vprecue == 2)
                        %  disp('pedestal ONE')
                        %   disp(time)
                        % Draw template
                        Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                        Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
                        %     Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                        %     Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                        Screen('DrawingFinished', win);
                        Screen('Flip', win);
                        
                        if (Vprecue == 0)
                            if (target_space == 1)
                                detect.bitsi.sendTrigger(100+(session*10)+1) % sends a trigger 111 marking a center template
                            else
                                detect.bitsi.sendTrigger(100+(session*10)+2)  % sends a trigger 112 marking a periphery template
                            end
                        end
                        
                        
                        if Vprecue == 2;
                            state = 'Play';
                        end % if this is the second time that the function draws the template, then it can go to play section
                        Vprecue = 1;
                    end
                    
                end
                
                
            case 'Play' % displays the target, the pedestal, and the second cue if necesary
                
                if (time >=  TRIAL_TARGET_2PED & Vflag == 0) % DISPLAYING VISUAL TARGET
                    %disp('pedestal & target')
                    %disp(time)
                    Screen('DrawTextures', win, gabortex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                    Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
                    %   Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                    %   Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                    Screen('DrawingFinished', win);
                    Vflag = 1; % Avoiding enter here and redraw the stimulus
                    
                    if time_calibration
                        Screen('FillOval', win, [0 0 0],centeredRect_calibration_cue, 400);
                    end
                    Screen('Flip', win);
                    
                    if (strcmp(stimulus_condition,'signal'))
                        detect.bitsi.sendTrigger(100+(session*10)+3) % sends a trigger marking the type of signal
                    else
                        detect.bitsi.sendTrigger(100+(session*10)+4) % sends a trigger marking the type of noise
                    end
                    
                end
                
                if (time >= TRIAL_TARGET_2PED + detect.STIMDUR) %  JUMP TO NEXT STAGE TO ERASE VISUAL TARGET
                    state = 'Pedestal_Two';
                end
                
            case 'Pedestal_Two'
                
                if (time >=  TRIAL_TARGET_2PED)
                    
                    if (time >=  TRIAL_CUE2  & second_cue == 0 )
                        Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                        Screen('FillOval', win,[0 0 0], centeredRect_fp, radious_fp);
                        
                        if target_space == 1 % for center oddball
                            if oddball_cfg(2) % this will only be 1 if it is an oddball detection block
                                %           Screen('FrameOval', win, [cue_color cue_color cue_color], in_VC_Rects, VC_thickness );
                                %           Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                                detect.bitsi.sendTrigger((session*10)+4) % sends a trigger marking the presence of a CUE
                                
                            else
                                %           Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                                %           Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                                detect.bitsi.sendTrigger((session*10)+5) % sends a trigger marking the presence of a CUE
                            end
                        else % for surround  oddball
                            if oddball_cfg(2) % this will only be 1 if it is an oddball detection block
                                %          Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                                %          Screen('FrameOval', win, [cue_color cue_color cue_color], out_VC_Rects, VC_thickness );
                                detect.bitsi.sendTrigger((session*10)+4) % sends a trigger marking the presence of a CUE
                            else
                                %         Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                                %         Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                                detect.bitsi.sendTrigger((session*10)+5) % sends a trigger marking the presence of a CUE
                                
                            end
                            
                        end
                        
                        Screen('DrawingFinished', win);
                        Screen('Flip', win);
                        
                        Vpostcue = 1; % avoids drawing the non cued stimulus
                        second_cue = 1;
                    end
                    
                    if (time >=  (TRIAL_CUE2 + detect.CUEDUR )) % remove the cue
                        Vpostcue = 2 ; % open again the function to draw for a second time the normal pedestal
                    end
                    
                    if (Vpostcue == 0 | Vpostcue == 2)
                        % Draw template
                        Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                        Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
                        %      Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                        %      Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                        Screen('DrawingFinished', win);
                        %disp('pedestal TWO')
                        Screen('Flip', win);
                        
                        if Vpostcue == 2;
                            state = 'Response_phase'; end % if this is the second time that the function draws the template, then it can go to response phas section
                        Vpostcue = 1;
                        % resuffle the spatial positions for the next phase of
                        % this trial
                    end
                    
                end
                
            case 'Response_phase'
                
                if (time >= TRIAL_RESPONSE)
                    
                    if  empty_resp == 0
                        detect.bitsi.clearResponses;  % empty bitsi buffer
                        empty_resp = 1;
                    end
                    
                    Screen('FillRect', win, gray);
                    Screen('TextSize',  win, 3*ppd);
                    
                    DrawFormattedText(win, response_options(1), wMid - (5*ppd), 'center',  [0,0,0]); % right response option
                    
                    DrawFormattedText(win, response_options(2), wMid + (3*ppd), 'center',  [0,0,0]); % left response option
                    Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
                    
                    Screen('Flip', win);
                    
                    if resp_trigger == 0 % only sends one trigger
                        detect.bitsi.sendTrigger((session*10)+6) % sends a trigger marking response phase
                        resp_trigger = 1;
                    end
                    
                    
                    % If the user press the key, record the response
                    
                    if ~keyIsDown % only enter if the characters has not been pressed. This allows to capture only one keypress
                        
                        [keyIsDown,secs,keycode] = KbCheck;
                        
                        % the option below is for mac when you have multiple
                        % keyboards (you have to specify what one)
                        % [keyIsDown,secs, keycode]= PsychHID('KbCheck', deviceNumber);
                        %RestrictKeysForKbCheck([KbName('ESCAPE'), KbName('SPACE'), KbName('UpArrow'),KbName('DownArrow')]); % select a group of keys to press. This
                        RestrictKeysForKbCheck([KbName('ESCAPE'), KbName('SPACE'), KbName('y'),KbName('u')]); % select a group of keys to press. This
                        
                        % WaitSecs(0.001); % delay to prevent CPU hogging
                        
                        [ triggerpressed] = getResponse(detect.bitsi,0.01,true); %
                        
                        if keyIsDown | triggerpressed > 0
                            % If the user holds down a key, KbCheck will report multiple events.
                            % To condense multiple 'keyDown' events into a single event, we wait until all
                            % keys have been released.
                            % KbReleaseWait;
                            response = 1; % responses confirmed
                            
                            
                            if MEG
                                secs =  absolute_time;
                                keyIsDown = 1;
                                %      detect.bitsi.clearResponses;
                            end
                            
                            
                            RT = ((secs - init_time)) - TRIAL_RESPONSE; % RT after the response options are shown
                            
                            % Note that we use find(keyCode) because keyCode is an array.
                            % See 'help KbCheck'
                            keypressed = KbName(keycode);
                            
                            % DONT FORGET TO REGISTER RTs
                            
                            %% Recoding button/key presses to response option
                            if (strcmp( keypressed , 'y') |  triggerpressed == 98 |  triggerpressed == 66 ) % sometimes the system takes 98 as 66
                                response_given =  response_options(1);% In this code, Y is left
                            end
                            if (strcmp( keypressed , 'u') |  triggerpressed == 99 |  triggerpressed == 67)
                                response_given =  response_options(2);% In this code, U is right (you can use button box commands)
                            end
                            
                            
                            
                            % disp( response_given )
                            
                            %% Cheking if the response is correct for the oddball task or visual detection task
                            if strcmp(procedure, 'oddball')
                                n_oddballs = sum(oddball_cfg);
                                if (n_oddballs == str2num(response_given))
                                    correct = 1;
                                    %disp(tilt);
                                    disp('correct!')
                                    detect.bitsi.sendTrigger(2) % sends a trigger marking a correct answer
                                    
                                    
                                    %disp([(secs - init_time),RT])
                                    
                                else
                                    correct = 0;
                                    disp('wrong!')
                                    detect.bitsi.sendTrigger(1) % sends a trigger marking a wrong answer
                                    RT = ((secs - init_time)) - TRIAL_RESPONSE; % RT after the response options are shown
                                    %disp([(secs - init_time),RT])
                                    
                                end
                                
                            else
                                if ((strcmp(stimulus_condition,'signal')  & strcmp( response_given , 'Y')) | (strcmp(stimulus_condition,'noise')  & strcmp( response_given , 'N')))
                                    correct = 1;
                                    %disp(tilt);
                                    detect.bitsi.sendTrigger(2) % sends a trigger marking a correct answer
                                    disp('correct!')
                                else
                                    correct = 0;
                                    detect.bitsi.sendTrigger(1) % sends a trigger marking a wrong answer
                                    disp('wrong!')
                                    
                                end
                                
                            end
                            
                            
                            if strcmp(procedure, 'detection')
                                if (strcmp(stimulus_condition,'signal')) % only update quest for signal trials
                                    
                                    if( strcmp(detect.visual_presentation(target_space),'center'))
                                        
                                        q_center_intersession = QuestUpdate(q_center_intersession,contrast_center,correct);
                                        % display the new contrast estimate
                                        disp(['intersessions new threshold = ', num2str(QuestQuantile(q_center_intersession))]);
                                        detect.q_center_intersession = q_center_intersession; % saving the quest
                                        detect.q_center_intersession_contrast = QuestQuantile(q_center_intersession); % saving final stimated contrast value
                                    else
                                        q_surround_intersession = QuestUpdate(q_surround_intersession,contrast_surround,correct);
                                        % display the new contrast estimate
                                        disp(['intersessions new threshold = ', num2str(QuestQuantile(q_surround_intersession))]);
                                        detect.q_surround_intersession = q_surround_intersession; % saving the quest
                                        detect.q_surround_intersession_contrast = QuestQuantile(q_surround_intersession);
                                    end
                                    
                                end
                                
                            end
                            
                            
                        end
                        
                        %disp( keypressed)
                        
                        if strcmp(KbName(keycode), 'ESCAPE');
                            Screen('Closeall')
                            % closing screen and audio
                            PsychPortAudio('Close', pahandle)
                            PsychPortAudio('Close', beephandle)
                            %break;
                        end
                        
                    end
                    
                    
                    if (time >=  TRIAL_FEEDBACK)
                        % disp('response')
                        %     disp(time)
                        state = 'Feedback';   % jump to the provide feedback phase and records trial information
                        
                        if  strcmp(response_given, 'NoResponse')
                            detect.bitsi.sendTrigger(3) % sends a trigger marking a that not response was registered
                        end
                        % s 1 session, 2 trial, 3stim time, 4 orientation,5postcue,6intensity, 7 resp, 8 RT, 9 correct
                        detect.responses_vector{iTrial,1} = session;
                        detect.responses_vector{iTrial,2} = iTrial;
                        detect.responses_vector{iTrial,3} = TRIAL_1PED; %% record postcue an all timing parameters for the MRI analysis
                        detect.responses_vector{iTrial,4} = TRIAL_CUE1; %% record postcue an all timing parameters for the MRI analysis
                        detect.responses_vector{iTrial,5} = TRIAL_TARGET_2PED; %% record postcue an all timing parameters for the MRI analysis
                        detect.responses_vector{iTrial,6} = TRIAL_CUE2; %% record postcue an all timing parameters for the MRI analysis
                        detect.responses_vector{iTrial,7} = TRIAL_RESPONSE; %% record postcue an all timing parameters for the MRI analysis
                        detect.responses_vector{iTrial,8} = char(stimulus_condition);
                        detect.responses_vector{iTrial,9} = char(modality_condition);
                        detect.responses_vector{iTrial,10} = char(detect.visual_presentation(target_space));
                        detect.responses_vector{iTrial,11} = rand_oddball; % if 1 is cue in 1 pedestal, 2 is cue in 2 pedestal, 3 is cue in both pedestals, 4 is nocue
                        detect.responses_vector{iTrial,12} = detect.orientations(target_definition);
                        detect.responses_vector{iTrial,13} =  contrast_trial; % saving contrast
                        detect.responses_vector{iTrial,14} = detect.cfg.freq; % saving contrast
                        detect.responses_vector{iTrial,15} = response_options;
                        detect.responses_vector{iTrial,16} = response_given;
                        detect.responses_vector{iTrial,17} = triggerpressed;
                        detect.responses_vector{iTrial,18} = correct;
                        detect.responses_vector{iTrial,19} = RT;
                        detect.responses_vector{iTrial,20} = char(initials);
                        
                        
                        %    Experimentalsession(session).quest = quest;
                        %     save(dataFile,'Experimentalsession');
                        
                        % calculating the proportion of correct responses at this
                        % trial
                        correct_responses = mean([detect.responses_vector{1:iTrial,18}]);
                        %find indexes for center and surround events
                        center_stimuli = find(detect.stim_combinations(1:iTrial,3)==1);
                        correct_center_responses = mean([detect.responses_vector{center_stimuli,18}]);
                        
                        surround_stimuli = find(detect.stim_combinations(1:iTrial,3)==2);
                        correct_surround_responses = mean([detect.responses_vector{surround_stimuli,18}]);
                        
                        
                        disp(['proportion of TOTAL correct responses = ',num2str(correct_responses)]);
                        disp(['proportion of correct responses center = ',num2str(correct_center_responses)]);
                        disp(['proportion of correct responses surround = ',num2str(correct_surround_responses)]);
                        
                        % printint test data Do not save in the experiment. Slows
                        
                        
                        fileID = fopen([dataadaptoutput],'wt'); % try to write on top
                        formatSpec = '%i %i %f %f %f %f %f %s %s %s %i %f %f %f %s %s %i %i %f %s \r\n';
                        formatSpec = '%i %i %f %f %f %f %f %s %s %s %i %f %f %f %s %s %i %i %f %s \r\n';
                        for row = 1 : nrows
                            fprintf(fileID, formatSpec, detect.responses_vector{row,:}) ;
                        end
                        fclose(fileID);
                    end
                end
                
            case 'Feedback' % after response time passes, provide feedback about response
                
                if response
                    if correct
                        color = [0,200 , 0]; % correct
                    else
                        color = [255, 0, 0]; % incorrect
                    end
                else
                    color = [255, 255 , 255]; % not response recorded
                end
                
                %        Screen('TextSize',  win, 2*ppd);
                %      keypressedis = ['triggerpress ', num2str(triggerpressed) ];
                %      DrawFormattedText(win, keypressedis , 'center', hMid-100, [200,0,0], max_n_ofcharperline,0,0);
                
                Screen('FillOval', win,  color, centeredRect_fp, radious_fp);
                Screen('Flip', win);
                % Screen('Close',  gabortex);
                
                if (time >=  TRIAL_ENDS)
                    %    disp('jump next TRIAL')
                    state = 'Initialize'; % lets prepare the next trial (going backwards)
                end
        end
        
        % Sound is played independently of the visual stimulus
        
        if (time >= (TRIAL_TARGET_2PED + AVdelay)  & Aflag == 0) % DISPLAYING SOUND
            if strcmp(modality_condition,'AV')
                % PsychPortAudio('Start', pahandle); % ensure aligment with oscilloscope - AVdelay variable
                PsychPortAudio('Start', beephandle); % ensure aligment with oscilloscope - AVdelay variable
                detect.bitsi.sendTrigger((session*10)+7) % sends a trigger marking the sound presentation
                Aflag = 1;
            else
                detect.bitsi.sendTrigger((session*10)+8) % sends a trigger marking no sound presentation
                Aflag = 1;
            end
            
        end
        
    end
    
    
    Experimentalsession(session).detect = detect;
    save(dataFile,'Experimentalsession');
    
    % save data in a txt file
    % printint test data
    
    fileID = fopen([dataadaptoutput],'wt'); % try to write on top
    
    formatSpec = '%i %i %f %f %f %f %f %s %s %s %i %f %f %f %s %s %i %i %f %s \r\n';
    % for windows
    %dataadaptoutput = [path,'\results\',initials,'_detect.txt'];
    
    for row = 1 : nrows
        fprintf(fileID, formatSpec, detect.responses_vector{row,:}) ;
    end
    
    fclose(fileID);
    
    
    %Experimentalsession(session).detect = detect;
    %save(dataFile,'Experimentalsession');
    
    toc
    Blockfinish;
    KbWait();
    
    disp(['intersessions new threshold = ', num2str(QuestQuantile(q_center_intersession))]);
    disp(['intersessions new threshold = ', num2str(QuestQuantile(q_surround_intersession))]);
    
    detect.bitsi.close();
    Screen('Closeall')
    % closing screen and audio
    PsychPortAudio('Close', pahandle)
    PsychPortAudio('Close',beephandle)
    
end
%fprintf('%5.2f	%4.1f	%5.2f\n',tActual,q.beta,q.gamma);


