function detect_AVE(task)

% to calibrate time use this and it will draw very big patches on the
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
%% Switch to Screen('Preference', 'SkipSyncTests', 1) if trouble starting program
%initimeexp = GetSecs; % get time to measure the total lenght of the experiment
%% Set current directory
%path = '/Volumes/ALEX_EXT/Experimentos/DondersExperiments' %
path = mfilename('fullpath'); % Finds current location of file
[pathstr] = fileparts(path); % Removes file information, saves directory
path = pathstr;

%% Collect subject info

% universal file separator across operative systems is filesep

if strcmp(task, 'demo')
    prompt = {'Contrast center (0.2):','Contrast surround (0.2):'};
    input = inputdlg(prompt);
    initials = task; % Collect subject initials
    contrast_center = input{1};
    contrast_center  = str2num(contrast_center);
    contrast_surround = input{2};
    contrast_surround  = str2num(contrast_surround);
    
else
    prompt = {'Enter subject initials: ','Contrast center (0.2):','Contrast surround (0.2):'};
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
else
    session = 1; % fill with 1 the session number just to not leave the field empty
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
Experimentalsession(session).date = date;


% creating the  subject file;
save(dataFile,'Experimentalsession');


init_AVE_parameters_v2;

detect = init;
% If MAC, select the external device
%[keyboardIndices, productNames] = GetKeyboardIndices
%deviceNumber = 2;

% PREPARING STIMULI COMBINATIONS

% reps =  12;
reps =  18;

if strcmp(task, 'demo')
reps =  4;

end

%detect.orientations = [90]; % degrees of possible gabor pathces
%detect.pre_cue = [90 90 315 45]; %orientation to be reported
%detect.post_cue = [315 45 315 45]; %orientation to be reported

% in the informative trials, the orientation is determined by the
detect.stim_conditions = {'signal','noise'};
%detect.combination_conditions = {'A_irrel_CCW','A_irrel_CW','A_rel_CCW','A_rel_CW'};
detect.modality_conditions = {'V','AV'}; % 'VC' % for this experiment I removed the VC condition
detect.visual_presentation = {'center','surround'};
% init parameters for the 1st pre adaptation test phase
[stim_rand_order,stim_combinations] = combinator([length(detect.orientations),length(detect.modality_conditions),length(detect.visual_presentation)],reps);
%randomize order
stim_combinations = stim_combinations(stim_rand_order,:);

detect.stim_combinations = stim_combinations;
detect.repetitions = reps;
detect.nTrials = length(stim_combinations);


%% Experiment timings
for t_idx = 1:detect.nTrials % consider randomizing a littlebit the ITI in each trial
    detect.TRIAL_TIMES(t_idx) = detect.INIT_DELAY +...
        (t_idx-1)*(detect.FIXATION + detect.PRE_STIM_CUE  +  detect.POST_STIM_CUE + detect.RESPONSETIME + detect.ITI); %
end

curr_sched =   detect.TRIAL_TIMES ; %Shift all the times by the INIT DELAY

%% RECORDING DATA IN THIS TEMPLATE MATRIX
detect.responses_vector = zeros(detect.nTrials,15); % vector to record conditions 1 session, 2 trial, 3stim time, 4 orientation,5postcue,6intensity, 7 resp, 8 RT, 9 correct
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
rect=[0 0 800 600]; % mini-screen.
%rect=[0 0 320 240]; % mini-screen.

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
audiodevice = 1;
pahandle = PsychPortAudio('Open',audiodevice, [], 1, fs, nrchannels);
beephandle = PsychPortAudio('Open',audiodevice, [], 1, fs, nrchannels);

% load memory buffers with auditory stimulus. Ready to fire!.
PsychPortAudio('FillBuffer', pahandle, white_noise); % Auditory stimulus
PsychPortAudio('FillBuffer',beephandle, beep_noise); % warning signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initializing some screen basic parameters

[w, h] = RectSize(winRect);
wMid = w / 2;
hMid = h / 2;

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
mask(1).filter = mask(1).filter - (sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)<=mask_radius/20); % add inner empty area
mask(2).filter = sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)>=mask_radius/2;  %surround in pixels
                
% the area of stimulation in the center is 4 times smaller than in the
% surround (fovea see 2 mm)
%sum(sum(mask(2).filter) )
%sum(sum(mask(1).filter) )
% precompute signal energy filters
%gaborang = pi/180*(70:1:110); % filter orientations (rad) % calculate only a small part of the spectrum. Then you can calculate the rest from the images.
%gaborang = pi/180*(60:1:120); % filter orientations (rad) % calculate only a small part of the spectrum. Then you can calculate the rest from the images.
gaborang = pi/180*(1:1:180); % filter orientations (rad) % calculate only a small part of the spectrum. Then you can calculate the rest from the images.

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
          %save(surround_filters_file,'cfg_surround_e');
          
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
    
        % presenting center or surround
    if  strcmp(detect.visual_presentation(target_space),'center')    
        spatial_filter = mask(1).filter;
        indexes_sf = find(spatial_filter == 0);
    else
        spatial_filter = mask(2).filter;
        indexes_sf = find(spatial_filter == 0);
    end
    
                                       
% Now I create the SANDWICH noise mask fat 0
% contrast for SF (for the target)
                
                  while true
                    cfg.gaborcon = 0;
                    noiseimg = make_gabor_and_noise(cfg);
                    noiseimg(indexes_sf ) = 0.5;
                    %[patchimg] = absoluteDifferenceBetweenWhiteAndGray * patchimg;
                    noiseimg = (gray * noiseimg)/0.5;
                    
                    if all(noiseimg(:) > black & noiseimg(:) < white)
                        break
                    end
                end
                            
     % assigning contrast values for the two visual field targets     
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
    
    cfg.gaborphi = rand; % Gabor unit phase (randomize in each trial)
    % clip patch luminance to [0,1] range to avoid ugly peaks of contrast in the
    % patch
    
    percetange_stimuli = (i/nimages)*100;
    Creating_stimuli_screen;

% Now I create the target stimilus              
    while true
        
     [img] = make_gabor_and_noise(cfg);
     %[patchimg] = absoluteDifferenceBetweenWhiteAndGray *
     %patchimg;\
    % making not stimulus part like grey background
       img(indexes_sf ) = 0.5;
   
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
    stimuli(i).img = img; % dont save in the same file, it is to heavy
    noises(i).noiseimg = noiseimg; % dont save in the same file, it is to heavy
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
KbWait(); % press a key to continue with the script

if strcmp(task, 'demo') 
examplestimuli;
KbWait(); % press a key to continue with the script
Screen('Close', signalexampletex)
Screen('Close', noisexampletex)
end

%% Initialize quest parameters THIS QUEST IS NOT VERY USEFUL BECAUSE IM JUST RUNNING ONE IN BOTH, CENTER AND SURROUND TOGETHER
%tGuess = 0.1;
%tGuessSd = 0.075

tGuess = cfg.gaborcon;
tGuessSd = 0.025;

pThreshold=0.72; % (previous was 0.75 and gamma 0.5)
beta=3.5;delta=0.01;
gamma=0.5; % gamma is the proportion correct of responses that are expected when stimulus is not present. slope, threshold and minimum value (it is a 3AFC, noise, 2 x tilt )
% the quest is designed only to update when there was a signal CW or CCW.
% Si el sujeto no ve nada, juzgara aleatoriamente, asi que puede acertara
% el 50% de las veces
see_quest = 0;
range = 0.25; % range below and above threshold to be tested
grain = 0.001;
q=QuestCreate(tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range,see_quest);
q.normalizePdf=1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.

tTest=tGuess; % The guess is equal to the contrast in this run -> assign here a fix contrast value and do not update the contrast (but update the quest)

% preparing the quest that runs across different runs
if (session == 1)
   % If this is the first run, initialize the intersession quest 
   q_intersession = q; 
else
    %If this is not the first run, continue using the quest from the
    %previous run
   q_intersession = q;
  % q_intersession = Experimentalsession(1).detect.q_intersession;
end

% thus q is the quest that it is applied in this run, and q_intersession is
% the quest across sessions.


% POSITION COORDINATES FOR CIRCULAR OBJECTS(FP AND INNER ANNULUS)

color_fp = [0 0 0];
radious_fp = 20;
baseRect_fp = [0 0 radious_fp  radious_fp ]; % create a fake box to set the center of ovals
centeredRect_fp = CenterRectOnPointd(baseRect_fp, w/2, h/2);

radious_center = 50;
baseRect_center = [0 0 radious_center  radious_center ]; % create a fake box to set the center of ovals
centeredRect_center = CenterRectOnPointd(baseRect_center, w/2, h/2);

% Calibration cue. Use to calibrate the timing wuith the oscilloscope
baseRect_calibration_cue = [0 0 500 500 ]; % create a fake box to set the center of ovals
centeredRect_calibration_cue = CenterRectOnPointd(baseRect_calibration_cue, w/2, h/2);

% Initialize loop parameters
iTrial = 0;
state = 'Initialize';
exit_loop = 0;

Screen('FillRect', win, gray);    % Gray background
vbl = Screen('Flip', win);
ifi = Screen('GetFlipInterval', win);

% use this line to draw and avoid tearing effects ->  vbl = Screen('Flip', win, vbl + 0.5 * ifi);

%tstart = vbl;
%count_frames = 0;

init_time = GetSecs;

while ~exit_loop % exit the loop when trials are finish
    
    % Updating time
    
    relative_time = GetSecs;
    
    time = relative_time - init_time; % measure absolute elapsed time
    
    switch state
        
        case 'Initialize'
            
            disp(iTrial)
            
            if (iTrial  >= detect.nTrials)
                exit_loop = 1;
                
            else % buffer stim and get time for next stimulus presentation
                
                iTrial = iTrial + 1;
                % Trial conditions
                %  target_definition = detect.stim_combinations(iTrial,1); % select orientation target
                %  detect.stim_conditions(target_definition);
                
                % select pre-postcue combinations and modality/
                stimulus_definition = detect.stim_combinations(iTrial,1); % select orientation target
                stimulus_condition = detect.stim_conditions(stimulus_definition);
                modality_definition = detect.stim_combinations(iTrial,2); %
                modality_condition = detect.modality_conditions(modality_definition);
                target_space = detect.stim_combinations(iTrial,3); 
                
                if target_space == 1
                  contrast_trial = contrast_center;
                else
                    contrast_trial = contrast_surround;
                    
                end
                
                % modality_condition = 'AV';
                %  pre_cue_orientation = detect.pre_cue(combination_definition); % select orientation target
                %  post_cue_orientation = detect.post_cue(combination_definition); % select orientation target                                
                % Select trial timings
                Fixation_point_time = curr_sched(iTrial);
                
                Pre_stim_cue_time = curr_sched(iTrial) + detect.FIXATION;
                Pre_stim_cue_time_off = Pre_stim_cue_time + detect.CUEDUR; % switch off the cue
                
                % jitter 500 ms [- and + 250 ms]
                visual_stimulus_jitter = detect.jitter_range + (-detect.jitter_range-detect.jitter_range).*rand(1,1); % visual stimulus can be jittered between -.2 and 0.2 s
                % check that this visual_stimulus_jitter value is not
                % smaller or larget than the pedestal presentation time.
                V_stimulus_time = Pre_stim_cue_time + detect.PRE_STIM_CUE ; % this is the time where the visual stimulus should be presented in this trial
                real_stimulus_time = V_stimulus_time + visual_stimulus_jitter;
                real_stimulus_time_off =  real_stimulus_time + detect.STIMDUR; % this is the time where the visual stimulus should be removel
                
                Post_stim_cue_time = V_stimulus_time + detect.POST_STIM_CUE; % after postcue, response should be collected
                Post_stim_cue_time_off = Post_stim_cue_time + detect.CUEDUR; % switch off the cue
                
                
                Feedback_time = Post_stim_cue_time + detect.RESPONSETIME; % provide feedback about correctness
                
                % here draw mask texture
                tilt = pi/180*detect.orientations(stimulus_definition);
                cfg.gaborang = tilt; % Gabor orientation (rad)
                cfg.gaborphi = rand; % pi/180*90; % Gabor unit phase CHECK THIS (radians?)
                cfg.gaborcon = 0; % Gabor Michelson contrast
                
                % clip patch luminance to [black,white] range to avoid black clusters at
                % high contrasts
                
                while true
                    img = make_gabor_noise(cfg);
                    %[patchimg] = absoluteDifferenceBetweenWhiteAndGray * patchimg;
                    img = (gray * img)/0.5;
                    if all(img(:) > black & img(:) < white)
                        break
                    end
                end
               
                masktext = Screen('MakeTexture', win,  img); %Texture generated
      
                % here prepare target texture 
                % all the target textures were pre-generatedthe at the begining). Save the matrix always
                % for regression analyses
                
                % Creating texture
                gabortex=Screen('MakeTexture', win,   stimuli(iTrial).img); %Texture generated
                noisetex = Screen('MakeTexture', win,  noises(iTrial).noiseimg); %Texture noise generated
                % imagesc(n)
                texrect = Screen('Rect', gabortex); %Extract information about texture size
                scale = 1;
                posx = [w/2];   %coordinates
                posy = [h/2];
                poscenter = [w/2, h/2];
                rotAngles = [0]; % we create the gratings already rotated, so no rotation parameter should be applied to the texture.                
                dstRects = CenterRectOnPoint(texrect*scale , posx, posy);
                
                %Initializing flags / We only want the functions to draw
                %each item one time to avoid delays
                Vflag = 0;
                Aflag = 0;
                Vprecue = 0;
                Vpostcue = 0;
                response = 0;
                keycode = [];
                correct = 0;
                RT = 0;
                keypressed = 'NoResponse';
                state = 'Fixationpoint';
                keyIsDown = 0;
                
            end
            
        case 'Fixationpoint'
            if (time >=  Fixation_point_time)
                % Display fixaxtion point
                
                %Screen('TextSize',  win, 30);
                %DrawFormattedText(win,'+', 'center','center', [240,240,240]);
                Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp); % new and better centered fixation point
                                
                if time_calibration
                    Screen('FillOval', win, [255 255 255],centeredRect_calibration_cue, 400);
                end
                
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                state = 'Precue';
            end
            
            
        case 'Precue'
            
            if (time >=  Pre_stim_cue_time & Vprecue == 0) % DRAW VISUAL CUE
                
                % Non attentional relevant
                %if pre_cue_o rientation == 90
                %   Screen('DrawLine', win, [255,255,255], (w/2)-16, (h/2),(w/2)+16,(h/2), 6);
                %Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                % end
                % Attentional relevant
                %  if pre_cue_orientation == 45
                %      Screen('DrawLine', win, [255,255,255], (w/2)-11, (h/2)+11,(w/2)+11,(h/2)-11, 6);
                %Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                %  end
                %  if pre_cue_orientation == 315
                %     Screen('DrawLine', win, [255,255,255], (w/2)-11, (h/2)-11,(w/2)+11,(h/2)+11, 6);
                %
                %  end
                 Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround                
                %DrawFormattedText(win,'-', 'center','center', [240,240,240]); % in quest only one type of cue
                % Screen('DrawLine', win, [255,255,255], (w/2)-16, (h/2),(w/2)+16,(h/2), 6);
                Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
                Screen('DrawingFinished', win);
                Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                Vprecue = 1;
                 state = 'Play';
            end
            
%             if (time >=  Pre_stim_cue_time_off) % ERASE VISUAL CUE,RESTORE NORMAL FIXATION AND MOVE TO STIM PLAYING.
%                 Screen('TextSize',  win, 30);
%                 Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
%                 vbl = Screen('Flip', win, vbl + 0.5 * ifi);
%                 state = 'Play';
%                 
%             end
            
            
        case 'Play'
                
            if (time >= real_stimulus_time & Vflag == 0) % DISPLAYING VISUAL TARGET
             %   gabortex=Screen('MakeTexture', win,   stimuli(12).img); %Texture generated
               
                Screen('DrawTextures', win, gabortex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                Screen('TextSize',  win, 30);
                % Screen('FillOval', win, gray, centeredRect_center, radious_center);
              
                
                if strcmp(modality_condition,'VC')
                  Screen('FillOval', win, [200 200 200] , centeredRect_fp, radious_fp); % color_fp    
                else
                  Screen('FillOval', win, color_fp , centeredRect_fp, radious_fp); %   
                 
                end
                 
                Screen('DrawingFinished', win);
                Vflag = 1; % Avoiding enter here and redraw the stimulus
              
                if time_calibration
                    Screen('FillOval', win, [0 0 0],centeredRect_calibration_cue, 400);
                end
                
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
            end
            
            if (time >= real_stimulus_time_off) % ERASE VISUAL TARGET and REDRAW ONLY NOISE
             %  Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround                
                %DrawFormattedText(win,'-', 'center','center', [240,240,240]); % in quest only one type of cue
                % Screen('DrawLine', win, [255,255,255], (w/2)-16, (h/2),(w/2)+16,(h/2), 6);
                Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
              %  Screen('DrawingFinished', win);
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                state = 'Postcue';
                
            end
            
        case 'Postcue'
            
            if (time >=  Post_stim_cue_time & Vpostcue == 0) % If both stimuli displayed, move to the next phase
                
                %DrawFormattedText(win,cue_sign{post_cue}, 'center','center', [240,240,240]);
                % if post_cue_orientation == 45
                %     Screen('DrawLine', win, [255,255,255], (w/2)-11, (h/2)+11,(w/2)+11,(h/2)-11, 6);
                %     %Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                %  end
                %  if post_cue_orientation == 315
                %      Screen('DrawLine', win, [255,255,255], (w/2)-11, (h/2)-11,(w/2)+11,(h/2)+11, 6);
                %
                %   end
                Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                                
                if time_calibration
                    Screen('FillOval', win, [255 255 255],centeredRect_calibration_cue, 400);
                end
                
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                Vpostcue = 1;
            end
            
            if (time <=  Feedback_time)
                
                % If the user press the key, record the response
                
                if ~keyIsDown % only enter if the characters has not been pressed. This allows to capture only one keypress
                    
                    [keyIsDown,secs,keycode] = KbCheck;
                    
                    % the option below is for mac when you have multiple
                    % keyboards (you have to specify what one)
                    % [keyIsDown,secs, keycode]= PsychHID('KbCheck', deviceNumber);
                    %RestrictKeysForKbCheck([KbName('ESCAPE'), KbName('SPACE'), KbName('UpArrow'),KbName('DownArrow')]); % select a group of keys to press. This
                    RestrictKeysForKbCheck([KbName('ESCAPE'), KbName('SPACE'), KbName('y'),KbName('u')]); % select a group of keys to press. This
                    
                    WaitSecs(0.001); % delay to prevent CPU hogging
                    
                    if keyIsDown
                        % If the user holds down a key, KbCheck will report multiple events.
                        % To condense multiple 'keyDown' events into a single event, we wait until all
                        % keys have been released.
                        % KbReleaseWait;
                        response = 1; % responses confirmed
                        
                        RT = ( (secs - init_time)) - V_stimulus_time;
                        %disp([(secs - init_time),RT])
                        
                        % Note that we use find(keyCode) because keyCode is an array.
                        % See 'help KbCheck'
                        keypressed = KbName(keycode);
                        % DONT FORGET TO REGISTER RTs
                        % correct or incorrect?
                        
                        if (strcmp(stimulus_condition,'signal')  & strcmp( keypressed , 'y') | strcmp(stimulus_condition,'noise')  & strcmp(  keypressed , 'u'))
                            correct = 1;
                            %disp(tilt);
                            disp('correct!')
                        else
                            correct = 0;
                            disp('wrong!')
                        end
                        
                        q_intersession = QuestUpdate(q_intersession,tTest,correct);
                       % display the new contrast estimate
                        disp(['intersessions new threshold = ', num2str(QuestQuantile(q_intersession))]);
                        detect.q_intersession = q_intersession; % saving the quest
                        detect.q_intersession_contrast = QuestQuantile(q_intersession);
                        
                        q = QuestUpdate(q,tTest,correct);
                        disp(['blockwise threshold = ', num2str(QuestQuantile(q))]);
                        detect.q = q; % saving the quest adn the estimated contrast
                        detect.q_contrast = QuestQuantile(q);
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
                
            end

            if (time >  Feedback_time)
                state = 'Feedback';   % jump to the provide feedback phase and records trial information
                % s 1 session, 2 trial, 3stim time, 4 orientation,5postcue,6intensity, 7 resp, 8 RT, 9 correct
                detect.responses_vector{iTrial,1} = session;
                detect.responses_vector{iTrial,2} = iTrial;
                detect.responses_vector{iTrial,3} = Pre_stim_cue_time; %% record postcue an all timing parameters for the MRI analysis
                detect.responses_vector{iTrial,4} = real_stimulus_time; %% record postcue an all timing parameters for the MRI analysis
                detect.responses_vector{iTrial,5} = Post_stim_cue_time; %% record postcue an all timing parameters for the MRI analysis
                detect.responses_vector{iTrial,6} = char(stimulus_condition);
                detect.responses_vector{iTrial,7} = char(modality_condition);
                detect.responses_vector{iTrial,8} = char(detect.visual_presentation(target_space));
                detect.responses_vector{iTrial,9} = detect.orientations(target_definition);
                detect.responses_vector{iTrial,10} =  contrast_trial; % saving contrast
                detect.responses_vector{iTrial,11} = detect.cfg.freq; % saving contrast
                detect.responses_vector{iTrial,12} = keypressed;
                detect.responses_vector{iTrial,13} = correct;
                detect.responses_vector{iTrial,14} = RT;
                detect.responses_vector{iTrial,15} = char(initials);
                %    Experimentalsession(session).quest = quest;
                %     save(dataFile,'Experimentalsession');
              
                % calculating the proportion of correct responses at this
                % trial
                correct_responses = mean([detect.responses_vector{:,13}]);
                %find indexes for center and surround events
                center_stimuli = find(detect.stim_combinations(:,3)==1);
                correct_center_responses = mean([detect.responses_vector{center_stimuli,13}]);
                
                surround_stimuli = find(detect.stim_combinations(:,3)==2);
                correct_surround_responses = mean([detect.responses_vector{surround_stimuli,13}]);
               
                % if the participants has perform a multiple of 15 trials
                % and performance is perfect, adjust level of contrast. A
                % beep signal marks the problem
%                 if( ~mod(iTrial,25) & correct_responses > 0.95)
%                      PsychPortAudio('Start',beephandle); % ensure aligment with oscilloscope - AVdelay variable
%               % Beeper(2000,3,0.5); AVOID BEEPER. INTERACTS WITH
%               % PSYCHPORTAUDIO BADLY
%                 end
%                               
%                 if( ~mod(iTrial,25) & correct_responses < 0.54)
%                     PsychPortAudio('Start',beephandle); 
%                 end
                
              %  PsychPortAudio('Start',beephandle); %
                       
                disp(['proportion of TOTAL correct responses = ',num2str(correct_responses)]);            
                disp(['proportion of correct responses center = ',num2str(correct_center_responses)]);
                disp(['proportion of correct responses surround = ',num2str(correct_surround_responses)]);
                

                dataadaptoutput = [path, filesep 'results',filesep initials,'_session_',num2str(session),'_AVE_detect.txt'];
                % printint test data Do not save in the experiment. Slows
                % down the program, possibly because it saves the data in
                % the intranet
          %      fileID = fopen([dataadaptoutput],'wt'); % try to write on top
          %      formatSpec = '%i %i %f %f %f %s %s %f %f %s %i %f %s \r\n';
                                
          %      for row = 1 : nrows
          %          fprintf(fileID, formatSpec, detect.responses_vector{row,:}) ;
          %      end
                
           %     fclose(fileID);
                
            end
            
        case 'Feedback' % after response time passes, provide feedback about response
            
            if (time >=  Feedback_time)
                
                if response
                    if correct
                        color = [0,200 , 0]; % correct
                    else
                        color = [255, 0, 0]; % incorrect
                    end
                else
                    color = [255, 255 , 255]; % not response recorded
                end
                
                 
                Screen('DrawTextures', win, masktext,[],dstRects,rotAngles); %, 2 is the adaptor surround
                Screen('TextSize',  win, 30);
                % Screen('FillOval', win, gray, centeredRect_center, radious_center);                
                FeedBack_size = [0, 0, 30, 30];
                square_rect = CenterRect(FeedBack_size, winRect);
                % Screen('FillRect',win, color, square_rect);
                Screen('FrameRect', win, color, square_rect, 2);
                Screen('FillOval', win, [0,0,255], centeredRect_fp, radious_fp);
                Screen('DrawingFinished', win);
                
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                
                Screen('Close',  masktext);
                Screen('Close',  gabortex);
                state = 'Initialize'; % lets prepare the next trial (going backwards)
                
            end

    end
       
    % Sound is played independently of the visual stimulus

    if (time >= (real_stimulus_time + AVdelay) & strcmp(modality_condition,'AV') & Aflag == 0) % DISPLAYING SOUND
       % PsychPortAudio('Start', pahandle); % ensure aligment with oscilloscope - AVdelay variable
       PsychPortAudio('Start', beephandle); % ensure aligment with oscilloscope - AVdelay variable
       
        Aflag = 1;
    end
    
    %count_frames = count_frames + 1;
end

Experimentalsession(session).detect = detect;
save(dataFile,'Experimentalsession');


% save data in a txt file
dataadaptoutput = [path, filesep 'results' filesep ,initials,'_session_',num2str(session),'_AVE_detect.txt'];
% printint test data
fileID = fopen([dataadaptoutput],'wt'); % try to write on top
formatSpec = '%i %i %f %f %f %s %s %s %f %f %f %s %i %f %s \r\n';  
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

Screen('Closeall')
% closing screen and audio
PsychPortAudio('Close', pahandle)
PsychPortAudio('Close',beephandle)

end

%fprintf('%5.2f	%4.1f	%5.2f\n',tActual,q.beta,q.gamma);
