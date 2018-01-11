function quest_AVE()

%% Alexis Pérez Bellido 20-5-2017
% this function run 2 possible task, although is specifically designed for QUEST and visual detection: 1 task is a
% visual grating detection task, where the participants have to detect visual
% grating. 2 task is a oddball counting task, where the participants have
% to attend to rings which are black.
% Stimulation can be visual.

%% taking pictures 
% imageArray = Screen('GetImage',win);
% imwrite(imageArray, [path, filesep 'results' filesep 'center1pedestal.jpg'])

% to calibrate time use this as 1 and it will draw very big patches on the
% screen that can be measure with an oscilloscope
time_calibration = 0;

% Main experiment. Use contrast from quest to adjust stimulus to threshold
% level. This function precreates gabor patches with particular
% characteristics and calculate their energy for each orientation.
% Provide feedback Use detect_AVE('demo') for demo or quest_AVE('Exp') for
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

    
    prompt = {'Enter subject initials: ','Contrast center (0.007):','Contrast surround (0.007):'};
    input = inputdlg(prompt);
    initials = input{1}; % Collect subject initials REMEMBER ATACCH ALL THE
    contrast_center = input{2};
    contrast_center  = str2num(contrast_center);
    contrast_surround = input{3};
    contrast_surround  = str2num(contrast_surround);


outputFolder = [path, filesep 'results' filesep ,initials,filesep];

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

dataFile = [outputFolder initials,'_AVE_quest.mat'];

if exist(dataFile, 'file')
    disp('This subject has been already tested, checking session');
    load(dataFile);
    n_sessions_runned = length(Experimentalsession);
    uiwait(msgbox(['Session tested previously = ', num2str(n_sessions_runned)],'modal'));
end


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


% folder to save the raw images
imgFile = [outputFolder 'imgfiles' filesep,'/'];

if ~exist(imgFile,'dir')
    mkdir(imgFile);
end

rawimg_fname = ['IMG_run_',num2str(session),'.mat']

Experimentalsession(session).date = datetime;


procedure = 'detection'; % or vquestion or oddball task (this is the quest, so there is no oddball task

% creating the  subject file;
save(dataFile,'Experimentalsession');

% loading general parameters
init_AVE_parameters_pedestal_meg;

quest = init;

% If MAC, select the external device
%[keyboardIndices, productNames] = GetKeyboardIndices
%deviceNumber = 2;

% PREPARING STIMULI COMBINATIONS

reps =  18;
%reps =  18;

% MEG or not MEG
MEG = 1;

% bitsi stuff and trigger coding
if MEG
quest.bitsi = Bitsi('com1');
quest.triggers = def_detect_AVE_MEG_Triggers();
else
 quest.bitsi = Bitsi('');  
  quest.triggers = def_detect_AVE_MEG_Triggers();
end


%quest.orientations = [90]; % degrees of possible gabor pathces
%quest.pre_cue = [90 90 315 45]; %orientation to be reported
%quest.post_cue = [315 45 315 45]; %orientation to be reported

% in the informative trials, the orientation is determined by the
quest.stim_conditions = {'signal','noise'};
%quest.combination_conditions = {'A_irrel_CCW','A_irrel_CW','A_rel_CCW','A_rel_CW'};
quest.modality_conditions = {'V'}; % ,'AV'
quest.visual_presentation = {'center','surround'};



% Combinations of conditions and randomization
allCombinations = combvec([1:length(quest.stim_conditions)],[1:length(quest.modality_conditions)],[1:length(quest.visual_presentation)]);
allCombinations = allCombinations'; % transposing all combinations matrix
sizeComb        = length(allCombinations);
nTrials         = reps*sizeComb;
allCombinations = repmat(allCombinations,reps,1);
%generic variable to randomize inside the allCombination

rand_trials_vector = 1 : nTrials;
rand_trials_vector = rand_trials_vector(randperm(nTrials));  % This will serve to select the corresponding arrow during the experiment.

%randomize order
stim_combinations = allCombinations(rand_trials_vector,:);

quest.stim_combinations = stim_combinations;
quest.repetitions = reps;
quest.nTrials = nTrials;


% creating vector with oddball trials for the attention experiment

if strcmp(procedure,'oddball')
    quest.proportion_oddball = 0.1;
else
    quest.proportion_oddball = 0.0; % no oddball in the visual question task
end

n_oddballs = round(quest.proportion_oddball*quest.nTrials); % 10 % of the trials are oddballs in the attention task
quest.oddballs_vector = [ones(1,n_oddballs),zeros(1,quest.nTrials-n_oddballs)];
quest.oddballs_vector = quest.oddballs_vector(randperm(quest.nTrials));
quest.oddtions = [1 0; 0 1; 1 1]; % possible oddball trials (it presents the possiblities that the 1st, 2nd or both cues change to white)


%% Experiment events timings TIMINGS FOR QUEST ARE A LITTLE BIT FASTER THAN FOR MAIN TASK
for t_idx = 1 : quest.nTrials % consider randomizing a little bit the ITI in each trial
    quest.INIT_DELAY = 2 % time before the first trial
  %  quest.FIXATION(t_idx) = 0.75 + (1.0-0.75)*rand(1,1); % intertrial time (fixation point)
    quest.FIXATION(t_idx) = 0.5 + (0.75-0.5)*rand(1,1); % intertrial time (fixation point)
    
    quest.PRE_STIM_TEMPLATE(t_idx) = 1.0 + (1.5-1.0)*rand(1,1); % when this time finishes the target is drawn
    quest.CUE1(t_idx) = 0.5 + (0.8-0.5)*rand(1,1); % time to flicker the first visual distractor
    
    quest.POST_STIM_TEMPLATE(t_idx) = 1.0; % after this, time to responde
    quest.CUE2(t_idx) = 0.5 + (0.75-0.5)*rand(1,1); % time to flicker the second visual distractor
    
    quest.RESPONSETIME(t_idx) = 1.5;
    quest.FEEDBACK(t_idx) = 0.25;
    
    if (t_idx== 1) % previously the function did not create the timings correctly, this is my modification
        quest.TRIAL_TIMES(t_idx) = quest.INIT_DELAY;
    else
        
        quest.TRIAL_TIMES(t_idx) = quest.TRIAL_TIMES(t_idx-1) + ...
            quest.FIXATION(t_idx) + quest.PRE_STIM_TEMPLATE(t_idx)  + ...
            quest.POST_STIM_TEMPLATE(t_idx) + quest.RESPONSETIME(t_idx) + ...
            quest.FEEDBACK(t_idx); %
        
    end
    
end

% quest.INIT_DELAY   + (t_idx-1)*(quest.FIXATION(t_idx) + quest.PRE_STIM_TEMPLATE(t_idx)  +  quest.POST_STIM_TEMPLATE(t_idx) + quest.RESPONSETIME(t_idx)+quest.FEEDBACK(t_idx)); %

%% RECORDING DATA IN THIS TEMPLATE MATRIX
quest.responses_vector = zeros(quest.nTrials,20); % vector to record conditions 1 session, 2 trial, 3stim time, 4 orientation,5postcue,6intensity, 7 resp, 8 RT, 9 correct
quest.responses_vector = num2cell(quest.responses_vector);
% these will be my variables
quest.responses_vector_variables = {'session','trial','Pre_time','Stim_time','grating_tilt','contrast','SF','visual_space','keypressed','correct','RT','subject'};
           
[nrows,ncols] = size(quest.responses_vector); % this makes very slow the process. Record the data at the end in txt file

%% setup keyboard basic parameters

KbName('UnifyKeyNames')

%will avoid subject pressing incorrect buttons. It can also improve
%computer performance

%% Video set up %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check script is running on Psychtoolbox-3

% PsychDebugWindowConfiguration % transparent window to debug

% Screen('Preference', 'SkipSyncTests', 1); % In lab I dont need to skip syn

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
% 
% %% Auditory stimuli %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASTIMDUR  = quest.ASTIMDUR; %sec.
% fs = quest.fs; %sampling fq
% % creating the whitenoise
% white_noise = (rand(1,round(ASTIMDUR*fs))*2-1); % LETS create some whitenoise.
% [white_noise] = linramp(white_noise,ASTIMDUR,fs,0.015); % apply ramp on/off
% 
% % create a warning beep (beeper works very bad so I create my own
% % stimulus
% freq = 1000; % hz
% beep_noise = sin(linspace(0, ASTIMDUR *freq*2*pi, round(ASTIMDUR*fs))); 
% [beep_noise] = linramp(beep_noise,ASTIMDUR,fs,0.005); % apply ramp on/of
% AVdelay = -0.02; % negative values sound before the visual, postive values after
% quest.AVdelay = AVdelay;
% InitializePsychSound;
% % Perform basic initialization of the sound driver:
% nrchannels = 1;
% 
% devices = PsychPortAudio('GetDevices')
% audiodevice = 6;
% pahandle   =  PsychPortAudio('Open',audiodevice, [], 1, fs, nrchannels);
% beephandle =  PsychPortAudio('Open',audiodevice, [], 1, fs, nrchannels);
% 
% % load memory buffers with auditory stimulus. Ready to fire!.
% PsychPortAudio('FillBuffer', pahandle, white_noise); % Auditory stimulus
% PsychPortAudio('FillBuffer',beephandle, beep_noise); % warning signal
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
cfg = quest.cfg

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


% I wont save the stimulus anymore. I simply calculate the sensitivities
% online


% Display instructions whereas stimuli are created
Detectinstructions;

% Wait until instructions have been readed
if MEG
    quest.bitsi.clearResponses;
    getResponse(quest.bitsi,1000,true); % seems that when it says false, even if you press a button it doesnt exit the function
else
    KbWait(); % press a key to continue with the script
end

%% Initialize quest parameters / only for the visual question experiment
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
        disp('loading previous session quest')
        q_center_intersession = Experimentalsession((session-1)).quest.q_center_intersession;
        q_surround_intersession = Experimentalsession((session-1)).quest.q_surround_intersession;
        % generate the new contrast for the new trial
        contrast_center =  QuestQuantile(q_center_intersession);
        contrast_surround  =  QuestQuantile(q_surround_intersession);
    end    
end
% thus q is the quest that it is applied in this run, and q_intersession is
% the quest across sessions.


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



% Opening file to save data
   
dataadaptoutput = [path, filesep 'results',filesep initials,'_session_',num2str(session),'_AVE_pedestal_meg_quest.txt'];


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
   
while ~exit_loop % exit the loop when trials are finish
    
    % Updating time
    
    absolute_time = GetSecs;
    
    time = absolute_time - init_time; % measures relative elapsed time
    
    switch state
        
        case 'Initialize'
         %   disp('Initialize')
            disp(iTrial)
            
            if (iTrial  >= quest.nTrials)
                exit_loop = 1; % ends the loop
                
            else % buffer stim and get time for next stimulus presentation
                
                iTrial = iTrial + 1;
                % Trial conditions

                % select center-surround combinations x modality x
                % signal/noise
                stimulus_definition = quest.stim_combinations(iTrial,1); % select orientation target
                stimulus_condition = quest.stim_conditions(stimulus_definition);
                modality_definition = quest.stim_combinations(iTrial,2); %
                modality_condition = quest.modality_conditions(modality_definition);
                target_space = quest.stim_combinations(iTrial,3);
                target_definition = quest.stim_combinations(iTrial,1); % orientation definition
                 
                 %% selecting whether the trials if an oddball or not
                % (in attention experiment there are oddball trials, not in the QUEST)
                oddball = quest.oddballs_vector(iTrial);
                % oddball = 1;
               
                if oddball
                    rand_oddball = randi([1 3]);
                    oddball_cfg =   quest.oddtions(rand_oddball,:);
                else
                    rand_oddball = 4; % if there is no oddball, name 4 the trial
                    oddball_cfg = [0 0];
                end
                
                
                tilt = pi/180*quest.orientations(stimulus_definition);
                
                cfg.gaborang = tilt; % Gabor orientation (rad)
                cfg.gaborphi = rand; % Gabor unit phase (randomize in each trial)
                
                if strcmp(stimulus_condition,'noise') % 1 is CCW, 2 is CW and 3 is noise
                    cfg.gaborcon = 0; % Gabor Michelson contrast
                    contrast_trial = 0;
                else
                    
                    % select contrast for center or surround
                    if  strcmp(quest.visual_presentation(target_space),'center')
                        
                        if contrast_center <0
                            contrast_center = 0; % contrast cant never be smaller than 0
                        end
                        cfg.gaborcon =   contrast_center; % Gabor contrast
                        
                        contrast_trial = contrast_center; % to save in the data matrix
                    else
                        if  contrast_surround <0
                            contrast_surround = 0; % contrast cant never be smaller than 0
                        end
                        cfg.gaborcon =   contrast_surround; % Gabor contras
                        contrast_trial = contrast_surround; 
                    end                    
                end
                
                % presenting center or surround
                if  strcmp(quest.visual_presentation(target_space),'center')
                    spatial_filter = mask(1).filter;
                    indexes_sf = find(spatial_filter == 0);
                else
                    spatial_filter = mask(2).filter;
                    indexes_sf = find(spatial_filter == 0);
                end
                                
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
                                
                quest.stimuli(iTrial).cfg = cfg;
                
                img = (gray * img)/0.5;
                noiseimg = (gray * noiseimg)/0.5;    % scale to rgb values
  
                               
                %% Here I prepare the targets texture
                % all the target textures were pre-generatedthe at the begining). Save the matrix always
                % for regression analyses
                
                % Creating textures
                gabortex = Screen('MakeTexture', win,   img); %Texture generated
                noisetex = Screen('MakeTexture', win,  noiseimg); %Texture noise generated
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
                keypressed = 'NoResponse';
                response_given = 'NoResponse';
                state =      'Fixationpoint';
                keyIsDown = 0;
                triggerpressed = 0;
                %quest.bitsi.clearResponses;
                empty_resp = 0;
                
            end

            %% SETTING UP TRIAL REFERENCE TIMINGS
            TRIAL_FP      = quest.TRIAL_TIMES(iTrial); % DRAWS FP
            TRIAL_1PED    = TRIAL_FP + quest.FIXATION(iTrial); % DRAW 1s PEDESTAL
            TRIAL_CUE1    = TRIAL_1PED + quest.CUE1(iTrial); % DRAW 1ST CUE
            
            TRIAL_TARGET_2PED  =  TRIAL_1PED + quest.PRE_STIM_TEMPLATE(iTrial); % DRAWS TARGET & PEDESTAL
            TRIAL_CUE2        = TRIAL_TARGET_2PED + quest.CUE2(iTrial); % DRAW TARGET
            TRIAL_RESPONSE    = TRIAL_TARGET_2PED + quest.POST_STIM_TEMPLATE(iTrial); %
            TRIAL_FEEDBACK    = TRIAL_RESPONSE + quest.RESPONSETIME(iTrial); %
            TRIAL_ENDS        = TRIAL_FEEDBACK + quest.FEEDBACK(iTrial);
            
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
                
                Screen('Flip', win);
                state = 'Pedestal_One';
            end
            
            
        case 'Pedestal_One'
             
            if (time >=  TRIAL_1PED) % DRAW VISUAL CUE
                
                if (time >=  TRIAL_CUE1  )
                    Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                    Screen('FillOval', win,[0 0 0], centeredRect_fp, radious_fp);
                    
                    if target_space == 1 % for center oddball
                        if oddball_cfg(1) % this will only be 1 if it is an oddball question block
                  %          Screen('FrameOval', win, [cue_color cue_color cue_color], in_VC_Rects, VC_thickness );
                  %          Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                        else
                  %          Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                  %          Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                        end
                    else % for surround  oddball
                        if oddball_cfg(1) % this will only be 1 if it is an oddball question block
                  %          Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                  %          Screen('FrameOval', win, [cue_color cue_color cue_color], out_VC_Rects, VC_thickness );
                        else
                  %          Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                  %          Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                        end
                        
                    end
                    
                    
                    Screen('DrawingFinished', win);
                    Screen('Flip', win);
                    Vprecue = 1; % avoids drawing the non cued stimulus
                end
                
                if (time >=  (TRIAL_CUE1 + quest.CUEDUR )) % remove the cue
                    Vprecue = 2 ; % open again the function to draw for a second time the normal pedestal
                end
                
                if (Vprecue == 0 | Vprecue == 2)
                    %  disp('pedestal ONE')
                    %   disp(time)
                    % Draw template
                    Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                    Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
                %   Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
                %   Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                    Screen('DrawingFinished', win);
                    Screen('Flip', win);
                    
                    if Vprecue == 2;
                        state = 'Play';
                    end % if this is the second time that the function draws the template, then it can go to play section
                    Vprecue = 1;
                end
                
            end
            
            
        case 'Play' % displays the target, the pedestal, and the second cue if necesary
           
            % imageArray = Screen('GetImage',win);
            % imwrite(imageArray, [path, filesep 'results' filesep 'pedestal1surround.jpg'])

            if (time >=  TRIAL_TARGET_2PED & Vflag == 0) % DISPLAYING VISUAL TARGET
                %disp('pedestal & target')
                %disp(time)
                Screen('DrawTextures', win, gabortex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
         %      Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
         %      Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                Screen('DrawingFinished', win);
                Vflag = 1; % Avoiding enter here and redraw the stimulus
                
                if time_calibration
                    Screen('FillOval', win, [0 0 0],centeredRect_calibration_cue, 400);
                end
                Screen('Flip', win);
            end
            
            if (time >= TRIAL_TARGET_2PED + quest.STIMDUR) %  JUMP TO NEXT STAGE TO ERASE VISUAL TARGET
                state = 'Pedestal_Two';
            end
                        
        case 'Pedestal_Two'
            
            if (time >=  TRIAL_TARGET_2PED)
                
                if (time >=  TRIAL_CUE2  )
                    Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                    Screen('FillOval', win,[0 0 0], centeredRect_fp, radious_fp);
                    
                     if target_space == 1 % for center oddball
                        if oddball_cfg(2) % this will only be 1 if it is an oddball question block
          %                 Screen('FrameOval', win, [cue_color cue_color cue_color], in_VC_Rects, VC_thickness );
          %                 Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                        else
          %                 Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
          %                 Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                        end
                    else % for surround  oddball
                        if oddball_cfg(2) % this will only be 1 if it is an oddball question block
          %                 Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
          %                 Screen('FrameOval', win, [cue_color cue_color cue_color], out_VC_Rects, VC_thickness );
                        else
          %                 Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
          %                 Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
                        end
                        
                    end
                                        
                    Screen('DrawingFinished', win);
                    Screen('Flip', win);
                    Vpostcue = 1; % avoids drawing the non cued stimulus
                end
                
                if (time >=  (TRIAL_CUE2 + quest.CUEDUR )) % remove the cue
                    Vpostcue = 2 ; % open again the function to draw for a second time the normal pedestal
                end
                
                if (Vpostcue == 0 | Vpostcue == 2)
                    % Draw template
                    Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                    Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
          %         Screen('FrameOval', win, [0 0 0], in_VC_Rects, VC_thickness );
          %         Screen('FrameOval', win, [0 0 0], out_VC_Rects, VC_thickness );
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
                    quest.bitsi.clearResponses;  % empty bitsi buffer
                    empty_resp = 1;
                end
                
                Screen('FillRect', win, gray);
                Screen('TextSize',  win, 3*ppd);
                
                DrawFormattedText(win, response_options(1), wMid - (5*ppd), 'center',  [0,0,0]); % right response option
                
                DrawFormattedText(win, response_options(2), wMid + (3*ppd), 'center',  [0,0,0]); % left response option
                Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
                
                Screen('Flip', win);
                
                % If the user press the key, record the response
                
                if ~keyIsDown % only enter if the characters has not been pressed. This allows to capture only one keypress
                    
                    [keyIsDown,secs,keycode] = KbCheck;
                    
                    % the option below is for mac when you have multiple
                    % keyboards (you have to specify what one)
                    % [keyIsDown,secs, keycode]= PsychHID('KbCheck', deviceNumber);
                    %RestrictKeysForKbCheck([KbName('ESCAPE'), KbName('SPACE'), KbName('UpArrow'),KbName('DownArrow')]); % select a group of keys to press. This
                    RestrictKeysForKbCheck([KbName('ESCAPE'), KbName('SPACE'), KbName('y'),KbName('u')]); % select a group of keys to press. This
                    
                    % WaitSecs(0.001); % delay to prevent CPU hogging

                      
                        [triggerpressed] = getResponse(quest.bitsi,0.01,true); %

                        if keyIsDown | triggerpressed > 0
                            % If the user holds down a key, KbCheck will report multiple events.
                            % To condense multiple 'keyDown' events into a single event, we wait until all
                            % keys have been released.
                            % KbReleaseWait;
                            response = 1; % responses confirmed
              
                            if MEG 
                               secs =  absolute_time; 
                               keyIsDown = 1;
                              % quest.bitsi.clearResponses;
                            end
                              
                            
                            RT = ((secs - init_time)) - TRIAL_RESPONSE; % RT after the response options are shown
                            %disp([(secs - init_time),RT])
                             
                             
                            % Note that we use find(keyCode) because keyCode is an array.
                            % See 'help KbCheck'
                            keypressed = KbName(keycode);
                            
                            % DONT FORGET TO REGISTER RTs
                            %disp(keypressed)
                            %% Recoding button/key presses to response option
                            if (strcmp( keypressed , 'y') |  triggerpressed == 98 |  triggerpressed == 66)
                                response_given =  response_options(1);% In this code, Y is left
                            end
                            if (strcmp( keypressed , 'u') | triggerpressed == 99 |  triggerpressed == 67)
                                response_given =  response_options(2);% In this code, U is right (you can use button box commands)
                            end                          

                            disp (quest.visual_presentation(target_space))
                            %% Cheking if the response is correct for the oddball task or visual question task
                            if strcmp(procedure, 'oddball')
                                n_oddballs = sum(oddball_cfg);
                                if (n_oddballs == str2num(response_given))
                                    correct = 1;
                                    %disp(tilt);
                                    disp('correct!')
                                else
                                    correct = 0;
                                    disp('wrong!')
                                end
                                
                            else
                                if ((strcmp(stimulus_condition,'signal')  & strcmp( response_given , 'Y')) | (strcmp(stimulus_condition,'noise')  & strcmp( response_given , 'N')))
                                    correct = 1;
                                    %disp(tilt);
                                    disp('correct!')
                                else
                                    correct = 0;
                                    disp('wrong!')
                                end
                                
                            end
                            
                        
                        if strcmp(procedure, 'detection')
                            if (strcmp(stimulus_condition,'signal')) % only update quest for signal trials
                                                 
                                if( strcmp(quest.visual_presentation(target_space),'center'))
                                    
                                    q_center_intersession = QuestUpdate(q_center_intersession,contrast_center,correct);
                                    % display the new contrast estimate
                                    disp(['intersessions new threshold = ', num2str(QuestQuantile(q_center_intersession))]);
                                    quest.q_center_intersession = q_center_intersession; % saving the quest
                                    contrast_center = QuestQuantile(q_center_intersession); 
                                    quest.q_center_intersession_contrast = contrast_center; % saving final stimated contrast value
                                else
                                    q_surround_intersession = QuestUpdate(q_surround_intersession,contrast_surround,correct);
                                    % display the new contrast estimate
                                    disp(['intersessions new threshold = ', num2str(QuestQuantile(q_surround_intersession))]);
                                    quest.q_surround_intersession = q_surround_intersession; % saving the quest
                                    contrast_surround = QuestQuantile(q_surround_intersession);
                                    quest.q_surround_intersession_contrast = contrast_surround; 
                                    
                                end
                                
                            end
                            
                        end
                                                  
                    end
                    
                    
                    %disp( keypressed)
                    
                    if strcmp(KbName(keycode), 'ESCAPE');
                        Screen('Closeall')
                        % closing screen and audio
                     %   PsychPortAudio('Close', pahandle)
                     %  PsychPortAudio('Close', beephandle)
                        %break;
                    end
                    
                end                      
                
                if (time >=  TRIAL_FEEDBACK)
             % disp('response')
              %     disp(time)
                    state = 'Feedback';   % jump to the provide feedback phase and records trial information
                    % s 1 session, 2 trial, 3stim time, 4 orientation,5postcue,6intensity, 7 resp, 8 RT, 9 correct
                    quest.responses_vector{iTrial,1} = session;
                    quest.responses_vector{iTrial,2} = iTrial;
                    quest.responses_vector{iTrial,3} = TRIAL_1PED; %% record postcue an all timing parameters for the MRI analysis
                    quest.responses_vector{iTrial,4} = TRIAL_CUE1; %% record postcue an all timing parameters for the MRI analysis
                    quest.responses_vector{iTrial,5} = TRIAL_TARGET_2PED; %% record postcue an all timing parameters for the MRI analysis
                    quest.responses_vector{iTrial,6} = TRIAL_CUE2; %% record postcue an all timing parameters for the MRI analysis
                    quest.responses_vector{iTrial,7} = TRIAL_RESPONSE; %% record postcue an all timing parameters for the MRI analysis
                    quest.responses_vector{iTrial,8} = char(stimulus_condition);
                    quest.responses_vector{iTrial,9} = char(modality_condition);
                    quest.responses_vector{iTrial,10} = char(quest.visual_presentation(target_space));
                    quest.responses_vector{iTrial,11} = rand_oddball; % if 1 is cue in 1 pedestal, 2 is cue in 2 pedestal, 3 is cue in both pedestals, 4 is nocue
                    quest.responses_vector{iTrial,12} = quest.orientations(target_definition);
                    quest.responses_vector{iTrial,13} =  contrast_trial; % saving contrast
                    quest.responses_vector{iTrial,14} = quest.cfg.freq; % saving contrast
                    quest.responses_vector{iTrial,15} = response_options; 
                    quest.responses_vector{iTrial,16} = response_given;
                    quest.responses_vector{iTrial,17} = triggerpressed;
                    quest.responses_vector{iTrial,18} = correct;
                    quest.responses_vector{iTrial,19} = RT;
                    quest.responses_vector{iTrial,20} = char(initials);
                    %    Experimentalsession(session).quest = quest;
                    %     save(dataFile,'Experimentalsession');
                  
                    % calculating the proportion of correct responses at this
                    % trial
                    correct_responses = mean([quest.responses_vector{1:iTrial,18}]);
                    %find indexes for center and surround events
                    center_stimuli = find(quest.stim_combinations(1:iTrial,3) == 1);
                    correct_center_responses = mean([quest.responses_vector{center_stimuli,18}]);
                    
                    surround_stimuli = find(quest.stim_combinations(1:iTrial,3) == 2);
                    correct_surround_responses = mean([quest.responses_vector{surround_stimuli,18}]);
                    
                    
                    disp(['proportion of TOTAL correct responses = ',num2str(correct_responses)]);
                    disp(['proportion of correct responses center = ',num2str(correct_center_responses)]);
                    disp(['proportion of correct responses surround = ',num2str(correct_surround_responses)]);
                    
                     disp(['intersessions center new threshold = ', num2str(QuestMean(q_center_intersession))]);
                     disp(['intersessions surround new threshold = ', num2str(QuestMean(q_surround_intersession))]);
                    % printint test data Do not save in the experiment. Slows
%% CONSIDER COMMENTING THIS BELOW
                    fileID = fopen([dataadaptoutput],'wt'); % try to write on top
                    formatSpec = '%i %i %f %f %f %f %f %s %s %s %i %f %f %f %s %s %i %i %f %s \r\n';
                     for row = 1 : nrows
                        fprintf(fileID, formatSpec, quest.responses_vector{row,:}) ;
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
                
          %      Screen('TextSize',  win, 2*ppd);
           %     keypressedis = ['triggerpress ', num2str(triggerpressed) ];
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
    
%     if (time >= (TRIAL_TARGET_2PED + AVdelay) & strcmp(modality_condition,'AV') & Aflag == 0) % DISPLAYING SOUND
%         % PsychPortAudio('Start', pahandle); % ensure aligment with oscilloscope - AVdelay variable
%         PsychPortAudio('Start', beephandle); % ensure aligment with oscilloscope - AVdelay variable        
%         Aflag = 1;
%     end
%     
end

  

Experimentalsession(session).quest = quest;

% lets plot the number of signal reports as a function of contrast

for visualfield_quest = 1 : 2
    figure;
    
    if visualfield_quest == 1
        selected_trials = center_stimuli;
    else
        selected_trials = surround_stimuli;
    end
    
    contrast_responses = Experimentalsession(session).quest.responses_vector(:,[13,16]);
    x = cell2mat(contrast_responses(selected_trials,1));
    y = strcmp(contrast_responses(selected_trials,2),'Y');
    fit = glmfit(x, y, 'binomial','link','probit');
    yfit = glmval(fit,x,'probit');
    plot(x,y,'o',x,yfit,'-');
    
end


save(dataFile,'Experimentalsession');

% save data in a txt file
% printint test data
                    
fileID = fopen([dataadaptoutput],'wt'); % try to write on top
 formatSpec = '%i %i %f %f %f %f %f %s %s %s %i %f %f %f %s %s %i %i %f %s \r\n';

for row = 1 : nrows
    fprintf(fileID, formatSpec, quest.responses_vector{row,:}) ;
end


fclose(fileID);


%Experimentalsession(session).quest = quest;
%save(dataFile,'Experimentalsession');

toc

Blockfinish;

    KbWait(); % press a key to continue with the script


Screen('Closeall')
% closing screen and audio
%PsychPortAudio('Close', pahandle)
%PsychPortAudio('Close',beephandle)

end
%fprintf('%5.2f	%4.1f	%5.2f\n',tActual,q.beta,q.gamma);


