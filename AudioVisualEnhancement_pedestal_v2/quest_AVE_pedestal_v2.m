function quest_AVE

% Quest function for AV integration experiment. This is the basic
% experiment without attentional manipulation. People has to discrimiante
% signal+noise from noise

% Provide feedback

tic
%% Switch to Screen('Preference', 'SkipSyncTests', 1) if trouble starting program
%initimeexp = GetSecs; % get time to measure the total lenght of the experiment
%% Set current directory
%path = '/Volumes/ALEX_EXT/Experimentos/DondersExperiments' %
path = mfilename('fullpath'); % Finds current location of file
[pathstr] = fileparts(path); % Removes file information, saves directory
path = pathstr;


%% Collect subject info
% prompt = {'Enter subject initials: ', 'Enter Age: ','Enter Gender(M/F): ', 'Guess contrast(0.1):','Guess sd(0.075):'};
% input = inputdlg(prompt);
% initials = input{1}; % Collect subject initials REMEMBER ATACCH ALL THE
% INFORMATION TO THE STRUCT VARIABLE
% Age = input{2};
% Gender = input{3};
% tGuess = input{4};
% tGuess = str2num(tGuess);
% tGuessSd = input{5};
% tGuessSd  = str2num(tGuessSd );

prompt = {'Enter subject ID initials: ', 'Guess center (0.007):', 'Guess surround (0.007):'};
input = inputdlg(prompt);
initials = input{1}; % Collect subject initials

center_tGuess = input{2};
center_tGuess = str2num(center_tGuess);
surround_tGuess = input{3};
surround_tGuess = str2num(surround_tGuess);
tGuessSd = 0.05;
%tGuessSd = input{5};
%tGuessSd  = str2num(tGuessSd);



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

Experimentalsession(session).initials = initials;
Experimentalsession(session).session = session;
Experimentalsession(session).date = date;



% creating the  subject file;
save(dataFile,'Experimentalsession');

init_AVE_parameters_pedestal_v2;

quest = init;

% If MAC, select the external device
%[keyboardIndices, productNames] = GetKeyboardIndices
%deviceNumber = 2;

% PREPARING STIMULI COMBINATIONS

reps =  10;
quest.stim_conditions = {'signal','noise'};
quest.visual_presentation = {'center','surround'};
% init parameters for the 1st pre adaptation test phase
[stim_rand_order,stim_combinations] = combinator([length(quest.orientations),length(quest.visual_presentation)],reps);

%randomize order
stim_combinations = stim_combinations(stim_rand_order,:);
quest.stim_combinations = stim_combinations;
quest.repetitions = reps;
quest.nTrials = length(stim_combinations);


%% Experiment timings
for t_idx = 1:quest.nTrials % consider randomizing a littlebit the ITI in each trial
    quest.TRIAL_TIMES(t_idx) = quest.INIT_DELAY +...
        (t_idx-1)*(quest.FIXATION + quest.PRE_STIM_CUE  +  quest.POST_STIM_CUE + quest.RESPONSETIME + quest.ITI); %
end

curr_sched =   quest.TRIAL_TIMES ; %Shift all the times by the INIT DELAY

%% RECORDING DATA IN THIS TEMPLATE MATRIX
quest.responses_vector = zeros(quest.nTrials,12); % vector to record conditions 1 session, 2 trial, 3stim time, 4 orientation,5postcue,6intensity, 7 resp, 8 RT, 9 correct
quest.responses_vector = num2cell(quest.responses_vector);
% these will be my variables
quest.responses_vector_variables = {'trial','session','targetpresentation','visualspace','grating_tilt','contrast','SF','keypressed','correct','RT'};

[nrows,ncols] = size(quest.responses_vector); % this makes very slow the process. Record the data at the end in txt file

%% setup keyboard basic parameters

KbName('UnifyKeyNames')

%will avoid subject pressing incorrect buttons. It can also improve
%computer performance

%%%%%%%%%%%%%%%%%%%% Video set up %%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check script is running on Psychtoolbox-3
%PsychDebugWindowConfiguration % transparent window
%Screen('Preference', 'SkipSyncTests', 1);
AssertOpenGL;

% Selects most peripheral screen
whichScreen = max(Screen('Screens'));
%rect=[0 0 640 480]; % mini-screen.
rect=[0 0 800 600]; % mini-screen.

HideCursor;

% Open graphics window on peripheral screen at certain resolution
% Our test screen in the lab will have dimensions 1024 x 768, so this is to
% test the look/feel on that size screen
% [win winRect] = Screen('OpenWindow', whichScreen,[],rect) %,rectx
[win winRect] = Screen('OpenWindow', whichScreen,[]) %,rectx

%HideCursor;
% Running on PTB-3? Abort otherwise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

% set noisy Gabor parameters
cfg = quest.cfg

% creating center and surround masks
patchsiz = cfg.patchsiz;
[rr cc] = meshgrid(1:patchsiz,1:patchsiz);
mask_radius = patchsiz/2;
                
mask(1).filter = sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)<=mask_radius/2;  %center in pixels
mask(1).filter = mask(1).filter - (sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)<=mask_radius/20); % add inner empty area
mask(2).filter = sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)>=mask_radius/2;  %surround in pixels
             
% Instructions
newQuestinstructions;
KbWait(); % press a key to continue with the script

%% Initialize quest parameters

pThreshold=0.72; % (previous was 0.75 and gamma 0.5)
beta=3.5;delta=0.01;
gamma=0.5; % gamma is the proportion correct of responses that are expected when stimulus is not present. slope, threshold and minimum value (it is a 3AFC, noise, 2 x tilt )
% the quest is designed only to update in noise and signal trials
% Si el sujeto no ve nada, juzgara aleatoriamente, asi que puede acertara
% el 50% de las veces
see_quest = 0;
range = 0.1; % range below and above threshold to be tested
grain = 0.001;

% lets initialize one quest for the center and another for the surround
qcenter =                  QuestCreate(center_tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range,see_quest);
qcenter.normalizePdf =      1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.
qcenter_tTest =             QuestQuantile(qcenter); % assign a first tTest value

qsurround =                 QuestCreate(surround_tGuess,tGuessSd,pThreshold,beta,delta,gamma,grain,range,see_quest);
qsurround.normalizePdf =    1; % This adds a few ms per call to QuestUpdate, but otherwise the pdf will underflow after about 1000 trials.
qsurround_tTest =           QuestQuantile(qsurround); % assign a first tTest value

%horizontal_cue =  [(w/2)-16, (h/2),(w/2)+16,(h/2)];
%leftward_cue = [(w/2)-11, (h/2)+11,(w/2)+11,(h/2)-11];
%rightward_cue = [ (w/2)-11, (h/2)-11,(w/2)+11,(h/2)+11];

% POSITION COORDINATES FOR CIRCULAR OBJECTS(FP AND INNER ANNULUS)

color_fp = [0 0 0];
radious_fp = 20;
baseRect_fp = [0 0 radious_fp  radious_fp ]; % create a fake box to set the center of ovals
centeredRect_fp = CenterRectOnPointd(baseRect_fp, w/2, h/2);

radious_center = 50;
baseRect_center = [0 0 radious_center  radious_center ]; % create a fake box to set the center of ovals
centeredRect_center = CenterRectOnPointd(baseRect_center, w/2, h/2);


% Initialize loop parameters
iTrial = 0;
state = 'Initialize';
exit_loop = 0;

Screen('FillRect', win, gray);    % Gray background
vbl = Screen('Flip', win);
ifi = Screen('GetFlipInterval', win);

% use this line to draw and avoid tearing effects ->  vbl = Screen('Flip', win, vbl + 0.5 * ifi);

init_time = GetSecs;

while ~exit_loop % exit the loop when trials are finish
    
    % Updating time
    
    relative_time = GetSecs;
    
    time = relative_time - init_time; % measure absolute elapsed time
    
    switch state
        
        case 'Initialize'
            
            disp(iTrial)
            
            if (iTrial  >= quest.nTrials)
                exit_loop = 1;
                
            else % buffer stim and get time for next stimulus presentation
                
                iTrial = iTrial + 1;
                
                % Trial conditions
                
               target_definition = quest.stim_combinations(iTrial); 
                                
               target_space = quest.stim_combinations(iTrial,2); % center surround
    
                %%
                if  strcmp(quest.visual_presentation(target_space),'center')
                    spatial_filter = mask(1).filter;
                    indexes_sf = find(spatial_filter == 0);
                else
                    spatial_filter = mask(2).filter;
                    indexes_sf = find(spatial_filter == 0);
                end
    
           %     if strcmp(quest.stim_conditions(target_definition),'noise') % 1 is CCW, 2 is CW and 3 and 4 is noise
           %        ngratings = 0;             
           %     else
           %      ngratings = 1;
                
                
               %% Updating quests values
                    qsurround_tTest =    QuestQuantile(qsurround);	% Recommended by Pelli (1987), and still our favorite.
                    qcenter_tTest   =    QuestQuantile(qcenter);	% Recommended by Pelli (1987), and still our favorite.               
                    % 	tTest=QuestMean(q);		% Recommended by King-Smith et al. (1994)
                    % 	tTest=QuestMode(q);		% Recommended by Watson & Pelli (1983)
                    
                    if qsurround_tTest < 0    % restrict quest to positive values
                        qsurround_tTest = 0;
                    end
       
                    if qcenter_tTest < 0    % restrict quest to positive values
                        qcenter_tTest = 0;
                    end
                   
              %  end
                
                % Select trial timings
                Fixation_point_time = curr_sched(iTrial);
                
                Pre_stim_cue_time = curr_sched(iTrial) + quest.FIXATION;
                Pre_stim_cue_time_off = Pre_stim_cue_time + quest.CUEDUR; % switch off the cue
                
                visual_stimulus_jitter = quest.jitter_range + (-quest.jitter_range-quest.jitter_range).*rand(1,1); % visual stimulus can be jittered between -.2 and 0.2 s
                
                V_stimulus_time = Pre_stim_cue_time + quest.PRE_STIM_CUE ; % this is the time where the visual stimulus should be presented in this trial
               
                 %adding jitter in a trial by trial basis 
                real_stimulus_time = V_stimulus_time + visual_stimulus_jitter;
                real_stimulus_time_off =  real_stimulus_time + quest.STIMDUR; % this is the time where the visual stimulus should be removel
                
                Post_stim_cue_time = V_stimulus_time + quest.POST_STIM_CUE; % after postcue, response should be collected
                Post_stim_cue_time_off = Post_stim_cue_time + quest.CUEDUR; % switch off the cue
                
                
                Feedback_time = Post_stim_cue_time + quest.RESPONSETIME; % provide feedback about correctness
                
                
                if strcmp(quest.stim_conditions(target_definition),'noise')
                    
                    contrast = 0 ;
                else  % select contrast for center or surround (depending on the type of trial)
                    if strcmp(quest.visual_presentation(target_space),'center')
                        contrast =  qcenter_tTest;             
                    else
                        contrast =  qsurround_tTest;
                    end
                end

                
                % here draw MASK texture
                tilt = pi/180*quest.orientations(target_definition);
                cfg.gaborang = tilt; % Gabor orientation (rad)
                cfg.gaborphi = rand; % pi/180*90; % Gabor unit phase
                cfg.gaborcon = 0; % Gabor Michelson contrast
                
                % clip patch luminance to [black,white] range to avoid black clusters at
                % high contrasts
                
                while true
                    img = make_gabor_and_noise(cfg);
                    %[patchimg] = absoluteDifferenceBetweenWhiteAndGray * patchimg;
                    img = (gray * img)/0.5;
                    
                    if all(img(:) > black & img(:) < white)
                        break
                    end
                end
               
                masktext = Screen('MakeTexture', win,  img); %Texture generated
                
                % AND BELOW I DRAW THE TARGET TEXTURE WITH A PARTICULAR
                % CONTRAST LEVEL SELECTED BY QUEST
                
                cfg.gaborphi = rand; % pi/180*90; % Gabor unit phase CHECK THIS (radians?)
                cfg.gaborcon = contrast; % Gabor Michelson contrast
                
                % clip patch luminance to [black,white] range to avoid black clusters at
                % high contrasts
                
                while true
                    [img noiseimg] = make_gabor_and_noise(cfg);
                    % filtering center or surround to pedestal and
                    % stimulus
                     img(indexes_sf ) = 0.5;
                     noiseimg(indexes_sf) = 0.5;
                    %[patchimg] = absoluteDifferenceBetweenWhiteAndGray * patchimg;
                    img = (gray * img)/0.5;
                    noiseimg = (gray * noiseimg)/0.5;
                    
                    if all(img(:) > black & img(:) < white)
                        break
                    end
                end
                
                % save all the information about the visual stimulus
                quest.Grating(iTrial).cfg = cfg;
              % quest.Grating(iTrial).img = img; % i dont save this because
              % the file gets too heavy
                
                gabortex = Screen('MakeTexture', win,  img); %Texture generated
                noisetex = Screen('MakeTexture', win,  noiseimg); %Texture generated
                % imagesc(n)
                texrect = Screen('Rect', gabortex); %Extract information about texture size
                scale = 1; % scaling artifitially the stimulus; Use carefully. Default 1
                posx = [w/2];   %coordinates
                posy = [h/2];
                poscenter = [w/2, h/2];
                rotAngles = [0]; % we create the gratings already rotated. This can be applied to textures
                
                dstRects = CenterRectOnPoint(texrect*scale , posx, posy);
                
                %Initializing flags / We only want the functions to draw
                %each item one time to avoid delays
                Vflag = 0;
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
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                state = 'Precue';
            end
            
            
        case 'Precue'
            
            if (time >=  Pre_stim_cue_time & Vprecue == 0) % DRAW VISUAL CUE
                
                %DrawFormattedText(win,'-', 'center','center', [240,240,240]); % in quest only one type of cue
                % Screen('DrawLine', win, [255,255,255], (w/2)-16, (h/2),(w/2)+16,(h/2), 6);
                % Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround                
                %DrawFormattedText(win,'-', 'center','center', [240,240,240]); % in quest only one type of cue
                % Screen('DrawLine', win, [255,255,255], (w/2)-16, (h/2),(w/2)+16,(h/2), 6);
                Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
                Screen('DrawingFinished', win);
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                Vprecue = 1;
                state = 'Play';
            end
%             
%             if (time >=  Pre_stim_cue_time_off) % ERASE VISUAL CUE,RESTORE NORMAL FIXATION AND MOVE TO STIM PLAYING.
%                 Screen('TextSize',  win, 30);
%                % Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
%                 Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround                
%                 %DrawFormattedText(win,'-', 'center','center', [240,240,240]); % in quest only one type of cue
%                 % Screen('DrawLine', win, [255,255,255], (w/2)-16, (h/2),(w/2)+16,(h/2), 6);
%                 Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
%                 Screen('DrawingFinished', win);
%                 vbl = Screen('Flip', win, vbl + 0.5 * ifi);
%                 state = 'Play';
%             end
%             
            
        case 'Play'
            
            if (time >= real_stimulus_time & Vflag == 0) % DISPLAYING VISUAL TARGET
                
                Screen('DrawTextures', win, gabortex,[],dstRects,rotAngles); %, 2 is the adaptor surround
                Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                Screen('DrawingFinished', win);
                Vflag = 1; % Avoiding enter here and redraw the stimulus
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
            end
            
            if (time >= real_stimulus_time_off) % ERASE VISUAL TARGET and REDRAW NOISE PATCH
                 % Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                Screen('DrawTextures', win, noisetex,[],dstRects,rotAngles); %, 2 is the adaptor surround                
                %DrawFormattedText(win,'-', 'center','center', [240,240,240]); % in quest only one type of cue
                % Screen('DrawLine', win, [255,255,255], (w/2)-16, (h/2),(w/2)+16,(h/2), 6);
                Screen('FillOval', win,color_fp, centeredRect_fp, radious_fp); %color_fp
                Screen('DrawingFinished', win);
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                state = 'Postcue';
                
            end
            
        case 'Postcue'
            
            if (time >=  Post_stim_cue_time & Vpostcue == 0) % If both stimuli displayed, move to the next phase
                
                %DrawFormattedText(win,cue_sign{post_cue}, 'center','center', [240,240,240]);
                %    if post_cue_orientation == 45
                %       Screen('DrawLine', win, [255,255,255], (w/2)-11, (h/2)+11,(w/2)+11,(h/2)-11, 6);
                %               end
                %    if post_cue_orientation == 315
                %        Screen('DrawLine', win, [255,255,255], (w/2)-11, (h/2)-11,(w/2)+11,(h/2)+11, 6);
                %         end
                  Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                
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
                    % RestrictKeysForKbCheck([KbName('ESCAPE'), KbName('SPACE'), KbName('UpArrow'),KbName('DownArrow')]); % select a group of keys to press. This
                    RestrictKeysForKbCheck([KbName('ESCAPE'), KbName('SPACE'), KbName('y'),KbName('u')]); % select a group of keys to press. This
                    
                    WaitSecs(0.001); % delay to prevent CPU hogging
                    
                    if keyIsDown
                        % If the user holds down a key, KbCheck will report multiple events.
                        % To condense multiple 'keyDown' events into a single event, we wait until all
                        % keys have been released.
                        %KbReleaseWait;
                        response = 1; % responses confirmed
                        
                        RT = ( (secs - init_time)) - Post_stim_cue_time;
                        % disp([(secs - init_time),RT])
                        
                        % Note that we use find(keyCode) because keyCode is an array.
                        % See 'help KbCheck'
                        keypressed = KbName(keycode);
                        % DONT FORGET TO REGISTER RTs
                        % correct or incorrect?
                        
                        if (strcmp(quest.stim_conditions(target_definition),'signal')  & strcmp( keypressed , 'y') | strcmp(quest.stim_conditions(target_definition),'noise')  & strcmp(  keypressed , 'u'))
                            correct = 1;
                            %  disp(tilt);disp(post_cue_orientation);
                            disp('correct!')
                        else
                            correct = 0;
                            disp('wrong!')
                        end
                        
                        %if strcmp(quest.stim_conditions(target_definition),'signal') % only update contrast if there was signal
                        
                        if strcmp(quest.visual_presentation(target_space),'center')
                            
                            qcenter = QuestUpdate(qcenter,qcenter_tTest,correct); % Add the new datum (actual test intensity and observer response) to the database.
                        else
                            qsurround = QuestUpdate(qsurround,qsurround_tTest,correct); % Add the new datum (actual test intensity and observer response) to the database.
                            
                        end
                               
                        quest.questparameters.qcenter = qcenter; % save parameters
                        quest.questparameters.qsurround = qsurround; % save parameters
                                                                
                        disp(['center threshold ' num2str(qcenter_tTest)]);
                        disp(['surround threshold ' num2str(qsurround_tTest)]);
                    end
                    
                    %disp( keypressed)
                    
                    if strcmp(KbName(keycode), 'ESCAPE');
                        Screen('Closeall')
                        break;
                    end
                    
                end
            end
           % ADD TWO QUEST AND ALL THE DIFFERENT INFORMATION TO BE SAVED
            if (time >  Feedback_time)
                state = 'Feedback';   % jump to the provide feedback phase and records trial information
                % s 1 session, 2 trial, 3stim time, 4 orientation,5postcue,6intensity, 7 resp, 8 RT, 9 correct
                quest.responses_vector{iTrial,1} = session;
                quest.responses_vector{iTrial, 2} = iTrial;
                quest.responses_vector{iTrial,3} = real_stimulus_time; %% record postcue an all timing parameters for the MRI analysis
                quest.responses_vector{iTrial,4} = quest.orientations(target_definition);
                quest.responses_vector{iTrial,5} = char((quest.stim_conditions(target_definition)));
                quest.responses_vector{iTrial,6} = char(quest.visual_presentation(target_space));
                quest.responses_vector{iTrial,7} = quest.cfg.freq; % saving contrast
                quest.responses_vector{iTrial,8} = contrast;
                quest.responses_vector{iTrial,9} = keypressed;
                quest.responses_vector{iTrial,10} = correct;
                quest.responses_vector{iTrial,11} = RT;
                quest.responses_vector{iTrial,12} = char(initials);
                %    Experimentalsession(session).quest = quest;
                %     save(dataFile,'Experimentalsession');
               
                dataadaptoutput = [path, filesep 'results' filesep ,initials,'_session_',num2str(session),'_AVE_quest.txt'];
                % printint test data
                %% saving the data everytrial slowsdown the program (possibly because it saves the data in the net)
            %    fileID = fopen([dataadaptoutput],'wt'); % try to write on top
             %   formatSpec = '%i %i %f %f %s %f %s %i %f %s \r\n';
                % for windows
                %dataadaptoutput = [path,'\results\',initials,'_Quest.txt'];
                
            %    for row = 1 : nrows
            %        fprintf(fileID, formatSpec, quest.responses_vector{row,:}) ;
            %    end
                
             %   fclose(fileID);
                
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
                
                FeedBack_size = [0, 0, 30, 30];
                square_rect = CenterRect(FeedBack_size, winRect);
                % Screen('FillRect',win, color, square_rect);
                % drawing mask and other elements
                Screen('DrawTextures', win, masktext,[],dstRects,rotAngles); %, 2 is the adaptor surround
                Screen('FrameRect', win, color, square_rect, 2);
                Screen('FillOval', win, [0,0,255], centeredRect_fp, radious_fp);
               Screen('DrawingFinished', win);
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                
                state = 'Initialize'; % lets prepare the next trial (going backwards)
                
                
            end
            
            % Provide our prior knowledge to QuestCreate, and receive the data struct "q".
            
            % plot last guess and sd of this subject
    end
    
    %count_frames = count_frames + 1;
end

% save data in a txt file
dataadaptoutput = [path, filesep 'results' filesep ,initials,'_session_',num2str(session),'_AVE_quest.txt'];
% printint test data
fileID = fopen([dataadaptoutput],'wt'); % try to write on top
 formatSpec = '%i %i %f %f %s %s %f %f %s %i %f %s \r\n';
% for windows
%dataadaptoutput = [path,'\results\',initials,'_Quest.txt'];

for row = 1 : nrows
    fprintf(fileID, formatSpec, quest.responses_vector{row,:}) ;
end

fclose(fileID);

t_center = QuestMean(qcenter)		% Recommended by Pelli (1989) and King-Smith et al. (1994). Still our favorite.
t_center_sd = QuestSd(qcenter);

t_surround = QuestMean(qsurround)		% Recommended by Pelli (1989) and King-Smith et al. (1994). Still our favorite.
t_surround_sd = QuestSd(qsurround);

quest.questparameters.center_threshold = t_center;
quest.questparameters.center_sd = t_center_sd;

quest.questparameters.surround_threshold = t_surround;
quest.questparameters.surround_sd = t_surround_sd;
Experimentalsession(session).quest = quest;

save(dataFile,'Experimentalsession');

disp(['Quest CENTER estimated threshold = ',num2str(t_center)])
disp(['Quest SURROUND estimated threshold = ',num2str(t_surround)])
%avgfps = count_frames / (vbl - tstart);
toc
Questfinish; 
KbWait();
Screen('Closeall')

end

%fprintf('%5.2f	%4.1f	%5.2f\n',tActual,q.beta,q.gamma);
