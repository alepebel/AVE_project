
%% setup keyboard basic parameters

KbName('UnifyKeyNames')

%% Video set up %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check script is running on Psychtoolbox-3
%Screen('Preference', 'SkipSyncTests', 1);
AssertOpenGL;
%HideCursor;
% Selects most peripheral screen
whichScreen = max(Screen('Screens'));
%rect=[0 0 640 480]; % mini-screen.
rect=[0 0 320 240]; % mini-screen.

% Removed
%HideCursor;

% Open graphics window on peripheral screen at certain resolution
% Our test screen in the lab will have dimensions 1024 x 768, so this is to
% test the look/feel on that size screen


[win winRect] = Screen('OpenWindow', whichScreen,[],rect) %,rectx

flipTime = Screen('Flip',win)

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
absoluteDifferenceBetweenWhiteAndGray = abs(white - gray);

%% Auditory stimuli%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_STIMDUR  = 0.04; %sec.
fs = 44100; %sampling fq
% creating the whitenoise
white_noise = (rand(1,A_STIMDUR*fs )*2-1); % LETS create some whitenoise.
[white_noise] = linramp(white_noise,A_STIMDUR,fs,0.01); % apply ramp on/off

InitializePsychSound
% Perform basic initialization of the sound driver:
nrchannels = 1;
pahandle = PsychPortAudio('Open', [], [], 1, fs, nrchannels);
% load memory buffer with stimulus. Ready to fire.
PsychPortAudio('FillBuffer', pahandle, white_noise );


%% Visual stimuli

sf = [2];  %1st column comparison. 2nd adaptor. 3rd comparisonadapted
tilt = [0]; %rename later
contrast = [0.1];
phase = [0];
rotAngles = [0];
radius = [500]; % i will apply the circular mask later to create the Composed grating. Use maximal radius now to avoid not including some pixels in the draw
posx = [w/2];   %coordinates x gabor 1, xgabor 2 and ygabor1&2 (divide W/4 and H/2)
posy = [h/2];    %the adaptor should be larger

n = gabor_creator(sf, tilt, contrast, black, white, 1, radius, phase); %gabor_creator(freq, tilt, contrast, black, white, gabor,radiusinpixels,[hase)
sigma = 1.;

wn =  0 + sigma.*randn(301,301);
wn = imgaussfilt(wn,2);
n = wn + n;
texture = gray + absoluteDifferenceBetweenWhiteAndGray * n;

gabortex=Screen('MakeTexture', win, texture); %Texture generated
% imagesc(n)
texrect = Screen('Rect', gabortex); %Extract information about texture size
scale = 1;
dstRects = CenterRectOnPoint(texrect*scale , posx, posy);

Screen('FillRect', win, gray);
Screen('DrawTextures', win, gabortex,[],dstRects,rotAngles); %, 2 is the adaptor surround

Screen('DrawingFinished', win);

vbl = Screen('Flip', win);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    PsychPortAudio('Start', pahandle);
    Screen(win,'Flip');
    
    PsychPortAudio('Close', pahandle)
    Screen('closeall')
