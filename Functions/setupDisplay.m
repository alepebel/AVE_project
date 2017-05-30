
function ds = setupDisplay(psychMode)
% Set Display parameters

if ~exist('psychMode', 'var'), psychMode = false; end 
PsychTweak('UseGPUIndex', 1);
AssertOpenGL;

% Get the list of Screens and choose the one with the highest screen
% number, as a best guess about where you want the stimulus displayed

% Uncomment for development purposes:
% Screen('Preference', 'SkipSyncTests',1);

% Calculate screen size
RTOD = 180/pi;
% Screen width for the projector and the small in-bore screen
if psychMode
  ds.VIEWINGDISTANCE = 65; % Typical view distance from a screen in cm
  ds.SCREENXCM = 51.2; % Size of an iMac 27"
else
  % In-bore scanner screen dimensions:
  ds.VIEWINGDISTANCE = 49;
  ds.SCREENXCM = 32.5;
end


ds.scrnNum = max(Screen('Screens')); 
ds.scrnNum = 0; %1;




ds.SCREENXDEG = RTOD*atan(ds.SCREENXCM/ds.VIEWINGDISTANCE);
ds.refreshFreq = 60;

% Define some special colors:
ds.white = WhiteIndex(ds.scrnNum);
ds.black = BlackIndex(ds.scrnNum);
ds.gray = (ds.white+ds.black)/2;

% Open double-buffered onscreen window with the requested stereo mode:
[ds.windowPtr, ds.windowRect] = Screen('OpenWindow', ds.scrnNum, ds.gray); %remove [0 0 640 480] for f 
ds.PIXPERDEG = ds.windowRect(3)/ds.SCREENXDEG;

ds.ltheta = 0.00*pi; % Screen rotation to adjust for mirrors
ds.scr_rot = 0; % Screen Rotation for opponency conditions

% Set up alpha-blending for smooth (anti-aliased) drawing of dots:
Screen('BlendFunction', ds.windowPtr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
Screen('TextSize', ds.windowPtr, 40);
% glEnable(GL.DEPTH_TEST);
%Screen('DrawText', ds.windowPtr, 'Press Any Key to Continue', 220, 220, ds.black);  

ds.t = Screen('Flip', ds.windowPtr);
 %disp('flag')

% Start audio:
InitializePsychSound;

end