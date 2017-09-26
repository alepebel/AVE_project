function [Tm Tdb] = REGONS_sensation_level()

% [Tm Tdb] = ptb_sensation_level(Sound,dbRange,nRep)
%
% Sound   - Sound vector; by default a 12 s, 1000 Hz sine tone is used
% Sf      - sampling frequency; default = 44100
% dbRange - the minimun and maximum decibel applied to the Sound; default = [-150 -70]
% nRep    - number of repetitions per descending and ascending Sound; default = 6
%
% Tm  - threshold in ratio/multiplier units; use this to multiply your Sound with
% Tdb - threshold in decibel units
%
% Description: The script calculates a hearing threshold (sensation level)
% for a given Sound. It applies a method of limits whereby the Sound
% descends in intensity until inaudible, and then ascends until audible.
% -----------------------------------------------------------------------------------
% B. Herrmann, Email: bjoern.herrmann@outlook.com, 2015-02-17

% get all the defaults going
% if nargin < 4, nRep    = 6; end
% % if nargin < 3, dbRange = [-150 -70]; end
% if nargin < 3, dbRange = [-120 -30]; end
% if nargin < 2, Sf = 44100; end
% if isempty(Sound)
nRep = 6;
dbRange = [-120 -30];
Sf = 44100;
centerfreq = [1000 1000];
freqrange  = [0 0];
ncomps = [30 30];
eventdur = [12 12];
rampdur  = .005;
lbs = centerfreq(1) - (freqrange(1) / 2); % lower bound
ubs = centerfreq(1) + (freqrange(1) / 2); % upper bound
compfreqs = (lbs + (ubs - lbs)*rand(1,ncomps(1)));
amps = 1 - (abs((compfreqs-centerfreq(1))/centerfreq(1)));
Sound = getstim(compfreqs,amps,eventdur(1),rampdur,Sf);


%% MEG or Behavioral?
MEG = 1;

if MEG
    sound_calib.bitsi = Bitsi('com1');    
else
    sound_calib.bitsi = Bitsi('');
end


%% Set current directory and ID of participant

path = mfilename('fullpath'); % Finds current location of file
[pathstr] = fileparts(path); % Removes file information, saves directory
path = pathstr;
% jump one level below
cd(path);
cd('..');
path =  pwd;
%% Collect subject info

% universal file separator across operative systems is filesep


prompt = {'Enter subject initials:'};
input = inputdlg(prompt);
initials = input{1}; % Collect subject initials REMEMBER ATACCH ALL THE

outputFolder = [path, filesep 'results' filesep ,initials,filesep];

if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

dataFile = [outputFolder initials,'_sound_calib.mat'];


% get Sound duration
% Sound = -1 + (1-(-1)).*rand(1,length(Sound));
Sound = Sound(:);
durSound = length(Sound)/Sf;

% get db multipliers
dbLin  = linspace(dbRange(1),dbRange(2),length(Sound));  % linear spacing between decibels
m_asc  = 10.^(dbLin/20);                                 % ascending multipliers
m_des  = 10.^(fliplr(dbLin)/20);                         % descending multipliers
db_asc = dbLin;                                          % ascending decibel
db_des = fliplr(dbLin);                                  % descending decibel

% get ascending and descending Sound
asound = Sound .* m_asc';
dsound = Sound .* m_des';

% initial params
PsychDefaultSetup(2);
InitializePsychSound;
pahandle = PsychPortAudio('Open',[],[],0,Sf,1);

sca;                                             % close all screens
screens      = Screen('Screens');                % get all screens available
screenNumber = max(screens);                     % present stuff on the last screen
Screen('Preference', 'SkipSyncTests', 1);
white        = WhiteIndex(screenNumber);         % get white given the screen used

% play once something to get the time stamps right later
PsychPortAudio('FillBuffer',pahandle,0);
initialTime = PsychPortAudio('Start',pahandle,1,0,1);

try
	bgcolor  = 0;                                    % background color 0-255
	txtcolor = round(white*0.8);                     % text color 0-255
	[shandle, windowRect] = PsychImaging('OpenWindow', screenNumber, bgcolor); % open window
	
	% Get the center and size of the on screen window
	[xCenter, yCenter] = RectCenter(windowRect);
	[screenXpixels, screenYpixels] = Screen('WindowSize', shandle);
	a_up     = arrow_matrix(100);
	a_down   = flipud(arrow_matrix(100));
	tmp      = [min([screenXpixels screenYpixels]) min([screenXpixels screenYpixels])] * 0.2;
	DestRect = CenterRectOnPointd([0 0 tmp],xCenter, yCenter);
	
	% Select specific text font, style and size:
	Screen('TextFont', shandle, 'courier new');
	Screen('TextSize', shandle, 24);
	Screen('TextStyle', shandle, 1);
	HideCursor;
	
	% Instruction screen
	arrows = [a_down zeros(size(a_up)) a_up];
	arrows = [zeros([size(a_down,1) size(arrows,2)]); arrows; zeros([size(a_down,1) size(arrows,2)])];
	tmp    = [min([screenXpixels screenYpixels]) min([screenXpixels screenYpixels])] * 0.6;
	DRect  = CenterRectOnPointd([0 0 tmp],xCenter, yCenter);
	tex = Screen('MakeTexture', shandle, arrows);
	Screen('DrawTextures', shandle, tex, [], DRect, [], [], [], [0 1 0]);
	DrawFormattedText(shandle, [['Downward arrow: The Sound starts loudly,\npress a button as soon as you cannot hear it anymore!\n\n'] ...
	                           ['Upward arrow: The Sound starts silently,\npress a button as soon as you can hear it!']], 'center', yCenter*0.3+yCenter, txtcolor);
	DrawFormattedText(shandle, 'Press a button to start!', 'center', yCenter*0.7+yCenter, txtcolor);
	Screen('Flip', shandle);
    
    if MEG
        sound_calib.bitsi.clearResponses;
        getResponse(sound_calib.bitsi,1000,true);
    else
        keyIsDown = KbCheck;
        while ~keyIsDown, keyIsDown = KbCheck; end
    end
    
    
    DrawFormattedText(shandle, '', 'center', 'center', txtcolor);
 	Screen('Flip', shandle);
	WaitSecs(2);
	
	[Tm Tdb] = deal([]);
 	for rr = 1 : nRep
 		for pp = 1 : 2  % pp = 1 --> descending; pp = 2 --> ascending
			if pp == 1
 				stim  = dsound;
 				arrow = a_down;
 			else
				stim  = asound;
				arrow = a_up;
			end
			
			% draw arrow
			tex = Screen('MakeTexture', shandle, arrow);
			Screen('DrawTextures', shandle, tex, [], DestRect, [], [], [], [0 1 0]);
			Screen('Flip', shandle);
			
			% play Sound
 			PsychPortAudio('FillBuffer',pahandle,stim');
 			starttime = PsychPortAudio('Start',pahandle,1,0,1);
            if MEG                                
                sound_calib.bitsi.clearResponses
                [r RT] =  getResponse(sound_calib.bitsi,1000,true)
            else
                
                [keyIsDown, endtime, keyCode] = KbCheck;
                
                while ~keyIsDown && GetSecs-starttime < durSound
                    [keyIsDown, endtime, keyCode] = KbCheck;
                end
                
            end
            
			PsychPortAudio('Stop', pahandle,[],1);
			DrawFormattedText(shandle, '', 'center', 'center', txtcolor);
			Screen('Flip', shandle);
			
			% get reaction time
            if ~MEG
			RT = endtime - starttime;
            end
            
			RTsamp = round(RT * Sf);
			if RTsamp > length(m_asc) || RT < 1        % exclude responses made too fast or late
				[Tm(rr,pp) Tdb(rr,pp)] = deal(NaN);
			else
				if pp == 1
					Tm(rr,pp)  = m_des(RTsamp);
					Tdb(rr,pp) = db_des(RTsamp);
				else
					Tm(rr,pp)  = m_asc(RTsamp);
					Tdb(rr,pp) = db_asc(RTsamp);
				end
			end
			
			% wait befor next tria;
			WaitSecs(1.5);
		end
	end
	sca;
catch
	sca;
	rethrow(lasterror);
end

if nnz(isnan(Tm(:))) == length(Tm)
	error('Error: No valid responses were made! Responses were either too early or too late.')
	[Tm Tdb] = deal([]);
elseif nnz(isnan(Tm(:)))/length(Tm) > 0.5
	disp('Info: More than half of the responses were invalid! Responses were either too early or too late.')
	Tm  = nanmean(Tm(:));
	Tdb = nanmean(Tdb(:));
    sound_calibration.Tm = Tm;
    sound_calibration.Tdb = Tdb
    save(dataFile,'sound_calibration');

else
	Tm  = nanmean(Tm(:));
	Tdb = nanmean(Tdb(:));
    sound_calibration.Tm = Tm;
    sound_calibration.Tdb = Tdb
    save(dataFile,'sound_calibration');

end



% get sine tone
function [y] = sine_tone(dur,Cf,Sf)
t = 0:1/Sf:(dur-1/Sf);
y = sin(2*pi*Cf*t);
	

% get rise and fall ramps
function [y] = wav_risefall(x,rf,Sf)
% get samples for rise and fall times
nsamp_rise = round(rf(1)*Sf);
nsamp_fall = round(rf(2)*Sf);

% get rise and fall vectors
rise = linspace(0,1,nsamp_rise)';
fall = linspace(1,0,nsamp_fall)';

% applied rise and fall vectors to x
x = x(:);
x(1:nsamp_rise) = x(1:nsamp_rise).*rise;
x(end-nsamp_fall+1:end) = x(end-nsamp_fall+1:end).*fall;
y = x;


% get the arrow image
function a = arrow_matrix(n)
a = [fliplr(tril(ones([n n]))) tril(ones([n n]))];
a = a(1:n/2,n/2:(n*2-n/2)-1);
a = [a; [zeros([n n/4]) ones([n n/2]) zeros([n n/4])]];
a = [zeros([n+n/2 n/4]) a zeros([n+n/2 n/4])];

