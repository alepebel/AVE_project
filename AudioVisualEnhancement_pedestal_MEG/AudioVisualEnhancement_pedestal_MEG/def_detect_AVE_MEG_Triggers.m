function triggers = def_detect_AVE_MEG_Triggers()

triggers = [];

triggers.detect.block = 1;                 % Range 1:6, for first, second... blocks

triggers.detect.resp = 8;                 % Range 9:10, 8 no response; for first, second button

triggers.detect.stimOnset = 10;              % Range 10:17 for 0:22.5:157.5
triggers.loc.stimOffset = 20;             % Range 20, 21 for NonOddball, Oddball
triggers.loc.dimOnset = 30;               

triggers.ET.stimOnset = 50;              % Range 50:82, for all ET positions

triggers.main.preFixOnset = 90;
triggers.main.stim1Onset = 100;           % Range 100:103, for 22.5, 67.5, 112.5, 157.5
triggers.main.stim1Offset = 110;          % Range 110:111, for Non-omission, Omission
triggers.main.stim2Onset = 120;           % Range 120:123, for 22.5, 67.5, 112.5, 157.5
triggers.main.stim2Offset = 130;          % Range 130:131, for Expected, Unexpected
triggers.main.stim3Onset = 140;           % Range 140:143, for 22.5, 67.5, 112.5, 157.5
triggers.main.stim3Offset = 150;          % Range 150:151, for contingeny: repetition expected, alternation expected
triggers.main.resp = 160;                 % Range 161:162, for first, second button
triggers.main.impulseOnset = 170;
triggers.main.blankOnset = 180;

triggers.id.preFixOnset = 200;              % Range 200:201, for oriSet
triggers.id.stim1Onset = 205;               % Range 205:208, for 22.5, 67.5, 112.5, 157.5
triggers.id.stim1Offset = 210;              % Range 210:211, for contingency         
triggers.id.stim2Onset = 215;               % Range 115:118, for 22.5, 67.5, 112.5, 157.5
triggers.id.stim2Offset = 220;              % Range 220:221, for Expected, Unexpected
triggers.id.resp = 225;                     % Range 226:227, for first, second button
triggers.id.blankOnset = 230;               % Range 230:232 (260: no response; 261: incorrect; 262: correct)

triggers.control.Interrupt = 254;
triggers.control.Resume = 255;

end

