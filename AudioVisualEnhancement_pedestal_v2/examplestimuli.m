% Drawing 2 example stimuli> signal and noise and displaying both
% sequentially

% here draw mask texture
 tilt = pi/180*detect.orientations(1); % signal example
cfg.gaborang = tilt; % Gabor orientation (rad)
cfg.gaborphi = rand; % pi/180*90; % Gabor unit phase CHECK THIS (radians?)
cfg.gaborcon = 0.2; % Gabor Michelson contrast

% clip patch luminance to [black,white] range to avoid black clusters at
% high contrasts

while true
    img = make_gabor_noise(cfg);
    % center example
    indexes_sf = find( mask(1).filter == 0);
    img(indexes_sf) = 0.5;
    %[patchimg] = absoluteDifferenceBetweenWhiteAndGray * patchimg;
    img = (gray * img)/0.5;
    if all(img(:) > black & img(:) < white)
        break
    end
end



noisexampletex = Screen('MakeTexture', win,  img); %Texture generated


cfg.gaborphi = rand; % pi/180*90; % Gabor unit phase CHECK THIS (radians?)
cfg.gaborcon = 0.2; % Gabor Michelson contrast

% clip patch luminance to [black,white] range to avoid black clusters at
% high contrasts

while true
    img = make_gabor_noise(cfg);
    % surround example
    indexes_sf = find( mask(2).filter == 0);
    img(indexes_sf ) = 0.5;
    img = (gray * img)/0.5;
    if all(img(:) > black & img(:) < white)
        break
    end
end



signalexampletex = Screen('MakeTexture', win,  img); %Texture generated


examplesrect = Screen('Rect', signalexampletex); %Extract information about texture size
scale = 1.0;
posx = [w/2];   %coordinates
posy = [h/2];
poscenternoise = [w*1/3, h/2];
poscentersignal = [w*2/3, h/2];
rotAngles = [0]; % we create the gratings already rotated. This can be applied to textures
noiseRects = CenterRectOnPoint(examplesrect*scale , 0.75*w/3, posy);
signalRects = CenterRectOnPoint(examplesrect*scale , 2.25*w/3, posy);

Screen('DrawTextures', win,noisexampletex,[],noiseRects,rotAngles); %, 2 is the adaptor surround
Screen('DrawTextures', win, signalexampletex,[],signalRects,rotAngles); %, 2 is the adaptor surround

Screen('TextSize',  win, 25);
linesep = 2;
max_n_ofcharperline = 30;

title= ['Small ring'];
DrawFormattedText(win,title , 0.65*w/3,'center', [0,0,0], max_n_ofcharperline,0,0,linesep);

title= ['Big ring'];
DrawFormattedText(win,title , 2.15*w/3,'center', [0,0,0], max_n_ofcharperline,0,0,linesep);

    thirdtext = ['These are noise and signal stimuli example. Press the spacebar to continue.' ];
    Screen('TextSize',  win, 16);
    DrawFormattedText(win,thirdtext , 'center', hMid+260, black, max_n_ofcharperline,0,0,linesep);
    
Screen('DrawingFinished', win);
vbl = Screen('Flip', win);


