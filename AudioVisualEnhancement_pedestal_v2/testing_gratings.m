ppd = 30; % pixels per degree of visual angle

% set noisy Gabor parameters
cfg          = [];
cfg.patchsiz = ppd*10.0; % patch size (pix)
cfg.patchenv = ppd*0.75; % patch spatial envelope s.d. (pix)
cfg.envscale = 10^3; % patch spatial envelope s.d. (pix)
cfg.patchlum = 0.5; % patch background luminance
cfg.gaborper = ppd*1.1; % Gabor spatial period (pix)
cfg.noisedim = cfg.gaborper/6; % noise dimension (pix)
cfg.noisecon = 0.3/3; % noise RMS contrast
cfg.diameter = cfg.patchsiz; % This way scale the gaussian evelop and I make the border very abrupt
cfg.gaborang = tilt; % Gabor orientation (rad)
                cfg.gaborphi = rand; % pi/180*90; % Gabor unit phase CHECK THIS (radians?)
                cfg.gaborcon = 0.4; % Gabor Michelson contrast
                
                while true
                    img = make_gabor_noise(cfg);
                    %[patchimg] = absoluteDifferenceBetweenWhiteAndGray * patchimg;
                    img = (gray * img)/0.5;
                    if all(img(:) > black & img(:) < white)
                        break
                    end
                end
                
                
                % save all the information about the visual stimulus
                %   quest.Grating(iTrial).texture = gray + absoluteDifferenceBetweenWhiteAndGray * imageMatrix;
                %   quest.Grating(iTrial).tilt = tilt;
                %   quest.Grating(iTrial).contrast = contrast;
                %   quest.Grating(iTrial).freq = freq;
                %   quest.Grating(iTrial).phase = phases(1); % there is only one grating presented.
                
                gabortex=Screen('MakeTexture', win,  img); %Texture generated
                % imagesc(n)
                texrect = Screen('Rect', gabortex); %Extract information about texture size
                scale = 1; % scaling artifitially the stimulus; Use carefully. Default 1
                posx = [w/2];   %coordinates
                posy = [h/2];
                poscenter = [w/2, h/2];
                rotAngles = [0]; % we create the gratings already rotated. This can be applied to textures
                
                dstRects = CenterRectOnPoint(texrect*scale , posx, posy);

color_fp = [255 0 0];
color_fp = [0 0 255]
Screen('DrawTextures', win, gabortex,[],dstRects,rotAngles); %, 2 is the adaptor surround
              %  Screen('FillOval', win, color_fp, centeredRect_fp, radious_fp);
                Screen('DrawingFinished', win);
                Vflag = 1; % Avoiding enter here and redraw the stimulus
                vbl = Screen('Flip', win, vbl + 0.5 * ifi);
                
                