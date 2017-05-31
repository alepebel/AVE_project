try

    
    Screen('FillRect', win, gray);
    Screen('TextSize',  win, 25);
    linesep = 2;
    max_n_ofcharperline = 90;

    
    title= ['VISUAL DISCRIMINATION TASK'];
    DrawFormattedText(win,title , 'center', hMid-260, [0,200,0], max_n_ofcharperline,0,0,linesep);
     
    Screen('TextSize',  win, 18);
    firstext = ['In this experiment you have to discriminate a noisy vertical grating from pure noise. ', ...
        'At the beginning of each trial, you will see a black small circle where you have to fixate your gaze during the experiment. ',...
        'Then, for a very brief period of time you will see a small or a big noise annulus patch. In some trials it will flicker containing a vertical grating blurred with visual noise. ',... 
        'After the stimulus presentation, is your time to report what you perceived in the annulus visual patch:',...
        'Use "Y" if you detected a grating and "U" if you only saw noise. ',...
        'At the end of each trial, the fixation circle will turn blue, and you can rest your eyes. You will receive feedback about your performance ',...
        '(green = correct, red = incorrect, white = not response). Do not forget to keep your eyes at the black fixation point during the trials'];
    DrawFormattedText(win,firstext , 'center', hMid-200, white, max_n_ofcharperline,0,0,linesep);
    
      
    Screen('TextSize',  win, 18);
    soundstext = ['During the experiment in some trials you will hear sounds or the fixation point will brieftly change to white. ',...
        'These sounds do not predict the presence of the grating signal and are not relevant for the task.'];
    DrawFormattedText(win, soundstext , 'center', hMid+180, white, max_n_ofcharperline,0,0,linesep);
    
    
    Screen('TextSize',  win, 18);
    secondtext = ['REMEMBER: Use "Y" if you detected a grating and "U" if you did not.' ];
    DrawFormattedText(win,secondtext , 'center', hMid+290, [255,0,0], max_n_ofcharperline,0,0,linesep);
    
    thirdtext = ['Press the spacebar to continue.' ];
    Screen('TextSize',  win, 16);
    DrawFormattedText(win,thirdtext , 'center', hMid+330, black, max_n_ofcharperline,0,0,linesep);
    
    %printﬂ
    Screen(win,'Flip');
end