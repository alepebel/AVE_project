try

    
    Screen('FillRect', win, gray);
    Screen('TextSize',  win, 25);
    linesep = 2;
    max_n_ofcharperline = 90;

    
    title= ['VISUAL DISCRIMINATION TASK'];
    DrawFormattedText(win,title , 'center', hMid-300, [0,200,0], max_n_ofcharperline,0,0,linesep);
     
    Screen('TextSize',  win, 18);
    firstext = ['In this experiment you have to discriminate a noisy vertical grating from pure noise. ', ...
        'At the beginning of each trial, you will see a black small circle where you have to fixate your gaze during the experiment. ',...
        'Then, for a very brief period of time you will see a big circular noise patch containing either a dim grating blurred with visual noise or only visual noise. ',... 
        'After the stimulus presentation, you have to report what you perceived in the visual patch:',...
        'Use "Y" if you detected a visual perturbation and "U" if you only saw static noise. ',...
        'At the end of each trial you will see a blue circle with a circular mask (you can rest your eyes here) and you will receive feedback about your performance ',...
        '(green = correct, red = incorrect, white = not response). Do not forget to keep your eyes at the black fixation point during the trials'];
    
    DrawFormattedText(win,firstext , 'center', hMid-240, white, max_n_ofcharperline,0,0,linesep);
    
    Screen('TextSize',  win, 18);
    secondtext = ['REMEMBER: Use "Y" if you detected a grating and "U" if you did not.' ];
    DrawFormattedText(win,secondtext , 'center', hMid+140, [255,0,0], max_n_ofcharperline,0,0,linesep);
    
    thirdtext = ['Press the spacebar to continue.' ];
    Screen('TextSize',  win, 16);
    DrawFormattedText(win,thirdtext , 'center', hMid+180, black, max_n_ofcharperline,0,0,linesep);
    
    %printﬂ
    Screen(win,'Flip');
end