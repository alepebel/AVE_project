try
    
    Screen('FillRect', win, gray);
    Screen('TextSize',  win, 25);
    linesep = 2;
    max_n_ofcharperline = 90;
    
    if strcmp(procedure,'oddball')
        title= ['COUNTING FLICKER EVENTS TASK'];
        DrawFormattedText(win,title , 'center', hMid-100, [200,0,0], max_n_ofcharperline,0,0,linesep);
        
        Screen('TextSize',  win, 18);
        firstext = ['In this part of the experiment you have to COUNT how many times the black rings flickered in each trial. ', ...
            ''];
        DrawFormattedText(win,firstext , 'center', hMid-50, white, max_n_ofcharperline,0,0,linesep);
        
        Screen('TextSize',  win, 18);
        secondtext = ['REMEMBER: Press the button located in the same side as the number of times "0, 1 or 2" that you saw the ring flicker.' ];
        DrawFormattedText(win,secondtext , 'center', hMid+15, [255,255,255], max_n_ofcharperline,0,0,linesep);
        
        thirdtext = ['If you are ready, press one of the buttons to continue..' ];
        Screen('TextSize',  win, 16);
        DrawFormattedText(win,thirdtext , 'center', hMid+100, black, max_n_ofcharperline,0,0,linesep);
        
    else
        
        title= ['GRATING DETECTION TASK'];
        DrawFormattedText(win,title , 'center', hMid-100, [0,20,200], max_n_ofcharperline,0,0,linesep);
        
        Screen('TextSize',  win, 18);
        firstext = ['In this part of the experiment you have to DETECT a noisy VERTICAL GRATING from pure noise. ', ...
            ''];
        DrawFormattedText(win,firstext , 'center', hMid-50, white, max_n_ofcharperline,0,0,linesep);
        
        Screen('TextSize',  win, 18);
        secondtext = ['REMEMBER: Press the button located in the same side as "Y" if you detected a grating and the button in the same side as "N" if you did not.' ];
        DrawFormattedText(win,secondtext , 'center', hMid+15, [255,255,255], max_n_ofcharperline,0,0,linesep);
        
        thirdtext = ['If you are ready, press one of the buttons to continue.' ];
        Screen('TextSize',  win, 16);
        DrawFormattedText(win,thirdtext , 'center', hMid+100, black, max_n_ofcharperline,0,0,linesep);
        
        
    end

    Screen(win,'Flip');
end