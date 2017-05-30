try
    Screen('FillRect', win, gray);
    Screen('TextSize',  win, 28);
    linesep = 2;
    max_n_ofcharperline = 90;
    
    
    title= ['BLOCK FINISHED: PROPORTION OF CORRECT RESPONSES = ', num2str(100*correct_responses),' %'];
    DrawFormattedText(win,title , 'center', 'center', [0,0,255], max_n_ofcharperline,0,0,linesep);
     correct_responses
 
    %firstext = ['your detection threshold is'];
    %DrawFormattedText(win,firstext , 'center', hMid-70, white, max_n_ofcharperline,0,0,linesep);
        
    thirdtext = ['Press the spacebar to continue.' ];
    Screen('TextSize',  win, 18);
    DrawFormattedText(win,thirdtext , 'center', hMid+120, black, max_n_ofcharperline,0,0,linesep);
    
    %printﬂ
    Screen(win,'Flip');
end