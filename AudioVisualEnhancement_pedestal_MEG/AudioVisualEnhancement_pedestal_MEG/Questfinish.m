try
    Screen('FillRect', win, gray);
    Screen('TextSize',  win, 30);
    linesep = 2;
    max_n_ofcharperline = 90;
    
    
    title= ['QUEST FINISHED'];
    DrawFormattedText(win,title , 'center', 'center', [0,0,255], max_n_ofcharperline,0,0,linesep);
    
    %firstext = ['your detection threshold is'];
    %DrawFormattedText(win,firstext , 'center', hMid-70, white, max_n_ofcharperline,0,0,linesep);
        
    thirdtext = ['Press the spacebar to continue.' ];
    Screen('TextSize',  win, 18);
    DrawFormattedText(win,thirdtext , 'center', hMid+120, black, max_n_ofcharperline,0,0,linesep);
    
    %print�
    Screen(win,'Flip');
end