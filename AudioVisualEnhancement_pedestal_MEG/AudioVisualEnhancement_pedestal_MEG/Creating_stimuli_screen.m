try
    Screen('FillRect', win, gray);
    Screen('TextSize',  win, 14);
    linesep = 2;
    max_n_ofcharperline = 60;
        
    title= ['Creating the stimuli. Please Wait. Percentage completed: ',num2str(percetange_stimuli),' %'];
    DrawFormattedText(win,title , 'center', 'center', [0,0,0], max_n_ofcharperline,0,0,linesep);
    
    %firstext = ['your detection threshold is'];
    %DrawFormattedText(win,firstext , 'center', hMid-70, white, max_n_ofcharperline,0,0,linesep);
        
    % Screen('TextSize',  win, 10);
    % DrawFormattedText(win,thirdtext , 'center', hMid+120, gray, max_n_ofcharperline,0,0,linesep);
    
    %printﬂ
    Screen(win,'Flip');
end