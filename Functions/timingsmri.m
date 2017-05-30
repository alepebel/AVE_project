
sompath =strcat('/Users/jaynecobb/Dropbox/BaylorExperiments/'); 
run = 2
datos = oPar.sched(run).FullSchedInfo.curr_sched;

conds = unique(datos(:,2));

for idx = 1 : length(conds)
    
    somfn       = strcat('condt',num2str(idx),'.txt');
    outfilename = fullfile(sompath,somfn);
    fid         = fopen(outfilename, 'w');
    
    ff = find(datos(:,2) == idx);
    
    t = datos(ff);
 
    fprintf(fid,'%f\t',t);    
    
    fclose(fid);
    
end