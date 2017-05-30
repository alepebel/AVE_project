%run('M:\Experiments\pathdef\pathdef.m')
%open('M:\Experiments\pathdef\pathdef.m')
path = mfilename('fullpath'); % Finds current location of file
[pathstr] = fileparts(path); % Removes file information, saves directory
path = pathstr;

modalities = {'V' 'AV'};

% loading data
%main_folder = '/Volumes/ALEX_EXT/Experimentos/DondersExperiments/AudioVisualEnhancement/'
%main_folder = 'M:\Experiments\AudioVisualEnhancement\'
main_folder = 'P:\3018029.02\AudioVisualEnhancement_pedestal_v2\'
cd(main_folder)

% here we will save all the output of energy calculations in a .mat file.
outputfolder_mat = [main_folder 'Analysis' filesep 'matlab_format_data' filesep]
% here we will save all the output of energy calculations in . txt file for R.
outputfolder_R = [main_folder  'Analysis' filesep 'R_format_data' filesep]

if ~exist(outputfolder_mat,'dir')
    mkdir(outputfolder_mat);
end
if ~exist(outputfolder_R,'dir')
    mkdir(outputfolder_R);
end

%subjects = {'001'}
subjects = {'024','025','026','027','028','029','030','031','032','033','034','035','036','037'}
subjects = {'038' '039' '040' '041' '042' '043' '044' '045' '046' '047' '048' '049' '050' '051' '052' '053'}
gaborang = pi/180*(1:1:180); % filter orientations (rad)
%gaborang = pi/180*(60:1:120); % filter orientations (rad)

for isubj = 1 : length(subjects)
    
    load([main_folder, 'results' filesep ,char(subjects(isubj)),filesep,char(subjects(isubj)),'_AVE_detect.mat'])
    
    n_sessions = size(Experimentalsession,2);
    
        all_session_energy =    [];
        all_session_responses = [];
        
        sessions = 1 : n_sessions;
        
        if strcmp(subjects(isubj),'024')
        sessions = [1 2 4]; % in subject 24 session 3 is missing
   
        end
        
        for isession = sessions 
      
            % calculating energies for each stimulus, orientation and mask          
            aux = Experimentalsession(1,isession).detect.energies; % energy.session(isession).mask(imask).stim(:,:);
            all_session_energy =  [all_session_energy; aux];
            % responses_recoded = zeros(length(Experimentalsession(isession).detect.responses_vector),1 )
            % responses_recoded(find(strcmp(Experimentalsession(isession).detect.responses_vector(:,10) ,'y'))) = 1;
            aux_resp = Experimentalsession(isession).detect.responses_vector;
            all_session_responses = vertcat(all_session_responses, aux_resp);
            
        end

        % In case of multiple responses, register just the second response
        for itrial = 1:size( all_session_responses,1)
            if length(all_session_responses{itrial,12}) == 2;
                  all_session_responses{itrial,12} = all_session_responses{itrial,12}{1,2};
            end
            
        end
        
       
        % collapsing all trials together
        data_modelling.all_session_energy = all_session_energy;
        data_modelling.all_session_trials_conditions = all_session_responses;
        data_modelling.all_session_trials_condition_names = Experimentalsession(isession).detect.responses_vector_variables;
        responses_recoded = zeros(size( all_session_responses ,1),1 );
        responses_recoded(find(strcmp( all_session_responses(:,12) ,'y'))) = 1;
        
        data_modelling.responses_recoded = responses_recoded;
        
        modality_indexes = [];
        
    
    % we save a mat file for each subject containing the experimental conditions matrix, the energy for each trial and orientation and the names of conditions
    % save_data_matrix en la carpeta analysis, crear un archivo para cada sujeto
    save([outputfolder_mat 'denergy_subject_' char(subjects(isubj))],'data_modelling');
    
    
    %% prepare data to export them to R
    data = [];
    %for isubj = 1:length(subjects)

        for j = 1 : length(gaborang)
            
            energy_orient = data_modelling.all_session_energy(:,j);
            current_orient = (gaborang(j)*ones(length(data_modelling.all_session_energy(:,j)),1))/pi*180;
            % current_mask = (imask *ones(length(data_modelling.mask(imask).all_session_energy(:,j)),1));
            experiment_vector = data_modelling.all_session_trials_conditions;
            recoded_resp = data_modelling.responses_recoded;
            
            % putting together conditions
            aux = [experiment_vector,   num2cell(current_orient), num2cell(energy_orient),num2cell( recoded_resp)];
            data = [data; aux];
        end

    
    % exporting data to calculate regression in R better
    
    % save data in a txt file
    dataoutput = [outputfolder_R 'denergy_subject_' char(subjects(isubj)),'.txt']
    
    %['C:\Users\alebel\Dropbox\DATA\AVE_experiment\Energy_correlations\first_analyses\detect_datos_energy.txt'];
    
    % printint test data
    fileID = fopen([dataoutput],'w'); % try to write on top
    formatSpec = '%i %i %f %f %f %s %s %s %i %f %f %s %i %f %s %i %f %i\r\n';
    
    % for windows
    %dataadaptoutput = [path,'\results\',initials,'_detect.txt'];
    [nrows,ncols] = size(data);
    for row = [1:nrows ] % 684, 687: nrows]
        fprintf(fileID, formatSpec, data{row,:}) ;
        row;
    end
    
    fclose(fileID);
    
   % size(data{1:row,12})
    %size(data{row,12},2)
    
end





