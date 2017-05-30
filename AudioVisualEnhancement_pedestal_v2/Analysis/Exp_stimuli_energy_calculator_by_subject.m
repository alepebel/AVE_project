%run('M:\Experiments\pathdef\pathdef.m')
%open('M:\Experiments\pathdef\pathdef.m')
path = mfilename('fullpath'); % Finds current location of file
[pathstr] = fileparts(path); % Removes file information, saves directory
path = pathstr;

modalities = {'V' 'AV'};

% loading data
%main_folder = '/Volumes/ALEX_EXT/Experimentos/DondersExperiments/AudioVisualEnhancement/'
%main_folder = 'M:\Experiments\AudioVisualEnhancement\'
main_folder = 'P:\3018029.02\AudioVisualEnhancement\'
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
subjects = {'001','002','003','004','005','006','007','008','009','010','011'}
gaborang = pi/180*(1:1:180); % filter orientations (rad)
%gaborang = pi/180*(60:1:120); % filter orientations (rad)

for isubj = 1 : length(subjects)
    
    load([main_folder, 'results' filesep ,char(subjects(isubj)),filesep,char(subjects(isubj)),'_AVE_detect.mat'])
    
    n_sessions = size(Experimentalsession,2);
    for imask = 1 : length(mask)
        
        for isession = 1 : n_sessions
            disp(['session number ' num2str(isession)])
            load([main_folder, 'results' filesep ,char(subjects(isubj)), filesep, 'imgfiles',filesep,'IMG_run_',num2str(isession),'.mat'])
            disp(['processing imgs -> ', char(subjects(isubj)), filesep, 'imgfiles',filesep,'IMG_run_',num2str(isession),'.mat']);
            n_stimuli = length(stimuli);
            
            
            if isession == 1 % only precompute weights for the first session
                screenwidhsize = Experimentalsession(isession).detect.screenwidhsize; % in cm Mac 13 = 33
                screen_width_resolution = Experimentalsession(isession).detect.screen_width_resolution; % in Mac 13 2560 pixels
                distance2screen = Experimentalsession(isession).detect.distance2screen; % in cm
                
                pixel_size = screenwidhsize/screen_width_resolution; %Pixels are squared, so we only need to know one dimension
                one_degree_length_incm = tan(deg2rad(1/2))*distance2screen*2;
                one_degree_length_in_px = one_degree_length_incm/pixel_size;
                %one_degree_length_in_px = 2*distance2screen*tan((1/2)*(pi/180))*(screen_width_resolution/screenwidhsize);
                
                % add convolve2 function to path
                % addpath('./convolve2');
                % set basic noisy Gabor parameters (define in each trial contrast and phase
                % during the loop(
                
                ppd = round(one_degree_length_in_px); % pixels per degree of visual angle
                
                
                % set noisy Gabor parameters
                cfg = Experimentalsession(isession).detect.cfg;
                
                
                % creating different masks
                patchsiz = cfg.patchsiz;
                
                [rr cc] = meshgrid(1:patchsiz,1:patchsiz);
                mask_radius = patchsiz/2;
                % no mask
                mask(1).filter = ones(patchsiz);
                % center mask
                mask(2).filter = sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)<=mask_radius/2;  %in pixels
                % periphery mask
                %mask(3).filter = (sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)<=mask_radius & sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)>=mask_radius/2);  %in pixels
                mask(3).filter = sqrt((rr-(patchsiz/2)).^2+(cc-(patchsiz/2)).^2)>=mask_radius/2;  %in pixels
                
                % precompute signal energy filters
                
                for j = 1:length(gaborang)
                    %  disp(j)
                    % precomputing grating weights
                    cfg_e{j} = get_patch_energy(setfield(cfg,'gaborang',gaborang(j)));
                    % applying mask
                    %  masked_weight(imask).cfg_e{j} = cfg_e{j};
                    % intersecting mask with stimulus
                    % masked_weight(imask).cfg_e{j}.gabor{1,1} = mask(imask).filter.*cfg_e{j}.gabor{1,1};
                    %  masked_weight(imask).cfg_e{j}.gabor{1,2} = mask(imask).filter.*cfg_e{j}.gabor{1,2};
                    cfg_e{j}.gabor{1,1} = mask(imask).filter.*cfg_e{j}.gabor{1,1};
                    cfg_e{j}.gabor{1,2} = mask(imask).filter.*cfg_e{j}.gabor{1,2};
                end
                % imagesc(masked_weight(3).cfg_e{50}.gabor{1,1}) imagesc(masked_weight(2).cfg_e{50}.gabor{1,1})
                
                
            end
            % calculating energies for each stimulus, orientation and mask
            
            
            
            for istimuli = 1:n_stimuli
                for j = 1:length(gaborang)
                    energy.session(isession).mask(imask).stim(istimuli,j) = get_patch_energy(cfg_e{j},stimuli(istimuli).img);
                end
            end
        end
        
    end
    
    
    clear stimuli % delete the very heavy stimuli variable
    %energies_file_save = [main_folder,'Analysis' filesep 'computed_stimuli_energies.mat'];
    %save(energies_file_save,'energy');
    
    
    % preparing data and calculating energy combining all sessions
    %for isubj = 1:length(subjects)
    
    load([main_folder, 'results' filesep ,char(subjects(isubj)),filesep,char(subjects(isubj)),'_AVE_detect.mat'])
    
    n_sessions = size(Experimentalsession,2);
    
    for imask = 1 : length(mask)
        all_session_energy =    [];
        all_session_responses = [];
        
        for isession = 1 : n_sessions
            
            % calculating energies for each stimulus, orientation and mask
            
            aux = energy.session(isession).mask(imask).stim(:,:);
            all_session_energy =  [all_session_energy; aux];
            % responses_recoded = zeros(length(Experimentalsession(isession).detect.responses_vector),1 )
            % responses_recoded(find(strcmp(Experimentalsession(isession).detect.responses_vector(:,10) ,'y'))) = 1;
            aux_resp = Experimentalsession(isession).detect.responses_vector;
            all_session_responses = vertcat(all_session_responses, aux_resp);
            
        end
        % collapsing all trials together
        data_modelling.mask(imask).all_session_energy = all_session_energy;
        data_modelling.mask(imask).all_session_trials_conditions = all_session_responses;
        data_modelling.mask(imask).all_session_trials_condition_names = Experimentalsession(isession).detect.responses_vector_variables;
        responses_recoded = zeros(size( all_session_responses ,1),1 );
        responses_recoded(find(strcmp( all_session_responses(:,10) ,'y'))) = 1;
        
        data_modelling.mask(imask).responses_recoded = responses_recoded;
        
        modality_indexes = [];
        
        %         for modal = 1 : length(modalities)
        %             % selecting V or AV stimuli for the correlation
        %             modality_indexes(modal,:) = find(strcmp(data_modelling.mask(imask).all_session_trials_responses(:,7) ,modalities(modal)));
        %             % calculating energy correlation trialwise
        %             for j = 1 : length(gaborang)
        %                 resp = data_modelling.mask(imask).responses_recoded(modality_indexes(modal,:))
        %                 energy_vector = data_modelling.mask(imask).all_session_energy(modality_indexes(modal,:),j);
        %                 data_modelling.subject(isubj).energy_sensitivity(imask).modality(modal).e(j,:) = glmfit(energy_vector,resp,'binomial','link','probit');
        %             end
        %
        %         end
        
    end
    
    
    % we save a mat file for each subject containing the experimental conditions matrix, the energy for each trial and orientation and the names of conditions
    % save_data_matrix en la carpeta analysis, crear un archivo para cada sujeto
    save([outputfolder_mat 'denergy_subject_' char(subjects(isubj))],'data_modelling');
    
    
    %% prepare data to export them to R
    data = [];
    %for isubj = 1:length(subjects)
    for imask = 1 : length(mask)
        
        for j = 1 : length(gaborang)
            
            energy_orient = data_modelling.mask(imask).all_session_energy(:,j);
            current_orient = (gaborang(j)*ones(length(data_modelling.mask(imask).all_session_energy(:,j)),1))/pi*180;
            current_mask = (imask *ones(length(data_modelling.mask(imask).all_session_energy(:,j)),1));
            
            experiment_vector = data_modelling.mask(imask).all_session_trials_conditions;
            
            recoded_resp = data_modelling.mask(imask).responses_recoded;
            
            % putting together conditions
            aux = [experiment_vector,   num2cell(current_mask), num2cell(current_orient), num2cell(energy_orient),num2cell( recoded_resp)];
            data = [data; aux];
        end
        
    end
    
    
    
    % exporting data to calculate regression in R better
    
    % save data in a txt file
    dataoutput = [outputfolder_R 'denergy_subject_' char(subjects(isubj)),'.txt']
    
    %['C:\Users\alebel\Dropbox\DATA\AVE_experiment\Energy_correlations\first_analyses\detect_datos_energy.txt'];
    
    % printint test data
    fileID = fopen([dataoutput],'w'); % try to write on top
    formatSpec = '%i %i %i %i %f %s %s %i %i %s %i %f %s %i %i %f %i\r\n';
    
    % for windows
    %dataadaptoutput = [path,'\results\',initials,'_detect.txt'];
    [nrows,ncols] = size(data);
    for row = 1 : nrows
        fprintf(fileID, formatSpec, data{row,:}) ;
    end
    
    fclose(fileID);
    
    
    
end






Experimentalsession(isession).detect.responses_vector
imagesc(data_modelling.subject(isubj).energy_sensitivity(imask).e(:,2))
plot((1:1:180),energy_sensitivity(1).e(:,2))
plot((gaborang)/pi*180,data_modelling.subject(1).energy_sensitivity(imask).modality(1).e(:,2));
        
plot((gaborang)/pi*180,data_modelling.subject(2).energy_sensitivity(imask).e(:,2))
plot((gaborang)/pi*180,data_modelling.subject(1).energy_sensitivity(imask).e(:,2))
imagesc(energy_sensitivity(1).e(:,1))
imagesc(energy_sensitivity(1).e(:,2))

for j = 1:length(gaborang)
    b = glmfit([e(:,j)],resp,'binomial','link','probit');
    energy.session(isession).mask(imask).stim(istimuli,j)
    eprime_hat(i,:) = b(2:3)*sqrt(2);
    xprime_hat(i,1) = b(4);
end

end

imagesc(energy.session(isession).mask(imask).stim)


% save results
save('test_gabor_noise.mat','n','gaborang','x','e','cfg');

energy.mask(imask).stim(istimuli,j)


