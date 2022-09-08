%%% Compute multiple functional profiles (alff/falff, Timescale, )

clear; clc;
path = '/Data_TimeCourses/HCP_S1200_rfMRI_REST_FIX'; % the timecourse directory
outpath = '/Data_HCP_S1200/REST1'; % the output metrics directory
sbj = importdata('/List/HCP_S1200_REST_4-session_excluBad25_sbjlist.txt'); % the subjects list
path_wb_command = '/usr/local/workbench-linux64-v1.4.2/bin_linux64/wb_command'; % workbench path

Surf_L = '/Templates/S1200.L.midthickness_MSMAll.32k_fs_LR.surf.gii';
Surf_R = '/Templates/S1200.R.midthickness_MSMAll.32k_fs_LR.surf.gii';
Mask_L = '/Templates/roi_fsLR_32k_L.shape.gii';
Mask_R = '/Templates/roi_fsLR_32k_R.shape.gii';

%% Parameters Setup
tr = 0.72;
alff_band = [0.01,0.08]; % frequency band
ReHo_neighbour = 2;
NetSparsity = 90:-2:70; % nework sparsity range
NetSpar_fname = '/List/Sparsity.txt';
dlmwrite(NetSpar_fname,NetSparsity,'delimiter','\n');
% Modulaity_Method = '1'; % '1' for a modified greedy optimization algorithm,  '2' for a spectral optimization algorithm;
Num_Rand = 1000;
% Load parcel mask
File_parcel = '/Templates/Glasser_360.dlabel.nii';
ParcelHeader = ciftiopen(File_parcel);
Parcel = ParcelHeader.cdata;
% Pscalar template load
Pscalar = ciftiopen('/Templates/sub-MSC01_MyelinMap_biascorrected_Glasser360.pscalar.nii');
    
%% Main Body

for sub = 1 : length(sbj)
    t1 = clock;
    disp(['...................', num2str(sbj(sub)),' Starting ...................'])
    disp(datestr(now)) 
    
    % TC file dir
    File_TC = [path, filesep, num2str(sbj(sub)), filesep, 'rfMRI_REST1_LRRL_Atlas_MSMAll_hp2000_clean_Neocortex.dtseries.nii'];    
    
    %% ALFF & fALFF Calculation
    disp(['----------  Computing ALFF & fALFF for ', num2str(sbj(sub)), '  ----------'])
    ALFFoutPath = [outpath, filesep, num2str(sbj(sub)), filesep, 'ALFF'];
    mkdir(ALFFoutPath)
    ALFFoutFile = cell(2,1);
    
    % Alff & fALFF Calculation
    ALFFoutFile{1} = [ALFFoutPath, filesep, num2str(sbj(sub)),'_NeoCortex_ALFF_32k_fsLR.dscalar.nii'];
    ALFFoutFile{2} = [ALFFoutPath, filesep, num2str(sbj(sub)),'_NeoCortex_fALFF_32k_fsLR.dscalar.nii'];
    FSN_alff_falff(File_TC, tr, alff_band(2), alff_band(1), ALFFoutFile); 
   
    % Parcel-wise ALFF & fALFF
    ALFF = ciftiopen(ALFFoutFile{1});
    ALFF = ALFF.cdata;
    fALFF = ciftiopen(ALFFoutFile{2});
    fALFF = fALFF.cdata;
    ALFF_pclwise = zeros(length(Parcel(any(unique(Parcel),2))),1);
    fALFF_pclwise = zeros(length(Parcel(any(unique(Parcel),2))),1);   
    for pcl = 1 : length(Parcel(any(unique(Parcel),2)))
        Index_pcl = find(Parcel==pcl);
        ALFF_pcl = mean(ALFF(Index_pcl,:));
        ALFF_pclwise(pcl) = ALFF_pcl;
        fALFF_pcl = mean(fALFF(Index_pcl,:));
        fALFF_pclwise(pcl) = fALFF_pcl;
    end    
    ALFF_pcl_cii = Pscalar;
    fALFF_pcl_cii = Pscalar;
    ALFF_pcl_cii.cdata = ALFF_pclwise;
    fALFF_pcl_cii.cdata = fALFF_pclwise;
    ALFF_pcl_fname = [ALFFoutPath, filesep, num2str(sbj(sub)),'_NeoCortex_ALFF_32k_fsLR_Glasser_360.pscalar.nii'];
    fALFF_pcl_fname = [ALFFoutPath, filesep, num2str(sbj(sub)),'_NeoCortex_fALFF_32k_fsLR_Glasser_360.pscalar.nii'];
    cifti_write(ALFF_pcl_cii, ALFF_pcl_fname)
    cifti_write(fALFF_pcl_cii, fALFF_pcl_fname)
    
    clear ALFFoutPath ALFFoutFile ALFF_pcl_fname fALFF_pcl_fname ALFF_pcl_cii fALFF_pcl_cii
    disp(['----------  ALFF & fALFF for ', num2str(sbj(sub)), ' Finished  ----------', newline, newline])
    
    %% Reho Calculation
    
    TCheader = cifti_read(File_TC);
    TC = TCheader.cdata;
    
    disp(['----------  Computing Reho for ', num2str(sbj(sub)), '  ----------'])
    RehooutPath = [outpath, filesep, num2str(sbj(sub)), filesep, 'Reho'];
    mkdir(RehooutPath)
    % Prerequisite Input File dirs
    RehoIn_L = zeros(32492, size(TC, 2));
    TC_L = TC(TCheader.diminfo{1, 1}.models{1, 1}.start : TCheader.diminfo{1, 1}.models{1, 1}.start + TCheader.diminfo{1, 1}.models{1, 1}.count - 1, :);
    RehoIn_L(TCheader.diminfo{1, 1}.models{1, 1}.vertlist  + 1, :) = TC_L;
    
    RehoIn_R = zeros(32492, size(TC, 2));
    TC_R = TC(TCheader.diminfo{1, 1}.models{1, 2}.start : TCheader.diminfo{1, 1}.models{1, 2}.start + TCheader.diminfo{1, 1}.models{1, 2}.count - 1, :);
    RehoIn_R(TCheader.diminfo{1, 1}.models{1, 1}.vertlist  + 1, :) = TC_L;
    
    RehooutFile_L = [RehooutPath,filesep, num2str(sbj(sub)),'_Reho_32k_fsLR_L.func.gii'];
    RehooutFile_R = [RehooutPath,filesep, num2str(sbj(sub)),'_Reho_32k_fsLR_R.func.gii'];

    % ReHo Calcaulation (L & R Hemi Seperately)
    FSN_reho_Surf_YM(RehoIn_L, ReHo_neighbour, Mask_L, RehooutFile_L, Surf_L);
    FSN_reho_Surf_YM(RehoIn_R, ReHo_neighbour, Mask_R, RehooutFile_R, Surf_R);
    
    % Merge both Hemi to single 32k_fs_LR cii
    Command = [path_wb_command,' -cifti-create-dense-scalar ', RehooutPath,filesep, num2str(sbj(sub)),'_Reho_32k_fsLR.dscalar.nii -left-metric ', RehooutFile_L, ' -roi-left ', Mask_L, ' -right-metric ', RehooutFile_R, ' -roi-right ', Mask_R];
    system(Command)
    % Parcel-wise Reho
    Reho = ciftiopen([RehooutPath,filesep, num2str(sbj(sub)),'_Reho_32k_fsLR.dscalar.nii']);
    Reho = Reho.cdata;
	Reho_pclwise = zeros(length(Parcel(any(unique(Parcel),2))),1);    
    for pcl = 1 : length(Parcel(any(unique(Parcel),2)))
        Index_pcl = find(Parcel==pcl);
        Reho_pcl = mean(Reho(Index_pcl,:));
        Reho_pclwise(pcl) = Reho_pcl;
    end    
    Reho_pcl_cii = Pscalar;
    Reho_pcl_cii.cdata = Reho_pclwise;
    Reho_pcl_fname = [RehooutPath, filesep, num2str(sbj(sub)),'_NeoCortex_Reho_32k_fsLR_Glasser_360.pscalar.nii'];
    ciftisavereset(Reho_pcl_cii,Reho_pcl_fname)
  
    clear command RehoIn_L RehoIn_R RehooutFile_L RehooutFile_R  ReHooutPath Reho_pcl_cii Reho_pcl_fname
    disp(['----------  Reho for ', num2str(sbj(sub)), ' Finished  ----------', newline, newline])

%% Network Graph Theory Calculation
    
%% Parce-wise Network Construction
    % FC Construction (Parcel-wise)
    % Load TimeCourse(TC)

    disp(['----------  TimeCourse of ', num2str(sbj(sub)), ' Loaded  ----------', newline, newline])
    
    % Get Parcel TC
    TC_pclwise = zeros(size(TC,2),length(Parcel(any(unique(Parcel),2))));
    for pcl = 1 : length(Parcel(any(unique(Parcel),2)))
        Index_pcl = find(Parcel==pcl);
        TC_pcl = mean(TC(Index_pcl,:));
        TC_pclwise(:,pcl) = TC_pcl';
    end
    % Parcel-wise Functional Connectivity
    FC_pcl = corr(TC_pclwise);
%     FCMat_fisherz = atanh(FCMat);

    % Posive & Absloutr Network under multi sparsity
    FC_pcl_pos = FC_pcl.*(FC_pcl > 0);
    FC_pcl_abs = abs(FC_pcl);
    disp(['----------  Parcel-wise FC of ', num2str(sbj(sub)), ' Constructed  ----------', newline, newline])
    % postive network
    FC_pcl_pos_multiThreshlds = cell(length(NetSparsity),1);
    parfor spars = 1: length(NetSparsity)
        Threshld = prctile(tril2vec(FC_pcl_pos),NetSparsity(spars));
        FC_pcl_pos_threshlded = FC_pcl_pos.*(FC_pcl_pos > Threshld);
        FC_pcl_pos_multiThreshlds{spars} = FC_pcl_pos_threshlded;
    end
    % absloute network
    FC_pcl_abs_multiThreshlds = cell(length(NetSparsity),1);
    parfor spars = 1: length(NetSparsity)
        Threshld = prctile(tril2vec(FC_pcl_abs),NetSparsity(spars));
        FC_pcl_abs_threshlded = FC_pcl_abs.*(FC_pcl_abs > Threshld);
        FC_pcl_abs_multiThreshlds{spars} = FC_pcl_abs_threshlded;
    end
    disp(['----------  Parcel-wise FC of ', num2str(sbj(sub)), ' Thresholded (Absloute & Postive) ----------', newline, newline])

    %% Network Features Calcaulation

    % Out Dir Creatation
    NetFeature_Dir = [outpath, filesep, num2str(sbj(sub)), filesep, 'NetworkFeatures'];
    mkdir(NetFeature_Dir)
        
    %% Node Degree
    
    disp(['----------  Computing Node Degree for ', num2str(sbj(sub)), '  ----------'])
    ki_multiThreshlds_pos = zeros(length(Parcel(any(unique(Parcel),2))),length(NetSparsity));
    ki_multiThreshlds_abs = zeros(length(Parcel(any(unique(Parcel),2))),length(NetSparsity));    
    parfor spars = 1: length(NetSparsity)
        [~, ki_pos] = gretna_node_degree_weight(FC_pcl_pos_multiThreshlds{spars});
        [~, ki_abs] = gretna_node_degree_weight(FC_pcl_abs_multiThreshlds{spars});
        ki_multiThreshlds_pos(:, spars) = ki_pos';
        ki_multiThreshlds_abs(:, spars) = ki_abs';
    end
    ki_pos_cii = Pscalar;
    ki_pos_cii.cdata =  ki_multiThreshlds_pos;
    ki_abs_cii = Pscalar;
    ki_abs_cii.cdata =  ki_multiThreshlds_abs; 
    ki_pos_fname = [NetFeature_Dir, filesep, num2str(sbj(sub)), '_pos_NodeDegree_32k_fsLR_Glasser_360.pscalar.nii'];
    ki_abs_fname = [NetFeature_Dir, filesep, num2str(sbj(sub)), '_abs_NodeDegree_32k_fsLR_Glasser_360.pscalar.nii'];
    ciftisavereset(ki_pos_cii,ki_pos_fname)
    ciftisavereset(ki_abs_cii,ki_abs_fname)
    % Set _Glasser_360.pscalar.nii map name
    command = [path_wb_command, ' -set-map-names ', ki_pos_fname, ' -name-file ',NetSpar_fname];
    system(command)
    command = [path_wb_command, ' -set-map-names ', ki_abs_fname, ' -name-file ',NetSpar_fname];
    system(command)    
    clear ki_pos_cii ki_abs_cii ki_pos_fname ki_abs_fname command
    disp(['----------  Node Degree for ', num2str(sbj(sub)), ' Finished  ----------', newline, newline])
    
    %% Shortest Path Length
    
    disp(['----------  Computing Shortest Path Length for ', num2str(sbj(sub)), '  ----------'])
    Lpi_multiThreshlds_pos = zeros(length(Parcel(any(unique(Parcel),2))),length(NetSparsity));
    Lpi_multiThreshlds_abs = zeros(length(Parcel(any(unique(Parcel),2))),length(NetSparsity));   
    parfor spars = 1: length(NetSparsity)
        [~, Lpi_pos] = gretna_node_shortestpathlength_weight(FC_pcl_pos_multiThreshlds{spars});
        [~, Lpi_abs] = gretna_node_shortestpathlength_weight(FC_pcl_abs_multiThreshlds{spars});
        Lpi_multiThreshlds_pos(:, spars) = Lpi_pos';
        Lpi_multiThreshlds_abs(:, spars) = Lpi_abs';
    end
    Lpi_pos_cii = Pscalar;
    Lpi_pos_cii.cdata =  Lpi_multiThreshlds_pos;
    Lpi_abs_cii = Pscalar;
    Lpi_abs_cii.cdata =  Lpi_multiThreshlds_abs; 
    Lpi_pos_fname = [NetFeature_Dir, filesep, num2str(sbj(sub)), '_pos_ShortestPathLength_32k_fsLR_Glasser_360.pscalar.nii'];
    Lpi_abs_fname = [NetFeature_Dir, filesep, num2str(sbj(sub)), '_abs_ShortestPathLength_32k_fsLR_Glasser_360.pscalar.nii'];
    ciftisavereset(Lpi_pos_cii,Lpi_pos_fname)
    ciftisavereset(Lpi_abs_cii,Lpi_abs_fname)
    % Set _Glasser_360.pscalar.nii map name
    command = [path_wb_command, ' -set-map-names ', Lpi_pos_fname, ' -name-file ', NetSpar_fname];
    system(command)
    command = [path_wb_command, ' -set-map-names ', Lpi_abs_fname, ' -name-file ', NetSpar_fname];
    system(command)     
    clear Lpi_pos_cii Lpi_abs_cii Lpi_pos_fname Lpi_abs_fname command
    disp(['----------  Shortest Path Length for ', num2str(sbj(sub)), ' Finished  ----------', newline, newline])

    %% Local Efficiency
    
    disp(['----------  Computing Local Efficiency for ', num2str(sbj(sub)), '  ----------'])
    locEi_multiThreshlds_pos = zeros(length(Parcel(any(unique(Parcel),2))),length(NetSparsity));
    locEi_multiThreshlds_abs = zeros(length(Parcel(any(unique(Parcel),2))),length(NetSparsity));   
    parfor spars = 1: length(NetSparsity)
        [~, locEi_pos] = gretna_node_local_efficiency_weight(FC_pcl_pos_multiThreshlds{spars});
        [~, locEi_abs] = gretna_node_local_efficiency_weight(FC_pcl_abs_multiThreshlds{spars});
        locEi_multiThreshlds_pos(:, spars) = locEi_pos';
        locEi_multiThreshlds_abs(:, spars) = locEi_abs';
    end
    locEi_pos_cii = Pscalar;
    locEi_pos_cii.cdata =  locEi_multiThreshlds_pos;
    locEi_abs_cii = Pscalar;
    locEi_abs_cii.cdata =  locEi_multiThreshlds_abs; 
    locEi_pos_fname = [NetFeature_Dir, filesep, num2str(sbj(sub)), '_pos_local_efficiency_32k_fsLR_Glasser_360.pscalar.nii'];
    locEi_abs_fname = [NetFeature_Dir, filesep, num2str(sbj(sub)), '_abs_local_efficiency_32k_fsLR_Glasser_360.pscalar.nii'];
    ciftisavereset(locEi_pos_cii,locEi_pos_fname)
    ciftisavereset(locEi_abs_cii,locEi_abs_fname)
    % Set _Glasser_360.pscalar.nii map name
    command = [path_wb_command, ' -set-map-names ', locEi_pos_fname, ' -name-file ', NetSpar_fname];
    system(command)
    command = [path_wb_command, ' -set-map-names ', locEi_abs_fname, ' -name-file ', NetSpar_fname];
    system(command)     
    clear locEi_pos_cii locEi_abs_cii locEi_pos_fname locEi_abs_fname command
    disp(['----------  Local Efficiency for ', num2str(sbj(sub)), ' Finished  ----------', newline, newline])

    %% Global Efficiency
    
    disp(['----------  Computing Global Efficiency for ', num2str(sbj(sub)), '  ----------'])
    gEi_multiThreshlds_pos = zeros(length(Parcel(any(unique(Parcel),2))),length(NetSparsity));
    gEi_multiThreshlds_abs = zeros(length(Parcel(any(unique(Parcel),2))),length(NetSparsity));   
    parfor spars = 1: length(NetSparsity)
        [~, gEi_pos] = gretna_node_global_efficiency_weight(FC_pcl_pos_multiThreshlds{spars});
        [~, gEi_abs] = gretna_node_global_efficiency_weight(FC_pcl_abs_multiThreshlds{spars});
        gEi_multiThreshlds_pos(:, spars) = gEi_pos';
        gEi_multiThreshlds_abs(:, spars) = gEi_abs';
    end
    gEi_pos_cii = Pscalar;
    gEi_pos_cii.cdata =  gEi_multiThreshlds_pos;
    gEi_abs_cii = Pscalar;
    gEi_abs_cii.cdata =  gEi_multiThreshlds_abs; 
    gEi_pos_fname = [NetFeature_Dir, filesep, num2str(sbj(sub)), '_pos_global_efficiency_32k_fsLR_Glasser_360.pscalar.nii'];
    gEi_abs_fname = [NetFeature_Dir, filesep, num2str(sbj(sub)), '_abs_global_efficiency_32k_fsLR_Glasser_360.pscalar.nii'];
    ciftisavereset(gEi_pos_cii,gEi_pos_fname)
    ciftisavereset(gEi_abs_cii,gEi_abs_fname)
    % Set _Glasser_360.pscalar.nii map name
    command = [path_wb_command, ' -set-map-names ', gEi_pos_fname, ' -name-file ', NetSpar_fname];
    system(command)
    command = [path_wb_command, ' -set-map-names ', gEi_abs_fname, ' -name-file ', NetSpar_fname];
    system(command)     
    clear gEi_pos_cii gEi_abs_cii gEi_pos_fname gEi_abs_fname command    
    disp(['----------  Global Efficiency for ', num2str(sbj(sub)), ' Finished  ----------', newline, newline])  
 
 
    %% Time Count
    t2 = clock;
    disp(datestr(now))  
    disp(['Elapsed ',num2str(etime(t2,t1)/60),' min'])
    disp(['...................',num2str(sbj(sub)),' Done ...................'])

end
