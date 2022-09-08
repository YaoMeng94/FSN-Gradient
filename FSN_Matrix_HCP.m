%% Initialization
clear; clc;
path = '/Data_HCP_S1200/REST2';
sbj = importdata('/List/HCP_S1200_REST_4-session_Bad25_sbjlist.txt');
% path_wb_command = '/usr/local/workbench-linux64-v1.4.2/workbench/bin_linux64/wb_command';

%% Setup
% Parcellation Template
File_parcel = '/Templates/Glasser_360.dlabel.nii';
ParcelHeader = ciftiopen(File_parcel);
Parcel = ParcelHeader.cdata;
% Features name
Features_abs_fname = {'ALFF';'fALFF';'Reho';'abs_NodeDegree';'abs_local_efficiency';'abs_global_efficiency'};
Features_pos_fname = {'ALFF';'fALFF';'Reho';'pos_NodeDegree';'pos_local_efficiency';'pos_global_efficiency'};
% .pconn.nii Header for cifti storage %should change for specific Parcellation Template
PconnHeader = ciftiopen('/Templates/sub-MSC01_MyelinMap_biascorrected_Glasser360.pconn.nii'); % .pconn.nii Header for cifti storage %should change for specific Parcellation Template

%% Circulation Body
for sub = 1 : length(sbj)
    t1 = clock;
    disp(['...................', num2str(sbj(sub)),' Starting ...................'])
    disp(datestr(now)) 
    
    %% Features Load
    % with abs network
    Features_abs_Matrix = zeros(length(Parcel(any(unique(Parcel),2))),length(Features_abs_fname));    
    for ftr = 1 : length(Features_abs_fname)
        fname = dir([path, filesep, num2str(sbj(sub)), filesep, 'AllFeatures', filesep, num2str(sbj(sub)), '_*', Features_abs_fname{ftr}, '_32k_fsLR_Glasser_360.pscalar.nii']);
        Ftr_fname = [fname(1).folder, filesep, fname(1).name];
        PscalarHeader = ciftiopen(Ftr_fname);
        Feature = PscalarHeader.cdata;
        Feature = nanmean(Feature, 2);
        Features_abs_Matrix(:,ftr) = zscore(Feature);   
    end    
    % with pos network
    Features_pos_Matrix = zeros(length(Parcel(any(unique(Parcel),2))),length(Features_pos_fname));
    for ftr = 1 : length(Features_pos_fname)
        fname = dir([path, filesep, num2str(sbj(sub)), filesep, 'AllFeatures', filesep, num2str(sbj(sub)), '_*', Features_pos_fname{ftr}, '_32k_fsLR_Glasser_360.pscalar.nii']);
        Ftr_fname = [fname(1).folder, filesep, fname(1).name];
        PscalarHeader = ciftiopen(Ftr_fname);
        Feature = PscalarHeader.cdata;
        Feature = nanmean(Feature, 2);
        Features_pos_Matrix(:,ftr) = zscore(Feature);   
    end 
    
    disp(['...................', num2str(sbj(sub)),' All Features Loaded ...................'])
     
     %% Similarity Matrix Consturction
     % Pearson Correlation
     FSN_pos = corr(Features_pos_Matrix');
     FSN_abs = corr(Features_abs_Matrix');
     % Cosine Distance
     FSN_pos_cosine = 1-squareform(pdist(Features_pos_Matrix,'cosine'));
     FSN_abs_cosine = 1-squareform(pdist(Features_abs_Matrix,'cosine'));
     % Normale Angle
     FSN_pos_angle = 1-acos(FSN_pos_cosine)/pi;
     FSN_abs_angle = 1-acos(FSN_abs_cosine)/pi; 
     
     %% Save matrix
     % correlation
     FSN_pos_cii = PconnHeader;
     FSN_pos_cii.cdata = FSN_pos;
     FSN_abs_cii = PconnHeader;
     FSN_abs_cii.cdata = FSN_abs;
     % cosine distance
     FSN_pos_cosine_cii = PconnHeader;
     FSN_pos_cosine_cii.cdata = FSN_pos_cosine;
     FSN_abs_cosine_cii = PconnHeader;
     FSN_abs_cosine_cii.cdata = FSN_abs_cosine;  
     % normal angle
     FSN_pos_angle_cii = PconnHeader;
     FSN_pos_angle_cii.cdata = FSN_pos_angle;
     FSN_abs_angle_cii = PconnHeader;
     FSN_abs_angle_cii.cdata = FSN_abs_angle;   
     % deatination dir
     corrdir = [path, filesep,  num2str(sbj(sub)), filesep, 'SimilarityMatrix', filesep, 'Correlation'];
     cosinedir = [path, filesep,  num2str(sbj(sub)), filesep, 'SimilarityMatrix', filesep, 'Cosine'];
     angledir = [path, filesep,  num2str(sbj(sub)), filesep, 'SimilarityMatrix', filesep, 'NormalAngle'];
     mkdir(corrdir)
     mkdir(cosinedir)
     mkdir(angledir)
     % save to .pconn.nii
     ciftisavereset(FSN_pos_cii, [corrdir, filesep, num2str(sbj(sub)), '_NeoCortex_pos_Correlation_32k_fsLR_Glasser_360.pconn.nii'])
     ciftisavereset(FSN_abs_cii, [corrdir, filesep, num2str(sbj(sub)), '_NeoCortex_abs_Correlation_32k_fsLR_Glasser_360.pconn.nii'])
     
     ciftisavereset(FSN_pos_cosine_cii, [cosinedir, filesep, num2str(sbj(sub)), '_NeoCortex_pos_Cosine_32k_fsLR_Glasser_360.pconn.nii'])
     ciftisavereset(FSN_abs_cosine_cii, [cosinedir, filesep, num2str(sbj(sub)), '_NeoCortex_abs_Cosine_32k_fsLR_Glasser_360.pconn.nii'])

     ciftisavereset(FSN_pos_angle_cii, [angledir, filesep, num2str(sbj(sub)), '_NeoCortex_pos_NormalAngle_32k_fsLR_Glasser_360.pconn.nii'])
     ciftisavereset(FSN_abs_angle_cii, [angledir, filesep, num2str(sbj(sub)), '_NeoCortex_abs_NormalAngle_32k_fsLR_Glasser_360.pconn.nii'])
     
     t2 = clock;
     disp(datestr(now))  
     disp(['Elapsed ',num2str(etime(t2,t1)/60),' min'])
     disp(['...................', num2str(sbj(sub)),' Similarity Matrix Constructed ...................'])

end