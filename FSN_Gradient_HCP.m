%% 
clear; clc;
path = '/Data_HCP_S1200/REST2';
sbj = importdata('/List/HCP_S1200_REST_4-session_excluBad25_sbjlist.txt');
% path_wb_command = '/usr/local/workbench-linux64-v1.4.2/workbench/bin_linux64/wb_command';

%% Parameter Setup
Thre_prctile = 0; % global threshold of Sim Mat
Dist = {'NormalAngle'};
Pscalar = ciftiopen('/Templates/sub-MSC01_MyelinMap_biascorrected_Glasser360.pscalar.nii');
parcels = Pscalar.cdata;
n_grad = 10;
nIterations = 10000;

%% Circulation Body

% Merge all sub's Sim Net (for average Sim Net and later Realign template)
All_SimNet = zeros(size(parcels, 1), size(parcels, 1), length(sbj)); 
All_Embedding = cell(length(sbj), 1);

% subject loop
for sub = 1 : length(sbj)
    t1 = clock;
    disp(['...................', num2str(sbj(sub)),' FSN Gradient Starting ...................'])
    disp(datestr(now)) 
       
    %% Postive Network   
        % load Similarity Matrix
        FSN_fname_pos = [path, filesep, num2str(sbj(sub)), filesep, 'SimilarityMatrix', filesep, Dist{1}, filesep, num2str(sbj(sub)), '_NeoCortex_pos_', Dist{1}, '_32k_fsLR_Glasser_360.pconn.nii'];
        FSNHeader = ciftiopen(FSN_fname_pos);
        FSN_pos = FSNHeader.cdata;
        % Threshold the matrix
        Threshld = prctile(tril2vec(FSN_pos),Thre_prctile);
        FSN_pos_threshlded = FSN_pos .* (FSN_pos > Threshld);
        All_SimNet(:, :, sub) = FSN_pos_threshlded;
        % diffusion map embedding
        [embedding, result] = mica_diffusionEmbedding(FSN_pos_threshlded, 'nComponents', n_grad);
        All_Embedding{sub} = embedding;
    % save embedded result 
    % .mat
    subdir = [path, filesep, num2str(sbj(sub)), filesep, 'FSN_gradient',filesep, 'pos'];
    mkdir(subdir);
    save([subdir, filesep, num2str(sbj(sub)), '_Embedding'], 'embedding')
    save([subdir, filesep, num2str(sbj(sub)), '_Result'], 'result')   
    clear embedding result FSN_pos FSN_pos_threshlded
    
    t2 = clock;
    disp(datestr(now))  
    disp(['Elapsed ',num2str(etime(t2,t1)/60),' min'])
    disp(['...................', num2str(sbj(sub)),' FSN Gradient Done ...................'])
    fprintf('\n\n')

end

% Average Sim Net
Mean_SimNet = mean(All_SimNet, 3);
[GrpMean_embedding, GrpMean_result] = mica_diffusionEmbedding(Mean_SimNet, 'nComponents', n_grad);
Grpdir = [path, filesep, 'Group', filesep, 'FSN_gradient',filesep, 'pos'];
mkdir(Grpdir);
save([Grpdir, filesep, 'GroupMeanSim_Embedding'], 'GrpMean_embedding')
save([Grpdir, filesep, 'GroupMeanSim_Result'], 'GrpMean_result')
RealignTemp = Pscalar;
RealignTemp.cdata = GrpMean_embedding;
RealignTemp_fname = [Grpdir, filesep, 'RealignTemplate_pos.pscalar.nii'];
ciftisavereset(RealignTemp, RealignTemp_fname)

[All_Embedding_Realigned, xfms]= mica_iterativeAlignment(All_Embedding, nIterations, GrpMean_embedding);
save([Grpdir, filesep, 'All_Embedding_Realigned'], 'All_Embedding_Realigned')
save([Grpdir, filesep, 'All_Embedding_xfms'], 'xfms')


for sub = 1 : length(sbj)
    
    SubGrad_cii = Pscalar;
    SubGrad_cii.cdata = squeeze(All_Embedding_Realigned(:, :, sub));
    subdir = [path, filesep, num2str(sbj(sub)), filesep, 'FSN_gradient',filesep, 'pos'];
    SubGrad_fname = [subdir, filesep, num2str(sbj(sub)), '_pos_gradient_32k_fsLR_Glasser_360.pscalar.nii'];
    ciftisavereset(SubGrad_cii, SubGrad_fname)
    
end
