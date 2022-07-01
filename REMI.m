% This script was used to integrate transcriptome and metabolome data for acidic stress conditions to make a context-specific model. 
% This script was modified from tutorial of REMI methods
% Ref.: Pandey V, Hadadi N, Hatzimanikatis V. 2019. Enhanced flux prediction by integrating relative expression 
% and relative metabolite abundance into thermodynamically consistent metabolic models. PLoS Comput Biol 15:e1007036.

addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio128')); % add your CPLEX path into matlab environment
initCobraToolbox

% load iQY486 model
load('E:\bioinformatics\iWH1038_h2o2.mat'); 

% Add relative expression constraint 
% (1) load relative expression data
delimiterIn = '\t';
[A,delimiterOut]=importdata('E:\bioinformatics\REMI_transcriptome.txt');
genedata = importdata('E:\bioinformatics\REMI_transcriptome.txt', delimiterIn);
genename = genedata.textdata(:,1);
ratios = genedata.data(:,1);
up_cutoff=1.42; % fold-change>2 was regraded as up regulated genes
down_cutoff=0.7;% fold-change < 1/2 was regraded as down regulated genes 
% we need to identify genes which are up and down regulated
% input argument : 1) model
%                  2) gene names
%                  3) gene ratios
%                  4) and opertaor (evalute based on GPR)
%                  5) or operator
% length of argument 2 and 3 should be same
[rxns,ratio]=evaluateGPR(Azorhizobium_caulinodans_ORS_571,genename,ratios,@min,@max);
% find up- and down- regulated reactions
indUP=find(ratio>up_cutoff);
ratioUP=ratio(indUP);
indDOWN=find(ratio<down_cutoff);
ratioDOWN=ratio(indDOWN);
regRxnRatio=[ratioUP;ratioDOWN];
regRxns=[rxns(indUP)';rxns(indDOWN)'];

% avoid numerical error (more than 100 fold is taken only 100)
regRxnRatio(regRxnRatio>100)=100;
regRxnRatio(regRxnRatio<1/100)=1/100;

% if we want to add relative constraint into TFA model then we need to add
% net flux variable to the model using 'addNetFluxVariablesNEW'
% and if one want to add into FBA model  then evalute scripts given below
mTmp=Azorhizobium_caulinodans_ORS_571;
mTmp.rev=ones(numel(mTmp.rev),1);
use_m=addUseVariablesDH(mTmp);
netM=addNetFluxVariablesNEW(use_m);
netM1=addNetFluxVariablesNEW(use_m);
netM2=addNetFluxVariablesNEW(use_m);

%% We are going to add constraints for only relative expression
% now we add relative expression constarint 
% input argument : 1) model1 represents condition 1
%                  2) model2 represents condition 2 
%                  3) regulated rxns 
%                  4) regulated reaction rations

% now we are build gex Model:  which integeartes relative gene expression
[gex]=addRelConsExpression(netM1,netM2,regRxns,regRxnRatio)


demand_reaction_names0 = {'demand_pydam_c', 'demand_camp_c', 'demand_C21593_c', 'demand_4h2kpi_c', 'demand_24dhhed_c', 'demand_C03827_c','demand_ivcoa_c', 'demand_glutcoa_c'}
demand_reaction_names1 = {'demand_f1p_c', 'demand_dha_c', 'demand_C02814_c', 'demand_frmd_c', 'demand_itacon_c', 'demand_C06554_c'}
demand_reaction_names2 = {'demand_thrnt_c', 'demand_5dh4dglcn_c', 'demand_oxam_c', 'demand_1acpc_c', 'demand_3hmp_c', 'demand_mmal_c'}
demand_reaction_names3 = {'demand_creat_c', 'demand_acmum_c', 'demand_fumpyr_c', 'demand_pydxn_c', 'demand_dgo15lac_c', 'demand_cyst__L_c'}
demand_reaction_names4 = {'demand_Nforglu_c', 'demand_2dh3dgal_c', 'demand_dma_c', 'demand_rb15bp_c', 'demand_cynt_c', 'demand_oca_c','demand_4hbzcoa_c','demand_fe2'}

demand_reaction_names = [demand_reaction_names0,demand_reaction_names1,demand_reaction_names2,demand_reaction_names3,demand_reaction_names4]

[~,ind_demand_reactions]=ismember(strcat('NF_',demand_reaction_names),gex.varNames);
[~,Pind_demand_reactions]=ismember(strcat('PERTURB_NF_',demand_reaction_names),gex.varNames);
gex.var_ub(ind_demand_reactions) = 5
gex.var_lb(ind_demand_reactions) = 0
gex.var_ub(Pind_demand_reactions) = 5
gex.var_lb(Pind_demand_reactions) = 0

sink_reaction_names0 = {'sink_preq1_c','sink_hspmd_c','sink_4hthr_c','sink_4abutn_c','sink_C04226_c','sink_ins_c','sink_dimp_c','sink_C04546_c'}
sink_reaction_names1 = {'sink_uri_c','sink_udpgalur_c','sink_udpgal_c','sink_uacmam_c','sink_uacgamo_c','sink_rnam_c','sink_2dr1p_c','sink_mocogdp_c'}
sink_reaction_names2 = {'sink_3hpp_c','sink_gal14lac_c','sink_gsn_c','sink_nicrns_c','sink_xtsn_c','sink_udpxyl_c','sink_thym_c','sink_2hyoxplac_c'}
sink_reaction_names3 = {'sink_C05165_c','sink_C21723_c','sink_C21724_c','sink_cdp4dh6doglc_c','sink_cytd_c','sink_dad_2_c','sink_dcyt_c','sink_dgsn_c'}
sink_reaction_names4 = {'sink_dialurate_c','sink_dtdprmn_c','sink_gdpdrhmn_c','sink_gdpfuc_c','sink_gdpper_c','sink_hemeA_c','sink_hgentis_c','sink_mococdp_c'}
sink_reaction_names5 = {'sink_ppgpp_c','sink_rhcys_c','sink_sheme_c','sink_thmpp_c'}
sink_reaction_names = [sink_reaction_names0,sink_reaction_names1,sink_reaction_names2,sink_reaction_names3,sink_reaction_names4,sink_reaction_names5]
[~,ind_sink_reactions]=ismember(strcat('NF_',sink_reaction_names),gex.varNames);
[~,Pind_sink_reactions]=ismember(strcat('PERTURB_NF_',sink_reaction_names),gex.varNames);
gex.var_ub(ind_sink_reactions) = 5
gex.var_lb(ind_sink_reactions) = 0
gex.var_ub(Pind_sink_reactions) = 5
gex.var_lb(Pind_sink_reactions) = 0

exchange_reaction_names0 = {'EX_fe3_e','EX_asn__L_e','EX_etha_e','EX_glu__L_e','EX_ncam_e','EX_glyc_e','EX_ade_e','EX_rhmn_e','EX_fum_e','EX_glyc3p_e'}
exchange_reaction_names1 = {'EX_cit_e','EX_btd_RR_e','EX_acon_C_e','EX_lac__L_e','EX_ac_e','EX_no3_e','EX_4hbz_e','EX_akg_e','EX_ala__L_e','EX_urea_e','EX_arg__L_e'}
exchange_reaction_names2 = {'EX_glcur_e','EX_no2_e','EX_mepn_e','EX_nh4_e','EX_hxan_e','EX_xan_e','EX_gua_e','EX_meoh_e','EX_h2s_e','EX_so4_e','EX_oxa_e','EX_34dhbz_e'}
exchange_reaction_names3 = {'EX_mobd_e','EX_bhb_e','EX_mal__L_e','EX_ppa_e','EX_gthrd_e','EX_starch_e','EX_mso3_e','EX_glyc__R_e','EX_for_e','EX_chol_e','EX_galct__D_e'}
exchange_reaction_names4 = {'EX_glyclt_e','EX_tsul_e','EX_co_e','EX_met__L_e','EX_h2_e','EX_urate_e','EX_btn_e','EX_ddca_e','EX_dtbt_e','EX_hdca_e','EX_hdcea_e'}
exchange_reaction_names5 = {'EX_galur_e','EX_ch4_e','EX_malon_e','EX_5dglcn_e','EX_ocdcea_e','EX_ocdca_e','EX_ttdcea_e','EX_ttdca_e','EX_asp__L_e','EX_galct__D_e','EX_agm_e','EX_n2_e'}

exchange_reaction_names = [exchange_reaction_names0,exchange_reaction_names1,exchange_reaction_names2,exchange_reaction_names3,exchange_reaction_names4,exchange_reaction_names5]
% exchange bounds for gas and microelements
exchange_reaction_micro = {'EX_pi_e','EX_o2_e','EX_h_e','EX_ni2_e','EX_na1_e','EX_cl_e','EX_mg2_e','EX_cobalt2_e','EX_k_e','EX_zn2_e','EX_mn2_e','EX_nh4_e','EX_co2_e','EX_pyr_e'}

[~,ind_ex_reactions]=ismember(strcat('NF_',exchange_reaction_names),gex.varNames);
[~,ind_micro]=ismember(strcat('NF_',exchange_reaction_micro),gex.varNames);
[~,Pind_ex_reactions]=ismember(strcat('PERTURB_NF_',exchange_reaction_names),gex.varNames);
[~,Pind_micro]=ismember(strcat('PERTURB_NF_',exchange_reaction_micro),gex.varNames);
gex.var_ub(ind_ex_reactions) = 20
gex.var_lb(ind_ex_reactions) = -20
gex.var_ub(ind_micro) = 100
gex.var_lb(ind_micro) = -100
gex.var_ub(Pind_ex_reactions) = 20
gex.var_lb(Pind_ex_reactions) = -20
gex.var_ub(Pind_micro) = 100
gex.var_lb(Pind_micro) = -100

[~,ind_ACONT]=ismember('NF_ACONT',gex.varNames);
[~,Pind_ACONT]=ismember('PERTURB_NF_ACONT',gex.varNames);
gex.var_ub(ind_ACONT) = 0
gex.var_lb(ind_ACONT) = 0
gex.var_ub(Pind_ACONT) = 0
gex.var_lb(Pind_ACONT) = 0

% Force directions
% second round: PGI
fw_reactions = {'FUM','ICDHyr','MDH','PGI','DDCAt','TTDCAt','HDCAt','TTDCEAt','HDCEAt','OCDCAt','OCDCEAt','SUCOAACTr'}
bw_reactions = {'SUCOAS','PYRt2','ALCD1'}
[~,ind_fw]=ismember(strcat('NF_',fw_reactions),gex.varNames);
[~,Pind_fw]=ismember(strcat('PERTURB_NF_',fw_reactions),gex.varNames);
[~,ind_bw]=ismember(strcat('NF_',bw_reactions),gex.varNames);
[~,Pind_bw]=ismember(strcat('PERTURB_NF_',bw_reactions),gex.varNames);
gex.var_lb(ind_fw) = 0
gex.var_lb(Pind_fw) = 0
gex.var_ub(ind_bw) = 0
gex.var_ub(Pind_bw) = 0

% Infutile cycles
% HACD1;HACD1y|KARA1;KARI_3hmoa;KARI_23dhmb|ACKr;SUCOAACTr;PTAr|ECOAH1R;HACD1R;ECOAH1|
loop1 = {'HACD1','ALATA_L','ACKr','ECOAH1R'}
loop1_ir = {'KARI_3hmoa','THRA2i'}
[~,ind_loop1]=ismember(strcat('NF_',loop1),gex.varNames);
[~,Pind_loop1]=ismember(strcat('PERTURB_NF_',loop1),gex.varNames);
[~,ind_loop1_ir]=ismember(strcat('NF_',loop1_ir),gex.varNames);
[~,Pind_loop1_ir]=ismember(strcat('PERTURB_NF_',loop1_ir),gex.varNames);
gex.var_lb(ind_loop1) = -15
gex.var_lb(Pind_loop1) = -15
gex.var_ub(ind_loop1) = 15
gex.var_ub(Pind_loop1) = 15
gex.var_ub(ind_loop1_ir) = 15
gex.var_ub(Pind_loop1_ir) = 15

% H2O2 invade

sol_gex=solveTFBAmodel(gex)

% maximum consistency score (MCS) will be objective value of the solution.
MCS=sol_gex.val; % this is the maximum consistency score of GeX model 

%% ALTERNATIVE ANALYSIS on a given MCS
% enumerate alternatives on expression variabels or for each MCS score
% there can be many alternative set of constraints and if one want to
% enumerate those alternatives one should use follwing script


path_save='E:\bioinformatics\0.01_alt.mat' % This path will save all alternative results
% all binary variables are stored as model.relExp.forB
%                  3) will be model.relExp.forB

% we want to find alternative at maximum consistency then we need to add
% maximum consistency to the lowerbound
coM2=gex;
coM2.var_lb(end)=MCS; % this will force maximum consistency can not be 
time=30;
% input argument : 1) numsol is number of alternatives
%                  2) model
%                  3) index for binary variables 
%                  4) path for saving the resultes
%                  5) time for solver
findAltCombi(time,coM2,coM2.relExp.forB,path_save);

%% ADD RELATIVE METABOLITE on M and Gex Model
%
% first load relative data
delimiterIn = '\t';
[A,delimiterOut]=importdata('E:\bioinformatics\REMI_h2o2_metabolite.txt');
metdata = importdata('E:\bioinformatics\REMI_h2o2_metabolite.txt', delimiterIn);
metname = metdata.textdata(:,1);
mratio = metdata.data(:,1);

indUP=find(mratio>1.2);
indDOWN=find(mratio<0.9);

regMetRatio=[mratio(indUP);mratio(indDOWN)];
regMets=[metname(indUP);metname(indDOWN)];
%% We are going to add constraints for only relative metabolites (M model)
% now we add relative expression constarint 
% input argument : 1) model1 represents condition 1
%                  2) model2 represents condition2 
%                  3) regulated mets 
%                  4) regulated metratios
%                  5) false (for only metaboliet), true(if expression is alreay intgerated)
[M]=addRelMetabolite(netM,netM,regMets,regMetRatio,false)
sol_M=solveTFBAmodel(M) 

% you can force consistency by putting MCS in the lower bound. 
% coM2.var_lb(end)=MCS; % this will force maximum consistency can not be 


% find alternatives

findAltCombi(2,M,M.metB,path_save,time);

%% We want to add relative expression and metabolites together (GexM model)

[gexM]=addRelMetabolite(gex,gex,regMets,regMetRatio,true);
sol=solveTFBAmodel(gexM)
% find alternatives
path_save='E:\bioinformatics\ALTERNATIVE_for_0.01mM_new_418.mat'
comIdx=[gexM.relExp.forB;gexM.metB];
findAltCombi(10,gexM,comIdx,path_save);

%% force model for one alternative solution

modelTMP=gexM;
sol1=solveTFBAmodel(modelTMP)
index1=comIdx(find(sol.x(comIdx)>0.98));
index0=comIdx(find(sol.x(comIdx)<0.08));
modelTMP.var_lb(index1)=1;
modelTMP.var_ub(index0)=0;
sol2=solveTFBAmodel(modelTMP);




