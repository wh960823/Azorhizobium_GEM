load('E:\bioinformatics\iWH1038_for_h2o2_1.mat'); 
model = readCbModel('E:\bioinformatics\new_Azo_for_h2o2_1.mat')
model.rev = Azorhizobium_caulinodans_ORS_571.rev
FBA = optimizeCbModel(Azorhizobium_caulinodans_ORS_571);
"""
[MinimizedFlux, modelIrrev]= minimizeModelFlux(model);
if isfield(FBA, 'obj')==1
    obj = FBA.obj;
else isfield(FBA,'f')== 1
    obj = FBA.f;
end
modelIrrev.lb(find(model.c)) = obj;
modelIrrev.ub(find(model.c)) = obj;
minflux = optimizeCbModel(modelIrrev,'min');
if isfield(minflux, 'full')==1
    vec = minflux.full;
else isfield(minflux,'x')== 1
    vec = minflux.x;
end
maxflux_val = max(abs(vec(1:(size(model.rxns,1)+numel(find(model.rev))))));
"""
[model_del, nochange_idx] = createDeltaModel(model, [914], [],1000)
model_delNet = addNetDeltaFlux(model_del, model.rxns)
delimiterIn = '\t';
[A,delimiterOut]=importdata('E:\bioinformatics\REMI_transcriptome.txt');
genedata = importdata('E:\bioinformatics\REMI_transcriptome.txt', delimiterIn);
genename = genedata.textdata(:,1);
ratios = genedata.data(:,2);
expressionData.gene = model.genes;
expressionData.value = zeros(numel(model.genes),1);
for i = 1:numel(model.genes)
    [~,ind]=ismember(model.genes{i},genename);
    if ind>0
        expressionData.value(i)=ratios(ind);
    end
end
RxnExpr = mapExpressionToReactions(model,expressionData);
RxnExpr(find(RxnExpr==(-1))) = 0;
[modelIrrev, matchRev, rev2irrev, irrev2rev] = convertToIrreversible(model);
de_exprrxns = mapExpressionToReactions(modelIrrev,expressionData);
de_rxns = modelIrrev.rxns;
sel = find(~(de_exprrxns==(-1) | de_exprrxns==0));
de_exprrxns = de_exprrxns(sel);
de_rxns = de_rxns(sel);
de_indUP = find(de_exprrxns>1.42)
de_indDOWN = find(de_exprrxns<0.7)
%cut_off = 0.25;
%[~,de_indUP]=maxk(de_exprrxns,round(cut_off*size(find(de_exprrxns>1),1)));
%[~,de_indDOWN]=mink(de_exprrxns,round(cut_off*size(find(de_exprrxns<1),1)));
regRxns=[de_rxns(de_indUP); de_rxns(de_indDOWN)];
regRxnsRatio = [de_exprrxns(de_indUP) ; de_exprrxns(de_indDOWN)];
%model_binary = createBinaryUseVariable(model_delNet, 2, 1e5, regRxns,regRxnRatio)
regRxnsRatio(find(abs(regRxnsRatio)<0.1)) = 0.1;
regRxnsRatio(find(abs(regRxnsRatio)>100)) = 100;
model_Bin = createBinaryUseVariable(model_delNet, 2, 1e5, regRxns, ones(numel(regRxns),1));
[MILPstructure, model_binary] = createRxnConsistencyObj(model_Bin,de_rxns,de_indUP, de_indDOWN, nochange_idx, 1, de_exprrxns, 0);
[~,NFind]=ismember(strcat('NF_',model.rxns),MILPstructure.varNames);


[~,NFind_NIT]=ismember('NF_NIT',MILPstructure.varNames);
[~,NFind_Symbiosis]=ismember('NF_Symbiosis',MILPstructure.varNames);
demand_reaction_names0 = {'demand_pydam_c', 'demand_camp_c', 'demand_C21593_c', 'demand_4h2kpi_c', 'demand_24dhhed_c', 'demand_C03827_c','demand_ivcoa_c', 'demand_glutcoa_c'}
demand_reaction_names1 = {'demand_f1p_c', 'demand_dha_c', 'demand_C02814_c', 'demand_frmd_c', 'demand_itacon_c', 'demand_C06554_c'}
demand_reaction_names2 = {'demand_thrnt_c', 'demand_5dh4dglcn_c', 'demand_oxam_c', 'demand_1acpc_c', 'demand_3hmp_c', 'demand_mmal_c'}
demand_reaction_names3 = {'demand_creat_c', 'demand_acmum_c', 'demand_fumpyr_c', 'demand_pydxn_c', 'demand_dgo15lac_c', 'demand_cyst__L_c'}
demand_reaction_names4 = {'demand_Nforglu_c', 'demand_2dh3dgal_c', 'demand_dma_c', 'demand_rb15bp_c', 'demand_cynt_c', 'demand_oca_c','demand_4hbzcoa_c','demand_fe2'}

demand_reaction_names = [demand_reaction_names0,demand_reaction_names1,demand_reaction_names2,demand_reaction_names3,demand_reaction_names4]
[~,NFind_demand_reactions]=ismember(strcat('NF_',demand_reaction_names),MILPstructure.varNames);

% To limit exchange bounds
exchange_reaction_names0 = {'EX_fe3_e','EX_asn__L_e','EX_etha_e','EX_glu__L_e','EX_ncam_e','EX_glyc_e','EX_ade_e','EX_rhmn_e','EX_fum_e','EX_glyc3p_e'}
exchange_reaction_names1 = {'EX_cit_e','EX_btd_RR_e','EX_acon_C_e','EX_lac__L_e','EX_ac_e','EX_no3_e','EX_4hbz_e','EX_akg_e','EX_ala__L_e','EX_urea_e','EX_arg__L_e'}
exchange_reaction_names2 = {'EX_glcur_e','EX_no2_e','EX_mepn_e','EX_nh4_e','EX_hxan_e','EX_xan_e','EX_gua_e','EX_meoh_e','EX_h2s_e','EX_so4_e','EX_oxa_e','EX_34dhbz_e'}
exchange_reaction_names3 = {'EX_mobd_e','EX_bhb_e','EX_pyr_e','EX_mal__L_e','EX_ppa_e','EX_gthrd_e','EX_starch_e','EX_mso3_e','EX_glyc__R_e','EX_for_e','EX_chol_e','EX_galct__D_e'}
exchange_reaction_names4 = {'EX_glyclt_e','EX_tsul_e','EX_co_e','EX_met__L_e','EX_h2_e','EX_urate_e','EX_btn_e','EX_ddca_e','EX_dtbt_e','EX_hdca_e','EX_hdcea_e'}
exchange_reaction_names5 = {'EX_galur_e','EX_ch4_e','EX_malon_e','EX_5dglcn_e','EX_ocdcea_e','EX_ocdca_e','EX_ttdcea_e','EX_ttdca_e','EX_asp__L_e','EX_galct__D_e','EX_agm_e','EX_n2_e'}

exchange_reaction_names = [exchange_reaction_names0,exchange_reaction_names1,exchange_reaction_names2,exchange_reaction_names3,exchange_reaction_names4,exchange_reaction_names5]
% exchange bounds for gas and microelements
exchange_reaction_micro = {'EX_pi_e','EX_o2_e','EX_h_e','EX_ni2_e','EX_na1_e','EX_cl_e','EX_mg2_e','EX_cobalt2_e','EX_k_e','EX_zn2_e','EX_mn2_e','EX_nh4_e','EX_co2_e'}


[~,NFind_ex_reactions]=ismember(strcat('NF_',exchange_reaction_names),MILPstructure.varNames);
[~,NFind_ex_reactions_micro]=ismember(strcat('NF_',exchange_reaction_micro),MILPstructure.varNames);
MILPstructure.ub(NFind_ex_reactions) = 10
MILPstructure.lb(NFind_ex_reactions) = -10
MILPstructure.ub(NFind_ex_reactions_micro) = 100
MILPstructure.lb(NFind_ex_reactions_micro) = -100

[~,NFind_demand_h2o2_c]=ismember('NF_demand_h2o2_c',MILPstructure.varNames);
[~,NFind_EX_o2_e]=ismember('NF_EX_o2_e',MILPstructure.varNames);
[~,NFind_EX_succ_e]=ismember('NF_EX_succ_e',MILPstructure.varNames);
[~,NFind_EX_h2o_e]=ismember('NF_EX_h2o_e',MILPstructure.varNames);
[~,NFind_NH3c]=ismember('NF_NH3c',MILPstructure.varNames);
[~,NFind_NO3t7]=ismember('NF_NO3t7',MILPstructure.varNames);
MILPstructure.ub(NFind_NO3t7) = 0
MILPstructure.lb(NFind_NO3t7) = 0


%To limit the fluxes of infutile cycles
[~,NFind_ADSK]=ismember('NF_ADSK',MILPstructure.varNames);
[~,NFind_ADNK1]=ismember('NF_ADNK1',MILPstructure.varNames);
[~,NFind_FTHFCL]=ismember('NF_FTHFCL',MILPstructure.varNames);
[~,NFind_AMPTASECG]=ismember('NF_AMPTASECG',MILPstructure.varNames);
[~,NFind_GTHRDH]=ismember('NF_GTHRDH',MILPstructure.varNames);
[~,NFind_TRE6PP]=ismember('NF_TRE6PP',MILPstructure.varNames);
[~,NFind_ASPO2y]=ismember('NF_ASPO2y',MILPstructure.varNames);
[~,NFind_P5CR]=ismember('NF_P5CR',MILPstructure.varNames);
[~,NFind_TPI]=ismember('NF_TPI',MILPstructure.varNames);
[~,NFind_GDPTPDP]=ismember('NF_GDPTPDP',MILPstructure.varNames);
[~,NFind_R01583]=ismember('NF_R01583',MILPstructure.varNames);
[~,NFind_MTHFR3]=ismember('NF_MTHFR3',MILPstructure.varNames);
[~,NFind_HCO3E]=ismember('NF_HCO3E',MILPstructure.varNames);
[~,NFind_FTHFD]=ismember('NF_FTHFD',MILPstructure.varNames);
[~,NFind_ATPHs]=ismember('NF_ATPHs',MILPstructure.varNames);
[~,NFind_GTPHs]=ismember('NF_GTPHs',MILPstructure.varNames);
[~,NFind_ADK2]=ismember('NF_ADK2',MILPstructure.varNames);
[~,NFind_GMPS2]=ismember('NF_GMPS2',MILPstructure.varNames);
[~,NFind_DPCOAPP]=ismember('NF_DPCOAPP',MILPstructure.varNames);
[~,NFind_PGK]=ismember('NF_PGK',MILPstructure.varNames);
[~,NFind_TMN]=ismember('NF_TMN',MILPstructure.varNames);
[~,NFind_NTPP4]=ismember('NF_NTPP4',MILPstructure.varNames);
[~,NFind_NTPP6]=ismember('NF_NTPP6',MILPstructure.varNames);
[~,NFind_NTPP2]=ismember('NF_NTPP2',MILPstructure.varNames);
[~,NFind_NTPP1]=ismember('NF_NTPP1',MILPstructure.varNames);
[~,NFind_ATHRDHr]=ismember('NF_ATHRDHr',MILPstructure.varNames);
[~,NFind_UDPGP]=ismember('NF_UDPGP',MILPstructure.varNames);
[~,NFind_TRPS1]=ismember('NF_TRPS1',MILPstructure.varNames);
[~,NFind_FORCT]=ismember('NF_FORCT',MILPstructure.varNames);
[~,NFind_ACODA]=ismember('NF_ACODA',MILPstructure.varNames);
[~,NFind_R11726]=ismember('NF_R11726',MILPstructure.varNames);
[~,NFind_ADC]=ismember('NF_ADC',MILPstructure.varNames);
[~,NFind_NADN]=ismember('NF_NADN',MILPstructure.varNames);
[~,NFind_PPDK]=ismember('NF_PPDK',MILPstructure.varNames);
[~,NFind_ECOAH1R]=ismember('NF_ECOAH1R',MILPstructure.varNames);
[~,NFind_HACD1R]=ismember('NF_HACD1R',MILPstructure.varNames);
[~,NFind_HACD1]=ismember('NF_HACD1',MILPstructure.varNames);
[~,NFind_HACD1y]=ismember('NF_HACD1y',MILPstructure.varNames);
[~,NFind_GLUR]=ismember('NF_GLUR',MILPstructure.varNames);
[~,NFind_ME1]=ismember('NF_ME1',MILPstructure.varNames);
[~,NFind_ASNS1]=ismember('NF_ASNS1',MILPstructure.varNames);
[~,NFind_THRD_L]=ismember('NF_THRD_L',MILPstructure.varNames);
[~,NFind_G5SADs]=ismember('NF_G5SADs',MILPstructure.varNames);
[~,NFind_CTPS2]=ismember('NF_CTPS2',MILPstructure.varNames);
[~,NFind_CTPS1]=ismember('NF_CTPS1',MILPstructure.varNames);
[~,NFind_SBP]=ismember('NF_SBP',MILPstructure.varNames);
[~,NFind_CYSS]=ismember('NF_CYSS',MILPstructure.varNames);
[~,NFind_XPPT]=ismember('NF_XPPT',MILPstructure.varNames);
[~,NFind_ADK1]=ismember('NF_ADK1',MILPstructure.varNames);
[~,NFind_AHSERL4]=ismember('NF_AHSERL4',MILPstructure.varNames);
[~,NFind_PPC]=ismember('NF_PPC',MILPstructure.varNames);
[~,NFind_UREA]=ismember('NF_UREA',MILPstructure.varNames);
[~,NFind_GLYCLTDy]=ismember('NF_GLYCLTDy',MILPstructure.varNames);
[~,NFind_ACKr]=ismember('NF_ACKr',MILPstructure.varNames);
[~,NFind_DAAD1]=ismember('NF_DAAD1',MILPstructure.varNames);
[~,NFind_HPYRI]=ismember('NF_HPYRI',MILPstructure.varNames);
[~,NFind_G3PD2]=ismember('NF_G3PD2',MILPstructure.varNames);
[~,NFind_KARI_3hmoa]=ismember('NF_KARI_3hmoa',MILPstructure.varNames);
[~,NFind_SUCOAS]=ismember('NF_SUCOAS',MILPstructure.varNames);
[~,NFind_PSP_L]=ismember('NF_PSP_L',MILPstructure.varNames);
[~,NFind_NO3t2]=ismember('NF_NO3t2',MILPstructure.varNames);
[~,NFind_ASR2]=ismember('NF_ASR2',MILPstructure.varNames);
[~,NFind_ACONT]=ismember('NF_ACONT',MILPstructure.varNames);
[~,NFind_CITt2r]=ismember('NF_CITt2r',MILPstructure.varNames);
[~,NFind_CAT]=ismember('NF_CAT',MILPstructure.varNames);
[~,NFind_R03105]=ismember('NF_R03105',MILPstructure.varNames);
[~,NFind_GLUN]=ismember('NF_GLUN',MILPstructure.varNames);
[~,NFind_Growth]=ismember('NF_Growth',MILPstructure.varNames);
[~,NFind_FACOAE180]=ismember('NF_FACOAE180',MILPstructure.varNames);
[~,NFind_FACOAE140]=ismember('NF_FACOAE140',MILPstructure.varNames);
[~,NFind_FACOAE120]=ismember('NF_FACOAE120',MILPstructure.varNames);
[~,NFind_FACOAE161]=ismember('NF_FACOAE161',MILPstructure.varNames);
[~,NFind_FACOAE181]=ismember('NF_FACOAE181',MILPstructure.varNames);
[~,NFind_FACOAE141]=ismember('NF_FACOAE141',MILPstructure.varNames);
[~,NFind_FACOAE160]=ismember('NF_FACOAE160',MILPstructure.varNames);
[~,NFind_MDH]=ismember('NF_MDH',MILPstructure.varNames);
[~,NFind_RNTR1]=ismember('NF_RNTR1',MILPstructure.varNames);
[~,NFind_RNTR2]=ismember('NF_RNTR2',MILPstructure.varNames);
[~,NFind_RNTR3]=ismember('NF_RNTR3',MILPstructure.varNames);
[~,NFind_RNTR4]=ismember('NF_RNTR4',MILPstructure.varNames);
[~,NFind_NADTRHD]=ismember('NF_NADTRHD',MILPstructure.varNames);
[~,NFind_SPODM]=ismember('NF_SPODM',MILPstructure.varNames);
[~,NFind_GLYCLTDx]=ismember('NF_GLYCLTDx',MILPstructure.varNames);
[~,NFind_NDPK2]=ismember('NF_NDPK2',MILPstructure.varNames);
[~,NFind_MCD]=ismember('NF_MCD',MILPstructure.varNames);
[~,NFind_NDPK1]=ismember('NF_NDPK1',MILPstructure.varNames);
[~,NFind_NDPK3]=ismember('NF_NDPK3',MILPstructure.varNames);
[~,NFind_PIt2r]=ismember('NF_PIt2r',MILPstructure.varNames);
[~,NFind_PRDX]=ismember('NF_PRDX',MILPstructure.varNames);
%,,,

direaction_reactions = {'ECOAH7_f','ECOAH6_f','ECOAH5_f','ECOAH4_f','ECOAH3_f','ECOAH2_f','ECOAH1_f','HXANt2r_b','ADEt2r_b','O2t_b','HACD1_f','HACD7_f','HACD6_f','HACD2_f','HACD3_f','FORCT_b','PIt2r_f','starcht2_b','BHBt2_b','L_LACt2r_b','GALCTt2r_b','ACt','ACACT6r_f','ACACT7r_f','ACACT3r_f','ACACT4r_f','CH4t_f','5DGLCNt2r_f','5DGLCNR_f','HPYRI_f','DB4PS_b','E4PD_b','PERD_b','OHPBAT_b','FBA3_f','GLUR_f','ALAR_b','MDH_b','GLA_b','GALURt2r_b','GLCURt2r_b','RHMNt2_b','AICART_b','T2DECAI_b','CYSTA_b','CRTNh_f','GDPPS_b','MGCH_b','PHETA1_f','G5SADs_b','TYRTA_f','ADK2_b','DADK_f','DURIPP_b','TMDPP_b','CYTK2_f','DAPE_b','PPAKr_b','PRPPS_b','MG2t_b','COBALT2t_b','CLt3r_2_b','N2t_b','UDPG4E_b','PGAMT_f','UAG2E_b','PTRCTA2_b','MALEI_b','ASPR_b','NH3c_f','R03925_b','MTI_b','PGMT_f','FBA_f','G1PCTYT_b','R09283_b','R05076_b','otnI_b','ACCAH_b','DKI_b','CITCIb_b','CITCIa_b','ATHRDHr_b','ORPT_f','DHORTS_f','GLYCt_b','SUCOAS_f','SDPTA_f','DURADx_f','A5PISO_b','MAN6PI_f','PMANM_f','GUI1_b','MANAO_f','ALAALAr_b','PSERT_b','GHMT2r_b','MMM2_f','NMGS_b','PRMICI_b','GUI2_b','TAGURr_f','AIRC3_f','R09979_b','ACTD_b','ENO_b','IPPMIb_f','IPPMIa_f','PGK_b','PSCVT_b','DDPGALA_b','MME_f','VALTA_f','IMPC_f','HSERTA_b','AHSERL4_b','FUM_b','ICDHyr_b','AGPR_f','ACOTA_f','ACONTa_b','ACONTb_b','ARGSL_b','ILETA_f'}
[~,NFind_dr]=ismember(direaction_reactions,MILPstructure.varNames);
MILPstructure.ub(NFind_dr) = 0
MILPstructure.lb(NFind_dr) = 0


[~,NFind_METSOXR2]=ismember('NF_METSOXR2',MILPstructure.varNames);
[~,NFind_METSOXR1]=ismember('NF_METSOXR1',MILPstructure.varNames);
[~,NFind_NAt3_1]=ismember('NF_NAt3_1',MILPstructure.varNames);
[~,NFind_MTHFD]=ismember('NF_MTHFD',MILPstructure.varNames);



MILPstructure.ub(NFind_Growth) = -0.001
%MILPstructure.lb(NFind_Growth) = -10
MILPstructure.ub(NFind_MTHFD) = 10
MILPstructure.lb(NFind_MTHFD) = -10
MILPstructure.ub(NFind_NAt3_1) = 20
MILPstructure.lb(NFind_NAt3_1) = -20
MILPstructure.ub(NFind_METSOXR2) = 20
MILPstructure.lb(NFind_METSOXR2) = -20
MILPstructure.ub(NFind_METSOXR1) = 20
MILPstructure.lb(NFind_METSOXR1) = -20
MILPstructure.ub(NFind_PRDX) = 20
MILPstructure.lb(NFind_PRDX) = -20
MILPstructure.ub(NFind_PIt2r) = 100
MILPstructure.lb(NFind_PIt2r) = -100
MILPstructure.ub(NFind_MCD) = 10
MILPstructure.lb(NFind_MCD) = -10
MILPstructure.ub(NFind_GLYCLTDx) = 10
MILPstructure.lb(NFind_GLYCLTDx) = -10
MILPstructure.ub(NFind_SPODM) = 100
MILPstructure.lb(NFind_SPODM) = -100
MILPstructure.ub(NFind_NADTRHD) = 200
MILPstructure.lb(NFind_NADTRHD) = -200
MILPstructure.ub(NFind_NDPK3) = 20
MILPstructure.lb(NFind_NDPK3) = -20
MILPstructure.ub(NFind_NDPK1) = 20
MILPstructure.lb(NFind_NDPK1) = -20
MILPstructure.ub(NFind_NDPK2) = 20
MILPstructure.lb(NFind_NDPK2) = -20
MILPstructure.ub(NFind_RNTR1) = 20
MILPstructure.lb(NFind_RNTR1) = -20
MILPstructure.ub(NFind_RNTR2) = 20
MILPstructure.lb(NFind_RNTR2) = -20
MILPstructure.ub(NFind_RNTR3) = 20
MILPstructure.lb(NFind_RNTR3) = -20
MILPstructure.ub(NFind_RNTR4) = 20
MILPstructure.lb(NFind_RNTR4) = -20
MILPstructure.ub(NFind_Symbiosis) = 0
MILPstructure.lb(NFind_Symbiosis) = 0
MILPstructure.ub(NFind_FACOAE180) = 5
MILPstructure.lb(NFind_FACOAE180) = -5
MILPstructure.ub(NFind_FACOAE140) = 5
MILPstructure.lb(NFind_FACOAE140) = -5
MILPstructure.ub(NFind_FACOAE120) = 5
MILPstructure.lb(NFind_FACOAE120) = -5
MILPstructure.ub(NFind_FACOAE161) = 5
MILPstructure.lb(NFind_FACOAE161) = -5
MILPstructure.ub(NFind_FACOAE181) = 5
MILPstructure.lb(NFind_FACOAE181) = -5
MILPstructure.ub(NFind_FACOAE141) = 5
MILPstructure.lb(NFind_FACOAE141) = -5
MILPstructure.ub(NFind_FACOAE160) = 5
MILPstructure.lb(NFind_FACOAE160) = -5

MILPstructure.ub(NFind_MDH) = 100
MILPstructure.lb(NFind_MDH) = -100

MILPstructure.ub(NFind_NTPP4) = 10
MILPstructure.lb(NFind_NTPP4) = -10

MILPstructure.ub(NFind_NTPP1) = 10
MILPstructure.lb(NFind_NTPP1) = -10
MILPstructure.ub(NFind_NTPP2) = 10
MILPstructure.lb(NFind_NTPP2) = -10
MILPstructure.ub(NFind_NTPP6) = 10
MILPstructure.lb(NFind_NTPP6) = -10
MILPstructure.ub(NFind_GLUN) = 50
MILPstructure.lb(NFind_GLUN) = -50
MILPstructure.ub(NFind_R03105) = 20
MILPstructure.lb(NFind_R03105) = -20
MILPstructure.ub(NFind_CAT) = 50
MILPstructure.lb(NFind_CAT) = -50
MILPstructure.ub(NFind_ACONT) = 0
MILPstructure.lb(NFind_ACONT) = 0
MILPstructure.ub(NFind_SUCOAS) = 100
MILPstructure.lb(NFind_SUCOAS) = -100
MILPstructure.ub(NFind_CITt2r) = 10
MILPstructure.lb(NFind_CITt2r) = -10
MILPstructure.ub(NFind_NO3t2) = 20
MILPstructure.lb(NFind_NO3t2) = -20
MILPstructure.ub(NFind_G3PD2) = 20
MILPstructure.lb(NFind_G3PD2) = -20
MILPstructure.ub(NFind_PSP_L) = 50
MILPstructure.lb(NFind_PSP_L) = -50
MILPstructure.ub(NFind_DAAD1) = 50
MILPstructure.lb(NFind_DAAD1) = -50
MILPstructure.ub(NFind_ACKr) = 50
MILPstructure.lb(NFind_ACKr) = -50
MILPstructure.ub(NFind_UREA) = 50
MILPstructure.lb(NFind_UREA) = -50
MILPstructure.ub(NFind_KARI_3hmoa) = 10
MILPstructure.lb(NFind_KARI_3hmoa) = -10
MILPstructure.ub(NFind_ASR2) = 5
MILPstructure.lb(NFind_ASR2) = -5
MILPstructure.ub(NFind_FTHFCL) = 5
MILPstructure.lb(NFind_FTHFCL) = -5
MILPstructure.ub(NFind_ADNK1) = 5
MILPstructure.lb(NFind_ADNK1) = -5
MILPstructure.ub(NFind_ADSK) = 5
MILPstructure.lb(NFind_ADSK) = -5
MILPstructure.ub(NFind_AMPTASECG) = 5
MILPstructure.lb(NFind_AMPTASECG) = -5
MILPstructure.ub(NFind_GTHRDH) = 5
MILPstructure.lb(NFind_GTHRDH) = -5
MILPstructure.ub(NFind_TRE6PP) = 5
MILPstructure.lb(NFind_TRE6PP) = -5
MILPstructure.ub(NFind_FTHFD) = 5
MILPstructure.lb(NFind_FTHFD) = -5
MILPstructure.ub(NFind_GLYCLTDy) = 5
MILPstructure.lb(NFind_GLYCLTDy) = -5

%HSST;

loop_reaction_level1 = {'NI2ti','PDX5POi','SSALx','BETALDHy','PYDXOR','NMNAT','HPPK2','HSST','ANS','POAACR','RBCh','GUAD','SBP','ACP1_FMN','NNATr','ALDD2y','GLNS','THRA2','G3PD1','HSDx','FMNRy','H4THDPRy','MECDPDH5','denD','ACYP_2','ACS2'} %level1: [-5,5]
[~,NFind_loop1]=ismember(strcat('NF_',loop_reaction_level1),MILPstructure.varNames);
MILPstructure.ub(NFind_loop1) = 5
MILPstructure.lb(NFind_loop1) = -5
MILPstructure.ub(NFind_AHSERL4) = 20
MILPstructure.lb(NFind_AHSERL4) = -20
MILPstructure.ub(NFind_CYSS) = 20
MILPstructure.lb(NFind_CYSS) = -20
MILPstructure.ub(NFind_CTPS1) = 20
MILPstructure.lb(NFind_CTPS1) = -20
MILPstructure.ub(NFind_CTPS2) = 20
MILPstructure.lb(NFind_CTPS2) = -20
MILPstructure.ub(NFind_G5SADs) = 20
MILPstructure.lb(NFind_G5SADs) = -20
MILPstructure.ub(NFind_THRD_L) = 20
MILPstructure.lb(NFind_THRD_L) = -20
MILPstructure.ub(NFind_PPC) = 20
MILPstructure.lb(NFind_PPC) = -20
MILPstructure.ub(NFind_ADK1) = 100
MILPstructure.lb(NFind_ADK1) = -100
MILPstructure.ub(NFind_XPPT) = 100
MILPstructure.lb(NFind_XPPT) = -100
MILPstructure.ub(NFind_ME1) = 100
MILPstructure.lb(NFind_ME1) = -100
MILPstructure.ub(NFind_ASNS1) = 20
MILPstructure.lb(NFind_ASNS1) = -20
MILPstructure.ub(NFind_GLUR) = 20
MILPstructure.lb(NFind_GLUR) = -20
MILPstructure.ub(NFind_ECOAH1R) = 50
MILPstructure.lb(NFind_ECOAH1R) = -50
MILPstructure.ub(NFind_TRPS1) = 0
MILPstructure.lb(NFind_TRPS1) = 0
MILPstructure.ub(NFind_HACD1y) = 100
MILPstructure.lb(NFind_HACD1y) = -100
MILPstructure.ub(NFind_HACD1) = 100
MILPstructure.lb(NFind_HACD1) = -100
MILPstructure.ub(NFind_HACD1R) = 100
MILPstructure.lb(NFind_HACD1R) = -100
MILPstructure.ub(NFind_HACD1R) = 100
MILPstructure.lb(NFind_HACD1R) = -100
MILPstructure.ub(NFind_HCO3E) = 100
MILPstructure.lb(NFind_HCO3E) = -100
MILPstructure.ub(NFind_PPDK) = 20
MILPstructure.lb(NFind_PPDK) = -20
MILPstructure.ub(NFind_ACODA) = 20
MILPstructure.lb(NFind_ACODA) = -20
MILPstructure.ub(NFind_TMN) = 20
MILPstructure.lb(NFind_TMN) = -20
MILPstructure.ub(NFind_PGK) = 20
MILPstructure.lb(NFind_PGK) = -20
MILPstructure.ub(NFind_NADN) = 10
MILPstructure.lb(NFind_NADN) = -10
MILPstructure.ub(NFind_HPYRI) = 10
MILPstructure.lb(NFind_HPYRI) = -10
MILPstructure.ub(NFind_R11726) = 10
MILPstructure.lb(NFind_R11726) = -10
MILPstructure.ub(NFind_UDPGP) = 10
MILPstructure.lb(NFind_UDPGP) = -10
MILPstructure.ub(NFind_FORCT) = 10
MILPstructure.lb(NFind_FORCT) = -10
MILPstructure.ub(NFind_ATHRDHr) = 10
MILPstructure.lb(NFind_ATHRDHr) = -10
MILPstructure.ub(NFind_DPCOAPP) = 10
MILPstructure.lb(NFind_DPCOAPP) = -10
MILPstructure.ub(NFind_GMPS2) = 10
MILPstructure.lb(NFind_GMPS2) = -10
MILPstructure.ub(NFind_ADK2) = 10
MILPstructure.lb(NFind_ADK2) = -10
MILPstructure.ub(NFind_ATPHs) = 10
MILPstructure.lb(NFind_ATPHs) = -10
MILPstructure.ub(NFind_GTPHs) = 10
MILPstructure.lb(NFind_GTPHs) = -10
MILPstructure.ub(NFind_MTHFR3) = 5
MILPstructure.lb(NFind_MTHFR3) = -5
MILPstructure.ub(NFind_R01583) = 5
MILPstructure.lb(NFind_R01583) = -5
MILPstructure.ub(NFind_ADC) = 5
MILPstructure.lb(NFind_ADC) = -5
MILPstructure.ub(NFind_GDPTPDP) = 5
MILPstructure.lb(NFind_GDPTPDP) = -5
MILPstructure.ub(NFind_TPI) = 20
MILPstructure.lb(NFind_TPI) = -20
MILPstructure.ub(NFind_P5CR) = 10
MILPstructure.lb(NFind_P5CR) = -10
MILPstructure.ub(NFind_ASPO2y) = 10
MILPstructure.lb(NFind_ASPO2y) = -10


%To limit the sink reaction
sink_reaction_names0 = {'sink_preq1_c','sink_hspmd_c','sink_4hthr_c','sink_4abutn_c','sink_C04226_c','sink_ins_c','sink_dimp_c','sink_C04546_c'}
sink_reaction_names1 = {'sink_uri_c','sink_udpgalur_c','sink_udpgal_c','sink_uacmam_c','sink_uacgamo_c','sink_rnam_c','sink_2dr1p_c','sink_mocogdp_c'}
sink_reaction_names2 = {'sink_3hpp_c','sink_gal14lac_c','sink_gsn_c','sink_nicrns_c','sink_xtsn_c','sink_udpxyl_c','sink_thym_c','sink_2hyoxplac_c'}
sink_reaction_names3 = {'sink_C05165_c','sink_C21723_c','sink_C21724_c','sink_cdp4dh6doglc_c','sink_cytd_c','sink_dad_2_c','sink_dcyt_c','sink_dgsn_c'}
sink_reaction_names4 = {'sink_dialurate_c','sink_dtdprmn_c','sink_gdpdrhmn_c','sink_gdpfuc_c','sink_gdpper_c','sink_hemeA_c','sink_hgentis_c','sink_mococdp_c'}
sink_reaction_names5 = {'sink_ppgpp_c','sink_rhcys_c','sink_sheme_c','sink_thmpp_c'}
sink_reaction_names = [sink_reaction_names0,sink_reaction_names1,sink_reaction_names2,sink_reaction_names3,sink_reaction_names4,sink_reaction_names5]


[~,NFind_sink_reactions]=ismember(strcat('NF_',sink_reaction_names),MILPstructure.varNames);


MILPstructure.lb(NFind_NIT) = -100
MILPstructure.ub(NFind_NIT) = 100
MILPstructure.ub(NFind_sink_reactions) = 1
MILPstructure.lb(NFind_sink_reactions) = -1
MILPstructure.ub(NFind_demand_reactions) = 1
MILPstructure.lb(NFind_demand_reactions) = -1
MILPstructure.ub(NFind_demand_h2o2_c) = 0
MILPstructure.lb(NFind_demand_h2o2_c) = 0
%MILPstructure.ub(NFind_EX_o2_e) = 100
%MILPstructure.lb(NFind_EX_o2_e) = -100
MILPstructure.ub(NFind_EX_succ_e) = 10
MILPstructure.lb(NFind_EX_succ_e) = -10
MILPstructure.ub(NFind_EX_h2o_e) = 200
MILPstructure.lb(NFind_EX_h2o_e) = -200
MILPstructure.ub(NFind_NH3c) = 200
MILPstructure.lb(NFind_NH3c) = -200



params.OutputFlag = 1;
params.DisplayInterval = 5;
%params.TimeLimit = cobraParams.timeLimit;
sol_Z_de = gurobi(MILPstructure, params);


repMILP = MILPstructure;
[~,num_vars] = size(repMILP.A);
repMILP.lb(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));
%repMILP.ub(find(repMILP.obj)) = sol_Z_de.x(find(repMILP.obj));
Qmat = sparse(num_vars,num_vars);
for i = 1:numel(NFind)
    Qmat(NFind(i), NFind(i)) = 1;
end

repMILP = rmfield(repMILP,'obj');
repMILP.Q = Qmat;
repMILP.modelsense = 'min';
sol_rep = gurobi(repMILP, params);


[MILP] = create1normSolutionMILP(MILPstructure, NFind, sol_Z_de.objval);
sol_rep2 = gurobi(MILP,params);
P_D2 = sol_rep2.x(NFind);
R_D = RxnExpr;
