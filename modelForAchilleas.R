##intermediate variables 
#==1/ - Bloc reel======================
#===1-1/ Production==================
VYDMNANet=YDMNANet/IPPRODMNA
VPRODMNA=(VYDMNAe+VINVNA)/(1-aPFATMNA)
VPRODA_variantes=VPRODA*(1-VARPROD_VIG)+VPRODA2*VARPROD_VIG
VPRODG=PRODG/IPPRODG

#===1-2/ PIB, valeur ajoutees==================
VPIB=omegaIPC*VCONS+omegaIPI*VITOT+omegaIPCG*VCG+omegaIPX*VXBS-omegaIPM*VMBS+omegaIPS_flux*VS_Flux
VVAMNA=VPRODMNA/omegaIPVANA-(omegaIPUIAMNA*VUIMNAA+omegaIPUINAMNA*(VUIMNAIA+VUIMNANAHIA))/omegaIPVANA
VVAAG=VPRODA_variantes/omegaIPVAAG-(omegaIPUIAA*VUIAA+omegaIPUINAA*(VUIAIA+VUIANAHIA))/omegaIPVAAG
VVAG=VPRODG/omegaIPVAG-(omegaIPUIAG*VUIGA+omegaIPUINAG*(VUIGIA+VUIGNAHIA))/omegaIPVAG
VVAT= omegaIPVAAGT*VVAAG+ omegaIPVAMNAT *VVAMNA+ omegaIPVAGT *VVAG

#===1-3/ Taxes sur les produits  nets des subventions======================================
#VTPTNS= VPIB/ omegaIPTNSPIB -omegaIPVAPIB*VVA/omegaIPTNSPIB

#===1-4/ Emploi et chômage===============================
#===1-4-1/ Emploi ===============================
L=LAG+LMNA+LG
LAG=VPRODA_variantes/(rhoAG+rhoAG7*DUM1719+rhoAG8*DUM09)
LMNA=VPRODMNA/rhoMNA
rhoMNA = rhoMNALT*(1-rhoMNA4*DUM08)*(1-rhoMNA6*DUM0910)
rhoMNA_Cr=(rhoMNA1+rhoMNA2*(VISOEI/VSKNA-TDEC-rhoMNA3))*(1-rhoMNA5*DUM0712)
LS=LSMNA+LSAG+LG
LSAG=tLSAG*LAG
LSMNA=tLSMNA*LMNA
LRNSMNA=tLRNSMNA*LMNA
LNRMNA=tLNRMNA*LMNA
#===1-4-2/ Chômage===============================
TCHOM=(1-L/PACT)*100

#===1-5/ Revenu des menages=================
VRDM=RDM/IPCF
VRDMA=RDMA/IPCF
VRDMNA=RDMNA/IPCF

#===1-6/ demande totale==========================
VDEMAeC=tetamba4*(VUIAA+VUIMNAA+VUIGA+VCGA+VCONSA+VIA+VXBA)

#===1-6-1/ utilisations intermediaires==================
VUIAA =aUIAA*VPRODA_variantes
VUIANA =aUIANA*VPRODA_variantes
VUIANAHIA=aUIANAHIA*VPRODA_variantes
VUIAIA=aUIAIA*VPRODA_variantes
VUINANA =aUINANA *VPRODMNA
VUINAA = aUINAA *VPRODMNA
VUIMNAA = aUIMNAA*VPRODMNA
VUIMNANAHIA=aUIMNANAHIA*VPRODMNA
VUIMNAIA=aUIMNAIA*VPRODMNA
VUINA = VUIANA + VUINANA
VUIGA=UIGA/IPUIA
VUIGIA=UIGIA/IPUINA
VUIGNAHIA=( UIGNAHIAHIC+UIG_SIF+UIG_SF )/IPUINA

#===1-6-2/consommation des menages===============
VCONS=VCONSA*omegaIPCA+VCONSIA*omegaIPCIA+VCONSNAHIA*omegaIPCNAHIA
VCONSA=(CONSAAut+epsilonA*(CONS-CONSAAut-CONSIAAut-CONSNAHIAAut))/IPCFA
VCONSNA=VCONSIA*omegaIPCIANA+VCONSNAHIA*omegaIPCNAHIANA
VCONSIA=(CONSIAAut+epsilonIA*(CONS-CONSAAut-CONSIAAut-CONSNAHIAAut))/IPCFIA
VCONSNAHIA=(CONS-VCONSIA*IPCFIA-VCONSA*IPCFA)/IPCFNAHIA
VCONSNAHIAHTS=VCONSNAHIA/((1+TTVANAHIAC0)*(1-TSPTNAHIA0)*(1+TTICNAHIA0)*(1+TAIPTNAHIA0))
VCONSIAHTS=VCONSIA/((1+TTVAIA0)*(1-TSPTIA0)*(1+TTICIA0)*(1+TAIPTIA0))
epsilonIA=(epsilon1IA-epsilon2IA)*(1/omegaIPCIA)^epsilon3IA
#===1-6-3/consommation des APU=========================
VCG=omegaIPCGM*VCGMNA+omegaIPCGM*VCGA+omegaIPCGG*VCGG
VCGA=CGA/IPCGM
VCGMNA=CGMNA/IPCGM
VCGG=CGG/IPCGG
#===1-6-4/ FBCF=====================================
VITOT=VISOEI+VILOG+VIG_variante
VISOEI=(VISOEILT)*(1+kappa11*DUM08)*(1+kappa12*DUM15)+VISOEI_LEAP
VISOEIC=VSKNA*(kappa1var+kappa2*(VVAMNA/VSKNA)+kappa3*(PROFITSOEIC/(VSKNA*IPI)))
kappa1var = 1/(1.0 + exp(-0.2*(t - 1)))*(kappa1 - kappa1*1.8) + kappa1*1.8

PROFITSOEIC=EBEMNAHIC-RNSMNA+IntSHIC_R-IntSHIC_V+HarmSHIC-RlyS_V-RassSHIC_V+RassSHIC_R
VILOG=VILOGLT
VILOGC=VRDMNA*(kappaVILOG1-kappaVILOG2*ratioDetteRevenuC)
# ratioDetteRevenuC=(kappaVILOG3+kappaVILOG4*TintD)*(1+kappaVILOG5*DUM08)*(1+kappaVILOG6*DUM10)
ratioDetteRevenuC=(kappaVILOG3+kappaVILOG4*CBIC_M/RDMNA)*(1+kappaVILOG5*DUM08)*(1+kappaVILOG6*DUM10)

VINA=VITOT/omegaIPINA-VIA*omegaIPIA/omegaIPINA 
VIA=(tVIA*VPRODA_variantes)*(1+tVIA1*DUM0714)
VIG_variante = VIG+varianteVIG*VPIB_observe
#=== 1-6-5/ Stocks======================================
VINVNA=upsilon2NA*(upsilon1NA*(VYDMNAe)-VSNA)
VS_Flux=S_Flux/IPS_Flux
VSA_Flux= SA_Flux/IPSA_Flux
VSNA_Flux= SNA_Flux/IPSNA_Flux
#===1-7/ commerce exterieur ============================
#===1-7-1/ import========================================
VMBS=omegaIPMBA*VMBA+omegaIPMHEA*VMBHEA+omegaIPME*VME+omegaIPMS*VMS
VMBA=sigmaMBA*VDEMAe
sigmaMBAC=(tetamba0+tetamba1*(VDEMAe-VPRODA_variantes)/(VPRODA_variantes*(1+tetamba2*IPMBA/IPPRODA)^tetamba3))
VMBHEA=(sigmaMC*VCONS+sigmaMI*(VITOT-VIA)+sigmaMX*VXBHAOCP+sigmaMUI*VUINA)
sigmaMCC=1/(1+iota1*IPMHEA/IPCFNAHIA)^iota2*((1+iota9*DUM1619)*(1+iota10*DUM0809))
# ((1/(resBaseline$sigmaMC[1]*(1+resMin['iota9'])))^(1/resMin['iota2'])-1)/(resBaseline$IPMHEA[1]/resBaseline$IPCFNAHIA[1])=iota1
sigmaMIC=1/(1+iota3*IPMHEA/IPI)^(iota4*(1+iota11*DUM1619))
# ((1/(resBaseline$sigmaMI[1]))^(1/(resMin['iota4']*(1+resMin['iota11'])))-1)/(resBaseline$IPMHEA[1]/resBaseline$IPI[1])
sigmaMXC=1/(1+iota5*IPMHEA/IPXHAOCP)^iota6 #sigmaMXC=1/(1+iota5*(1+AUGM1950*parsX)*IPMHEA/IPXHAOCP)^iota6
# ((1/(resBaseline$sigmaMX[1]))^(1/resMin['iota6'])-1)/(resBaseline$IPMHEA[1]/resBaseline$IPXHAOCP[1])
sigmaMUIC=1/(1+iota7*IPMHEA/IPUINA)^(iota8*(1+iota12*DUM08))
# ((1/(resBaseline$sigmaMUI[1]))^(1/resMin['iota8'])-1)/(resBaseline$IPMHEA[1]/resBaseline$IPUINA[1])
VME=(sigmaME*(VCONS+VITOT-VIA+VXBS+VCGMNA)^iotaME1)*(1+iotaME4*DUM08)*(1+iotaME5*DUM1219)*(1+iotaME6*DUM1214)
sigmaMEC=1/(1+iotaME2*IPME/IPPRODMNA)^iotaME3
VMS=(sigmaMS1*VMBHEA+sigmaMS2*VME+sigmaMS3*VMBA)*(1+sigmaMS4*DUM0708)*(1+sigmaMS5*DUM19)
#===1-7-2/ export========================================
VXBS=omegaIPXBA*VXBA+omegaIPXHAOCP*VXBHAOCP+omegaIPXOCP*VXOCP+omegaIPXSHVTR*VXSHVTR+omegaIPXSTR*VXSTR+omegaIPXSV*VXSV
# VXBAC=tetamba6*VDMHOCP*VPRODA_variantes^tetamba7
VXBAC=1.03*VXBA #TODOANTOINE
VXBHAOCP=(1+nuXBHAOCP5*DUM09)*sigmaX*(nuXBHAOCP2*VDMHOCP)^nuXBHAOCP3
sigmaXC=(1+nuXBHAOCP4Temp*19)*nuXBHAOCP0*(IPCET/(IPPRODMNA*TCEN))^nuXBHAOCP1
nuXBHAOCP4Temp=((1/(1.0 + exp(-elastExp*(t - initExp))))*(UBExp - LBExp) + LBExp)
VXSTRC = nuVXSTR0*VXBHAOCP + nuVXSTR1*ARRTOUR
# VXSHVTRC = nuVXSHVTR0*VXBHAOCP + nuVXSHVTR1*VXSTR
VXSHVTRC=1.03*VXSHVTR #TODOANTOINE
VXSV=XSV/IPXSV
#======================================================
#======================================================
#==2/- Boucle prix salaires======================
#===2-1/Taux de remuneration salarial et non salarial==================
TSAG_BC=VTSAG_B*IPCF
TSMNA= (1+MuMNA)*TSMNA_B
TSMNA_B=VTSMNA_B*IPCFNAHIALT*(1+omega3*DUM0815)*(1+omega4*DUM1415)*(1+omega5*DUM12)
TNSMNA=tTNSMNA*TSMNA
SMIG=VSMIG*IPCF
#===2-2/Prix à la production==========
IPPRODMNALT_Infl=betaIPVANA*(IPPRODMNAC-IPPRODMNA)
IPPRODMNA=IPPRODMNALT+phiIPVANA0*DUM08+phiIPVANA1*DUM10+phiIPVANA2*DUM15+phiIPVANA3*DUM17
IPPRODMNAC=(1+muVAMNA)*(TSMNA/(muVANA2*rhoMNA)+(IPUIA*(1+muVANA3*DUM13)*(1+muVANA4*DUM0809)*(1+muVANA5*DUM1516))*aUIMNAA+IPUINA*aUIMNAIA+IPUINA*aUIMNANAHIA)
muVAMNAC=muVANA0-muVANA1*VSNA/VPRODMNA
IPPRODG=IPCGG
#===2-3/Prix du PIB et prix de  valeurs  ajoutees ==================
IPPIB=PIB/VPIB
IPVAMNA=VAMNA/VVAMNA
IPVAAG=VAAG/VVAAG
IPVAG=VAG/VVAG
IPVAT= VAT/VVAT

#=== 2-4/ Prix de la demande ==========================
#===2-4-1/ Prix utilisations intermediaires==================
IPUIA= IPUIALT
IPUINA= IPUINALT
IPUIANA= IPUIANALT*(1+phiIPUIANA4*DUM0911)
IPUINANA= IPUINANALT

#===2-4-2/Prix consommation des menages===============
#infl=(IPCF-IPCFLag)/IPCFLag*iotaIPCFLag0
IPCF=CONS/VCONS
IPCFNA=CONSNA/VCONSNA
# IPCFA = (phiIPCFA4+phiIPCFA5*IPPRODA+phiIPCFA6*IPMBA+phiIPCFA7*IPME+phiIPCFA8*(1-phiIPCFA5-phiIPCFA6-phiIPCFA7)*IPCFALT)*(1+phiIPCFA9*DUM0910)
IPCFA = IPCFALT
# IPCFAC = phiIPCFA0+phiIPCFA1*IPPRODA+phiIPCFA2*IPMBA+phiIPCFA3*(1-phiIPCFA1-phiIPCFA2)*IPME
IPCFAC = 1.015*IPCFA #TODOANTOINE
IPCFNAHIA=IPCFNAHIALT
IPCFNAHIAHTS=IPCFNAHIA*(1+TTVANAHIAC0)*(1-TSPTNAHIA0)*(1+TAIPTNAHIA0)*(1+TTICNAHIA0)/((1+TTVANAHIAC)*(1-TSPTNAHIA)*(1+TAIPTNAHIA))-TTICNAHIA
IPCFIA=IPCFIALT
IPCFIAHTS=IPCFIA*(1+TTVAIA0)*(1-TSPTIA0)*(1+TAIPTIA0)*(1+TTICIA0)/((1+TTVAIA)*(1-TSPTIA)*(1+TAIPTIA))-TTICIA

#===2-4-3/ Prix consommation des APU=========================
IPCG=CG/VCG
IPCGM= IPCGMLT
# IPCGGC=muG0+muG1*TSG_B
IPCGGC = 1.015*IPCGG #TODOANTOINE

#===2-4-4/ Prix FBCF========================================
IPI=IPILT
# inflIPI=(IPI-IPILag)/IPILag*iotaIPILag0
IPIA = (phiIPIA4+phiIPIA5*IPPRODA+phiIPIA6*IPMBA+phiIPIA7*IPME+phiIPIA8*(1-phiIPIA5-phiIPIA6-phiIPIA7)*IPIALT)*(1+phiIPIA9*DUM11)*(1+phiIPIA10*DUM14)
# IPIAC = phiIPIA0+phiIPIA1*IPPRODA+phiIPIA2*IPMBA+phiIPIA3*(1-phiIPIA1-phiIPIA2)*IPME
IPIAC = 1.015*IPIA#TODOANTOINE
IPINA=INA/VINA

#===2-4-5/ Prix commerce exterieur et prix etrangers et taux de change===================
#===2-4-5-1/ Prix import========================================
IPM=MBS/VMBS
IPMBA=IPMBALT+alphaIPMBA2*DUM0809+alphaIPMBA3*DUM13+alphaIPMBA4*DUM1112+alphaIPMBA5*DUM10
# IPMBAC =alphaIPMBA0+alphaIPMBA1*(Pble_dolr*TCDOLR/(Pble_dolr_2007*TCDOLR_2007))
IPMBAC=1.015*IPMBA#TODOANTOINE
IPME=IPMELT+alphaIPME2*DUM08+alphaIPME3*DUM16
# IPMEC =alphaIPME0+alphaIPME1*(PPET*TCDOLR/(PPET_2007*TCDOLR_2007))
IPMEC=1.015*IPME#TODOANTOINE
# resBaseline$IPME[1]-resMin['alphaIPME1']*(resBaseline$PPET[1]*resBaseline$TCDOLR[1]/(resBaseline$PPET_2007[1]*resBaseline$TCDOLR_2007[1]))
IPMHEA=IPMHEALT+alphaIPMHEA1*DUM1516+alphaIPMHEA2*DUM08
# IPMHEAC=alphaIPMHEA0*(IPCET/TCEN)
IPMHEAC=1.015*IPMHEA#TODOANTOINE
IPMS=IPMSLT
# IPMSC= alphaIPMS*(IPCET/TCEN)*(1+alphaIPMS1*DUM0710)*(1+alphaIPMS2*DUM1417)
IPMSC=1.015*IPMS#TODOANTOINE
IPME_Infl=betaIPME*(IPMEC-IPMELT)
IPMBA_Infl=betaIPMBA*(IPMBAC-IPMBALT)
IPMHEA_Infl=betaIPMHEA*(IPMHEAC-IPMHEALT)

#===2-4-5-2/ Prix export========================================
IPX=XBS/VXBS
# IPXBA = (phiIPXBA4+phiIPXBA5*IPPRODA+phiIPXBA6*IPMBA+phiIPXBA7*IPME+phiIPXBA8*(1-phiIPXBA5-phiIPXBA6-phiIPXBA7)*IPXBALT)*(1+phiIPXBA9*DUM1112)*(1+phiIPXBA10*DUM19)
# IPXBAC = phiIPXBA0+phiIPXBA1*IPPRODA+phiIPXBA2*IPMBA+phiIPXBA3*(1-phiIPXBA1-phiIPXBA2)*IPME
IPXBAC = 1.015*IPXBA #TODOANTOINE
IPXOCP= alphaIPXOCP1* (PXPH_dolr *TCDOLR/(PXPH_dolr_2007*TCDOLR_2007))+ alphaIPXOCP2*(PXENG_dolr*TCDOLR/(PXENG_dolr_2007* TCDOLR_2007))+ alphaIPXOCP3* (PXAPH_dolr*TCDOLR/( PXAPH_dolr_2007*TCDOLR_2007))
IPXHAOCPC=1.015*IPXHAOCP #TODOANTOINE
# IPXHAOCPC=alphaIPXHAOCP1*IPPRODMNA+ alphaIPXHAOCP2tr*(IPCET/TCEN)
# alphaIPXHAOCP2tr=1/(1.0 + exp(0.5*(t - 5)))*(alphaIPXHAOCP2*0.8 - alphaIPXHAOCP2*1.27) + alphaIPXHAOCP2*1.27
IPXSHVTRC=1.015*IPXSHVTR #TODOANTOINE
# IPXSHVTRC=(alphaIPXSHVTR1*IPPRODMNA+ alphaIPXSHVTR2tr*(IPCET/TCEN))*(1+alphaIPXSHVTR3*DUM0710)*(1+alphaIPXSHVTR4*DUM1314)
# alphaIPXSHVTR2tr=1/(1.0 + exp(0.5*(t - 5)))*(alphaIPXSHVTR2*0.2 - alphaIPXSHVTR2*2.4) + alphaIPXSHVTR2*2.4
IPXSTRC=1.015*IPXSTR #TODOANTOINE
# IPXSTRC =(alphaIPXSTR0*IPPRODMNA+ alphaIPXSTR1tr*(IPCET/TCEN))*(1+alphaIPXSTR2*DUM1114)
# alphaIPXSTR1tr=1/(1.0 + exp(0.5*(t - 5)))*(alphaIPXSTR1*0.5 - alphaIPXSTR1*2.86) + alphaIPXSTR1*2.86
IPXSVC=1.015*IPXSV #TODOANTOINE
# IPXSVC=(alphaIPXSV0*IPPRODMNA+ alphaIPXSV1tr*(IPCET/TCEN))*(1+alphaIPXSV2*DUM1215)
# alphaIPXSV1tr=1/(1.0 + exp(0.5*(t - 5)))*(alphaIPXSV1*0.5 - alphaIPXSV1*2) + alphaIPXSV1*2
#=== 2-4-5-3/Prix etrangers et taux de change============================
IPCETC= (alphaIPCET0+alphaIPCET1*IPCUE+ alphaIPCET2*IPCUSA)
TCENC=(alphaTCEN0+alphaTCEN1*(TCEURO_2007/TCEURO)+alphaTCEN2*(TCDOLR_2007/TCDOLR))*(1+alphaTCEN3*DUM1619)*(1-varianteTCEN)

#===2-5/ Taux d’interets========================
#CMRBQC=VCMRBQ+iotaIPCFLag1*inflLT+(1-iotaIPCFLag1)*infl

#======================================================
#======================================================
#====3/- Bloc nominal======================
#===3-1/ Production et produit ==================
YDMNAe=VYDMNAe*IPPRODMNA
YDMNANet=DEMNA-MTCNA-TAXESNA+SPTNA-MNA 
PRODMNA=IPPRODMNA*VPRODMNA
PG=CGG
PRODA_variantes=PRODA*(1-VARPROD_VIG)+PRODA2*VARPROD_VIG
PRODG=PG+PFATG
PRODMNAHIC=PRODMNA-PRODIC
PRODIC=(SIFIM+SF_IC)-(MTCIC+TAXESIC+M_SIF)+PFATIC
PRODBAM=PRODBAM_NV+PFATBAM
PFATA=aPFATA*PRODA_variantes
#PFATMNA=aPFATMNA*PRODMNA
PFATG=aPFATG*PG
PFATMNA=-PFATG-PFATA
PFATIC=aPFATIC*((SIFIM+SF_IC)-(MTCIC+TAXESIC+M_SIF))
PFATMNAHIC=PFATMNA-PFATIC-PFATBAM

#===3-2/ PIB, valeur ajoutees==================
PIB = CONS+INA+IA+CG+XBS-MBS+S_Flux
VAT=VAAG+VAMNA+VAG
VAAG=PRODA_variantes-UIAA-UIAIA-UIANAHIAHIC-UIA_SIF-UIA_SF
VAMNA=VAMNAHIC+VAIC+VABAM
VAMNAHIC=(UIIA+UINAHIAHIC+CONSIA+CONSNAHIAHIC+CGMNAHIC+INA+SNA_Flux+XNAHIC)-(MTCNAHIC+TAXESNAHIC-SPTNA)+PFATMNAHIC-MNAHIC-(UIMNAA+UIMNAIA+UIMNANAHIAHIC+UIMNAHIC_SIF+UIMNAHIC_SF)
VAIC=PRODIC-UIICNAHIA-PRODBAM_NV
VABAM=PRODBAM-UIBAMNAHIA
VAG=PRODG-UIGA-UIGIA-UIGNAHIAHIC-UIG_SIF-UIG_SF
#===3-3/ Marges commerce et transport======================================
MTCA=-MTCNA
MTCNA= aMTCNA*DEMNA
MTCNAHIC =MTCNA-MTCIC
MTCIC=aMTCIC*(SIFIM+SF_IC)
#===3-4/ Taxes sur les produits et subventions sur produits=============================
TPTNS=TAXES-SPT
TAXES=TAXESA+TAXESNA
TAXESA=DDA
TAXESNA = DDNA + TICNA + TVANA + AIPTNA
TAXESNAHIC=TAXESNA-TAXESIC
TAXESIC=TVAIC+AIPTIC

#===3-4-1/ Droits de douanes ===============================
DDA=TDDA*MBA
DDNA=DDIA+DDNAHIA
DDIA=TDDIA*MBIA
DDNAHIA=TDDNAHIA*(MBHEAIA+ME+MS)
#===3-4-2/ TIC =========================================
TICNA=TICIA+TICNAHIA
TICIA=TTICIA*VCONSIAHTS
TICNAHIA=TTICNAHIA*VCONSNAHIAHTS
#===3-4-3/ TVA =========================================
TVANA=TVAIA+TVANAHIAC+TVANAF
TVAIA = TTVAIA*CONSIA/(1+TTVAIA)
TVANAHIAC = TTVANAHIAC*CONSNAHIA/(1+TTVANAHIAC)
TVANAF = TTVANAF*INA/(1+TTVANAF)
TVAIC=TTVAIC/(1+TTVAIC)*(CONS_SIF+CONS_SF)
#===3-4-4/ Autres impots ==================================
AIPTNA=AIPTIA+AIPTNAHIA
AIPTIA = TAIPTIA*CONSIA/((1+TAIPTIA)*(1-TSPTIA)*(1+TTVAIA))
AIPTNAHIA = TAIPTNAHIA*CONSNAHIA/((1+TAIPTNAHIA)*(1-TSPTNAHIA)*(1+TTVANAHIAC))
AIPTIC=TAIPTIC/(1+TAIPTIC)*(CONS_SIF+CONS_SF)
#===3-4-5/  subventions sur produits======================
SPT=SPTA+SPTNA
SPTA = TSPTA*CONSA/(1-TSPTA)
SPTNA=SPTIA+SPTNAHIA
SPTIA = TSPTIA*CONSIA/((1-TSPTIA)*(1+TTVAIA))
SPTNAHIA = TSPTNAHIA*CONSNAHIA/((1-TSPTNAHIA)*(1+TTVANAHIAC))

#====3-5/ Partage de valeur ajoutee================
# ===3-5-1/ Remuneration salariés ===================
# ===3-5-1-1/ Remuneration salaries brut ===================
RSTOT_B=RSMNA_B+RSA_B+RSG_B
RSA_B= TSAG_B*LSAG
RSMNA_B= TSMNA_B*LSMNA
RSG_B= TSG_B*LG
RSMNAHIC_B= RSMNA_B-RSIC_B-RSBAM_B
RSIC_B=aRSIC*RSMNA_B-RSBAM_B
#===3-5-1-2/ Cotisations employeurs===================
CSETOT=CSEA+CSEMNA+CSEG
CSEA = MuA*RSA_B 
CSEMNA= MuMNA*RSMNA_B
CSEMNAHIC= CSEMNA-CSEIC-CSEBAM
CSEIC=MuIC*RSIC_B
CSEBAM=MuIC*RSBAM_B
CSEG = MuG*RSG_B

#===3-5-1-3/ Cotisations employes=========================
CSSM_V= TauH*(RSTOT_B+RNSMNA)
CSSG_R=tCSSG*CSSM_V
CSSW_R=tCSSW*CSSM_V
CSSSHIC_R=tCSSSHIC*CSSM_V
CSSIC=tCSSIC*CSSM_V

#===3-5-2/ EBE et revenue mixte =============================
EBEMNA=EBEMNAHIC+EBEIC+EBEBAM
EBEMNAHIC=VAMNAHIC-RSMNAHIC_B-CSEMNAHIC-IPDMNAHIC+SPDMNAHIC
EBEIC=VAIC-RSIC_B-CSEIC-IPDIC+SPDIC
EBEBAM=VABAM-RSBAM_B-CSEBAM-IPDBAM
EBEG=VAG-RSG_B-IPDG-CSEG
EBEMG=RNSA+RNSMNA
RNSA=VAAG-RSA_B-CSEA-IPDA

# RNSMNA=TNSMNA*LRNSMNA+0.2*max(0,0.1*VSKNA*IPI-LNDSHIC_N)
RNSMNA=0.15*PRODMNAHIC#TODOANTOINE

#===3-6/ Taxes  et subventions sur la production===================
#===3-6-1/ Taxes sur la production===================
IPD=IPDA+IPDMNAHIC+IPDIC+IPDBAM+IPDG
IPDA = TIPDA/(1+TIPDA)*PRODA_variantes
IPDMNAHIC = TIPDMNAHIC/(1+TIPDMNAHIC)*PRODMNAHIC
IPDIC = TIPDIC/(1+TIPDIC)*PRODIC
IPDG = TIPDG/(1+TIPDG)*PRODG

#=== 3-6-2/ Subventions sur la production=========================
SPD=SPDMNAHIC+SPDIC
SPDMNAHIC = TSPDMNAHIC*PRODMNAHIC/((1+TIPDMNAHIC)*(1-TSPDMNAHIC))
SPDIC = TSPDIC*PRODIC/((1+TIPDIC)*(1-TSPDIC))

#===3-7/Interets=================================
#===APU
IntG_V=IntDg_N+IntCBIC_G+IntCw_G+IntBg_N+IntBgW
IntDg_N=TintB*Dg_N/100
IntCBIC_G=TintB*CBIC_G/100
IntCw_G=TintRfw_cd*Cw_G/100
IntBg_N=TintBg*Bg_N/100
IntBgW=TintBgW*BgW/100

IntG_R=IntDBG_IC
IntDBG_IC=TintB*DBG_IC/100

#===Menages
IntM_V = IntCBIC_M
IntCBIC_M= TintB*CBIC_M/100

IntM_R= IntDBM_IC
IntDBM_IC = TintB*DBM_IC/100
#===SHIC
IntSHIC_V=IntLNDSHIC_N+IntCwIC_SHIC+IntCwW_SHIC+IntBshicW
IntLNDSHIC_N =max(0,TintB* LNDSHIC_N/100)
IntCwIC_SHIC= TintB *CwIC_SHIC/100
IntCwW_SHIC= TintRfw_cd *CwW_SHIC/100
IntBshicW=TintBshicW*BshicW/100

IntSHIC_R = IntDBSHIC_IC + IntBgSHIC + IntDgSHIC+IntDwSHIC_IC
IntDBSHIC_IC = TintB * DBSHIC_IC/100
IntBgSHIC=TintBg*BgSHIC/100
IntDgSHIC= TintB*DgSHIC/100
IntDwSHIC_IC= TintB *DwSHIC_IC/100

#=====BAM
IntBAM_V= IntDBIC_BAM
IntDBIC_BAM = TintB* DBIC_BAM/100

IntBAM_R= IntDwBAM_W+IntBgBAM+IntBwBAM+IntCBBAM_IC+IntORDTS
IntDwBAM_W = TintAOR *DwBAM_W/100
IntBgBAM = TintBg*BgBAM/100
IntBwBAM = TintAOR*BwBAM/100
#IntCBBAM_IC=TintB*CBBAM_IC/100
IntORDTS = TintAOR *ORDTS/100
#===IC (intermediation financiere)
IntIC_V=IntDB_IC_N+IntDw_IC_N+IntCBBAM_IC+IntCwW_IC
IntDB_IC_N=TintB*DB_IC_N/100
IntDw_IC_N= TintB *Dw_IC_N/100
IntCBBAM_IC=TintB*CBBAM_IC/100
IntCwW_IC= TintRfw_cd *CwW_IC/100

IntIC_R= IntDBIC_BAM +IntDgIC+intDwIC_W+ IntLNDIC_N +IntCwIC_N+IntBgIC +IntBwIC
IntDgIC=TintB*DgIC/100
intDwIC_W= TintAOR *DwIC_W/100
IntLNDIC_N =TintB* LNDIC_N/100
IntCwIC_N= TintB *CwIC_N/100
IntBgIC=TintBg*BgIC/100
IntBwIC = TintAOR*BwIC/100
#===RDM
IntW_V = IntDw_W+IntBw+IntCwIC_W+IntORDTS
IntDw_W= TintAOR *Dw_W/100
IntBw = TintAOR*Bw/100
IntCwIC_W = TintB *CwIC_W/100

IntW_R=IntDwW_IC+IntCwW+IntBshicgW+IntBfW_SHIC
IntDwW_IC= TintB *DwW_IC/100
IntCwW= TintRfw_cd *CwW/100
IntBfW_SHIC=TintB*BfW_SHIC/100
IntBshicgW=IntBshicW+IntBgW
#===3-8/Dividendes====================================
DivS=DivSHIC+DivIC+DivBAM

DivIC=tDivIC*max(0,(EBEIC+IntIC_R-IntIC_V+CSSIC-PSIC-TRFIC-RassIC_V)-DivBAM)

DivSHIC=0.023*PRODMNAHIC#TODOANTOINE
										
# DivSHIC=tDivSHICtr*(PRODMNAHIC-(UIMNAA+UIMNAIA+UIMNANAHIAHIC+UIMNAHIC_SIF+UIMNAHIC_SF)-RSMNAHIC_B-CSEMNAHIC-IPDMNAHIC+SPDMNAHIC+IntSHIC_R-IntSHIC_V+CSSSHIC_R-PSSHIC-TRFSHIC-RlyS_V-RassSHIC_V+BrIDE)
# 
# tDivSHICtr=0.06335131-0.08+0.4*(PRODMNAHIC-(UIMNAA+UIMNAIA+UIMNANAHIAHIC+UIMNAHIC_SIF+UIMNAHIC_SF))/VSKNA*IPI

DivM_R = tDivM*DivS
DivG = tDivG*DivS
DivW = tDivW*DivS

#===3-9/Loyers====================================
Rly_V=RlyS_V+RlyM_V+RlyW_V
RlyS_V=aRlyS_V*PRODMNAHIC
RlyM_V=aRlyM_V*PRODMNA
RlyW_V=aRlyW_V*PRODMNA
RlyM=RlyM_R-RlyM_V
RlyM_R=aRlyM_R*Rly_V
RlyG_R=aRlyG_V*Rly_V

#===3-10/Revenus attribues aux assures====================================
RassSHIC_V=aRassSHIC_V*PRODMNAHIC
RassSHIC_R=aRassSHIC_R*PRODMNAHIC
RassIC_R=aRassIC_R*RassSHIC_V
RassM_R=aRassM*RassSHIC_V
RassG_R=aRassG*RassSHIC_V
RassW_R=aRassW*RassSHIC_V

#===3-11/Bénéfices réinvestis IDE====================================
BrIDE=BrIDE_R-BrIDE_V
BrIDE_R=aBrIDE_R*S_IDE_S
BrIDE_V=aBrIDE_V*S_IDE_E

#=== 3-12/ Prestations===============================
PSM=TPSM*(LS*1000+RETR)/1000000
PSSHIC=PSSHIC_STR*PSM/100
PSIC=PSIC_STR*PSM/100
PSG= PSG_STR*PSM/100
PSW= PSW_STR*PSM/100

#===3-13/ Revenu des menages=================
RDM=RDMNA+RDMA
RDMNA=(RSMNA_B+RSG_B+RNSMNA+IntM_R-IntM_V+HarmM+DivM_R+TRF_MRE+TRFM_A+PSM+RlyM+RassM_R-CSSM_V)*(1-TTRPMNA )
RDMA=VAAG-CSEA-IPDA
#===3-14/ demande totale hors variation des stocks==========================
YDA=CGA+XBA+IA+CONSA+UIA
DEMNA= UINA+CONSNA+CGMNA+INA+XNA
#=== 3-14-1/ utilisations intermediaires==================
#=== UI en produit agricole
UIA=UIAA+UIMNAA+UIGA
UIAA =IPUIA* VUIAA
UIMNAA = IPUIA*VUIMNAA
UIGA=aUIGA*PRODG
#=== UI en produit non agricole
UINA=UIANA+UINANA
UIANA =IPUIANA* VUIANA
UINANA =IPUINANA* VUINANA
#=== UI en produit alimentaire
UIIA=UIAIA+UIMNAIA+UIGIA
UIAIA=IPUINA*VUIAIA
UIMNAIA=IPUINA*VUIMNAIA
UIGIA=aUIGIA*PRODG 
#=== UI en produit non agri hors alimentaire et hors IC 
UINAHIAHIC=UIANAHIAHIC+UIMNANAHIAHIC+UIICNAHIA+UIBAMNAHIA+UIGNAHIAHIC
UIANAHIAHIC=UIANAHIA-UIA_SIF-UIA_SF
UIANAHIA=IPUINA*VUIANAHIA
UIMNANAHIAHIC=UIMNANAHIA-UIICNAHIA-UIBAMNAHIA-UIMNAHIC_SIF-UIMNAHIC_SF
UIMNANAHIA=IPUINA*VUIMNANAHIA
UIICNAHIAHIC=aUIICNAHIAHIC*PRODIC
UIGNAHIAHIC=aUIGNAHIAHIC*PRODG
#=== UI en SIFIM IC
UI_SIF=UIA_SIF+UIMNAHIC_SIF+UIG_SIF
UIA_SIF=pA_SIF*UISHIC_SIF
UIMNAHIC_SIF=UISHIC_SIF-UIA_SIF
UISHIC_SIF= SIFIMCBSHIC+SIFIMDBSHIC
# UIIC_SIF= SIFIMCBIC+SIFIMDBIC
UIG_SIF=pCIG_SIF*SIFIMG
#=== UI en services financiers IC
UI_SF= UIA_SF+UIMNAHIC_SF+UIG_SF
UIA_SF=aUIA_SF*UIA_SIF
UIG_SF=aUIG_SF*UIG_SIF
UIMNAHIC_SF=aUIMNAHIC_SF*UIMNAHIC_SIF

#====3-14-2/ SIFIM et commission sur services financiers========================================
SIFIM= UI_SIF+CONS_SIF+CG_SIF+X_SIF
SF_IC=UI_SF+CONS_SF+CG_SF+X_SF
SIFIMG= SIFIMCBG+SIFIMDBG
SIFIMCBM=(TintD-TintB)*CBIC_M/100
SIFIMDBM=(TintB-TintC)*DBM_IC/100
SIFIMCBSHIC=(TintD-TintB)*(max(0,LNDSHIC_N)+CwIC_SHIC)/100
SIFIMDBSHIC=(TintB-TintC)*(DBSHIC_IC+DwSHIC_IC)/100
# SIFIMCBIC=(TintD-TintB)*CBIC_IC/100
# SIFIMDBIC=(TintB-TintC)*DBIC_IC/100
SIFIMCBG=(TintD-TintB)*CBIC_G/100
SIFIMDBG=(TintB-TintC)*DBG_IC/100
SIFIMCBW=(TintD-TintB)*CwIC_W/100
SIFIMDBW=0

#===3-14-3/ consommation des menages===============
#CONS=CONSA+CONSNA
CONSA=VCONSA*IPCFA
CONSNA=CONSNAHIA+CONSIA
CONSNAHIA=VCONSNAHIA*IPCFNAHIA
CONSIA=VCONSIA*IPCFIA
CONSNAHIAHIC= CONSNAHIA-CONS_SIF-CONS_SF
CONS_SIF= SIFIMCBM+SIFIMDBM
CONS_SF=aCONS_SF*CONS_SIF

CONST=(alpha1*RDMNA+alpha2*RDMA)*(1+alpha8*DUM08+alpha9*DUM09)
alpha1=alpha10-alpha11*RDMNA/POPU
alpha2=alpha20-alpha21*RDMA/POPR
epsilonA=max(0,(epsilon1A-epsilon2A*(RDMNA/POP))*(1/omegaIPCA)^epsilon3A)
CONSAAut=consamin*POP*IPCFA+alpha3*PRODA
CONSNAHIAAut=consnahiamin*POP*IPCFNAHIA
CONSIAAut=(consiamin+DUM10*alpha5)*POP*IPCFIA+alpha4*PRODA

#===3-14-4/ consommation des APU=========================
CG= CGA+CGNA
CGNA= CGMNA+CGG
CGMNAHIC=CGMNA-CG_SIF-CG_SF
CG_SIF=(1-pCIG_SIF)*SIFIMG
CG_SF=aCG_SF*CG_SIF
#===3-14-5/ FBCF=====================================
ITOT=ISOEI+ILOG+IG
ISOEI=VISOEI*IPI
ILOG=VILOG*IPI
INA=ILOG+IG+ISOEINA
ISOEINA=ISOEI-IA
IA=VIA*IPIA
IG=IPI*VIG_variante
#===3-14-6/ Stocks===================================
S_Flux= SA_Flux+ SNA_Flux 
SNA_Flux=PRODMNA+MTCNA+TAXESNA-SPTNA-PFATMNA+MNA-DEMNA
# SNA_Flux=PRODMNA+(MTCNAHIC+TAXESNAHIC-SPTNA)-PFATMNAHIC+MNAHIC-DEMNA
# SSNF_Flux=aSSNF_Flux*SNA_Flux
SA_Flux=PRODA_variantes-MTCNA+TAXESA-SPTA+MBA-YDA-PFATA
# SM_Flux=1997#aSM_Flux*SA_Flux
#===3-15/ commerce exterieur ============================
BCRDM=MBS-XBS
#===3-15-1/ import========================================
MBS=MBA+MNA
MBA=IPMBA*VMBA
MNA = MBHEA+ME+MS
MBHEA=VMBHEA*IPMHEA
ME = VME *IPME
MS = VMS *IPMS
MBIA=MBHEA*tMBIA
MBHEAIA=MBHEA*(1-tMBIA)
MNAHIC = MNA -M_SIF
M_SIF= ((TintCDgW-TintRfw_cd)*Cw_G +(TintCDshicW-TintRfw_cd)*(CwW_SHIC+CwW_IC)+(TintB-TintC)*DwW_IC)/100
#===3-15-2/ Export========================================
XBS=XBA+XNA
XBA=IPXBA*VXBA
XNA =XBHAOCP+XOCP+XSHVTR+XSTR+XSV
XBHAOCP = VXBHAOCP*IPXHAOCP
XOCP = VXOCP*IPXOCP
XSV =XSVpartouriste * ARRTOUR/1000000
XSHVTR =VXSHVTR * IPXSHVTR
XSTR =VXSTR * IPXSTR
XBSHIC=XBS-X_SIF-X_SF
XNAHIC= XNA-X_SIF-X_SF
X_SIF=SIFIMCBW+SIFIMDBW
X_SF=aX_SF*X_SIF
#======================================================
#===3-16/  Balance des paiements====================
#===Compte Courant
SCC = XBS - MBS + EXM_BPCN + SRPM + SRSE
SRSE=RSE_R-RSE_V
RSE_R = TRF_MRE + RSE_RHMRE
SCCPIB=SCC/PIB*100
#===3-17/Autres flux et Capacite de financement=================================
#===3-17-1/  Autres transferts courants nets hors MRE et Autres flux============================
TRFM_A=aTRFM_A*(TRFSHIC+TRFIC+TRFBAM+TRFW_A)
TRFG=(1-aTRFM_A)*(TRFSHIC+TRFIC+TRFBAM+TRFW_A)
TRFSHIC=tTRFSHIC*PRODMNAHIC
TRFIC=tTRFIC*PRODIC 

AutreSHIC=tAutreSHIC*PRODMNAHIC
AutreIC=tAutreIC*PRODIC 
AutreM=tAutreM*PIB
AutreW=tAutreW*PIB
AutreG=-AutreSHIC-AutreIC-AutreM-AutreW
#===3-17-2/ Revenu disponible brut (RD), Epargne et Capcité de Financement Capacite de financement=================================
#===  la différence entre le revenu disponible brut des comptes natioanux et celui calculé ici est les autres transferts courants (hors transferts MRE pour les ménages)==
#=== On suppose que les sociétés financières ne font pas de FBCF (car faible) et que la FBCF des EI est faite par le MNA===
RDSOEIHIC=(EBEMNAHIC-RNSMNA+IntSHIC_R-IntSHIC_V+HarmSHIC-DivSHIC-RlyS_V-RassSHIC_V+RassSHIC_R+BrIDE+CSSSHIC_R-PSSHIC-TRFSHIC)*(1-TTRPSHIC)
EPARSOEIHIC= RDSOEIHIC
CAPFSOEIHIC= EPARSOEIHIC +AutreSHIC-ISOEINA-S_Flux
RDIC=(EBEIC+IntIC_R-IntIC_V+HarmIC-DivIC+RassIC_R+CSSIC-PSIC-TRFIC)*(1-TTRPIC)
EPARIC= RDIC
CAPFIC= EPARIC +AutreIC
RDBAM=(EBEBAM+IntBAM_R-IntBAM_V+HarmBAM)*(1-TTRPBAM)-DivBAM-TRFBAM
EPARBAM=RDBAM
CAPFBAM=EPARBAM
EPARM= RDM- CONS
CAPFM= EPARM +AutreM-ILOG-IA
RDG=EBEG+TAXES-SPT+IPD-SPD+IntG_R-IntG_V+HarmG+DivG+RlyG_R+RassG_R+TRPS+TRPMNA+CSSG_R+ CSETOT -PSG+TRFG
EPARG= RDG-CG
CAPFG= EPARG +AutreG -IG
EPARW=BCRDM+IntW_R-IntW_V+HarmW+DivW-RlyW_V+RassW_R-BrIDE-TRF_MRE+CSSW_R-PSW-TRFW_A
CAPFW=EPARW+AutreW


#===3-17-2/ Harmonisation===========================================
HarmSHIC =-HarmIC-HarmBAM-HarmM-HarmG-HarmW
#===================================================================
#===================================================================
#==== 3-18/ Finances publiques ===============================================
#==== 3-18-1/ Recettes  =================================================
RG_TR  = RG_TRYCL- RTVA_TR*0.3
RG_TRYCL  = RIR_TR  + RIS_TR  + RAID_TR + DDA+DDNA + RTVA_TR  + RTIC_TR + AURF_TR+RNF_TR+RCST_TR
#===impots sur les produits==
RTVA_TR  = aRTVA_TR *TVANA 
#===impots sur le revenu et le patrimoine==
RIR_TR  = aRIR_TR * TRPMNA 
TRPMNA =TTRPMNA*RDMNA/(1-TTRPMNA)
RIS_TR  = aRIS_TR *TRPS 
TRPS=TRPSHIC+TRPIC+TRPBAM
TRPSHIC=TTRPSHIC*EPARSOEIHIC/(1-TTRPSHIC)
TRPIC=TTRPIC*EPARIC/(1-TTRPIC)
TRPBAM=TTRPBAM*(EBEBAM+IntBAM_R-IntBAM_V+HarmBAM)
#====3-18-2/  Depenses ============================================
DG_TR  = DBSG_TR+ DCOMP_TR + DINTG_TR  + DIG_TR  
DBSG_TR = DPERG_TR + DABSG_TR
DPERG_TR = aDPERG_TR *RSG
DABSG_TR = aDABSG_TR *(CGA+CGMNA+PSG)
DCOMP_TR=aDCOMP_TR*(SPTA+SPTNA)
DIG_TR  = aDIG_TR *  IG 
DINTG_TR  = aDINTG_TR* (IntG_V+SIFIMG)
#==== 3-18-3/ Solde budgetaire et financement =============================
SBUDG  = RG_TR  - DG_TR  + SCST  
SBUDGPIB= SBUDG/PIB*100
#================================================================
#================================================================
#===3-18-4/ Encours credits et depots==========================================
# Menages
#============
#==Numeraire et Depots et opcvm non monetaires
MFIDMdot=betaMFIDM*(MFIDMC-MFIDM) #Dirhams
MFIDMC = lambdaMFIDM*RDM #Dirhams
DBM_ICdot=CAPFM-MFIDMdot+CBIC_Mdot+OMdot #Dirhams
#==Credits
CBIC_Mdot=betaCBIC_M*(CBIC_MC-CBIC_M) #Dirhams
# CBIC_MC=(100/TintD_IMB)*rhoCBIC_M*RDMNA*(1+rhoCBIC_M2*TND0709)*(1+rhoCBIC_M3*DUM1011) #Dirhams
CBIC_MC=(rhoCBIC_M+rhoCBIC_M2*(rhoCBIC_M3-4*TintD/100))*RDMNA
#==Autres	
# OMdot=lambdaOM*RDMNA #Dirhams
OMdot=((1/(1.0 + exp(-elastOM*(t - initOM))))*(UBOM - LBOM) + LBOM)*PIB

#============
# APU
#============
#============
#== dette APU: totale, domestique et exterieure
DetteTotG=DetteDG+DetteFXG  #Dirhams
DetteDGdot=Bg_Ndot+CBIC_Gdot  #Dirhams
DetteFXG=FXDetteG*TCEN  #Dirhams
FXDetteGC=max(thetaFXDetteG*DetteDG/TCEN,FXDetteG+(phiAoR*MBS-AoR)/TCEN)  #Devise
#==Numeraire et Depots et opcvm non monetaires
DBG_ICdot=betaDBG_IC*(DBG_ICC-DBG_IC) #Dirhams
DBG_ICC= lambdaDBG_IC*VAG #Dirhams
#==Titres autres qu actions (obligations APU)
Bg_Ndot=-CAPFG+DBG_ICdot-Dg_Ndot-BgWdot-CBIC_Gdot-Cw_Gdot+E_BAMdot+EG_ICdot+EG_SHICdot-OGdot  #Dirhams
BgWdot=betaBgW*(BgWC-BgW) #Dirhams
BgWC=DetteFXG-Cw_G #Dirhams
#==Credits
Cw_Gdot=betaCw_G*(Cw_GC-Cw_G) #Dirhams
Cw_GC=FXCw_G*TCEN #Dirhams
FXCw_Gdot=betaFXCw_G*(FXCw_GC-FXCw_G) #Devises
FXCw_GC=thetaFXCw_G*FXDetteG #Devises
CBIC_Gdot=betaCBIC_G*(CBIC_GC-CBIC_G)*(1+lambdaCBIC_G3*DUM19) #Dirhams
CBIC_GC=thetaCBIC_G*PRODG*(1+lambdaCBIC_G1*DUM08)*(1+lambdaCBIC_G2*DUM19) #Dirhams
#==Actions titres de participation et opcvm non monetaires
#E_BAMdot EG_ICdot EG_SHICdot sont exogenes  #Dirhams
E_BAMdot = betaEBAM*(lambdaEBAM*VSKNA*IPI-E_BAM)
EG_ICdot = betaEGIC*(lambdaEGIC*VSKNA*IPI-EG_IC)
# EG_SHICdot = betaEGSHIC*(lambdaEGSHIC*VSKNA*IPI-EG_SHIC)
EG_SHICdot = betaEGSHIC*(0.35*E_SHIC_N-EG_SHIC)
#==Autres
OGdot=((1/(1.0 + exp(-elastOG*(t - initOG))))*(UBOG - LBOG) + LBOG)*PIB
Dg_Ndot = DgSHICdot+DgICdot
#================
# Reste du Monde
#================
#==Avoirs de reserves etrangeres
ARdot=-CAPFW+DwW_ICdot+BshicgWdot+BfW_SHICdot-CwIC_Wdot+CwWdot+EW_Ndot-OWdot #Dirhams
#==Or monetaire et DTS
# ORDTSdot=betaORDTS*(ORDTSC-ORDTS) #Dirhams
# ORDTSC=lambdaORDTS*AR*(1+lambdaORDTS1*TND0913)*(1+lambdaORDTS2*DUM15)*(1+lambdaORDTS3*DUM0708) #Dirhams
#==Numeraire et Depots et opcvm non monetaires
Dw_Wdot=ARdot-ORDTSdot-Bwdot  #Dirhams
DwW_ICdot = betaDwW_IC*(DwW_ICC-DwW_IC) #Dirhams
# DwW_ICC=lambdaDwW_IC*(XBS+MBS)*(1+lambdaDwW_IC1*DUM0711)*(1+lambdaDwW_IC2*DUM1819) #Dirhams
DwW_ICC=lambdaDwW_IC*(XBS+MBS)
#==Titres autres qu actions (obligations RDM)
BshicgWdot=BgWdot+BshicWdot
Bwdot = betaBw*(BwC-Bw) #Dirhams
# BwC=AR*gammaBw1*(1+gammaBw3*DUM1819)*(1+gammaBw4*DUM08)*((1+TintBw/100)/(1+TintAR/100))^gammaBw2 #Dirhams
BwC=AR*gammaBw1*(1+gammaBw3*DUM07)*(1+gammaBw4*DUM12)*(1+gammaBw5*DUM1517)*((1+TintBw/100)/(1+TintAOR/100))^gammaBw2 #Dirhams
#==Credits et obligations ?mises par les SHIC d?tenues par le RDM
BfW_SHICC=lambdaBfW_SHIC*max(0,LNDSHIC_N)*(1+lambdaBfW_SHIC1*DUM07)*(1+lambdaBfW_SHIC2*DUM1012)*(1+lambdaBfW_SHIC3*DUM1415) #Dirhams
BfW_SHICdot=betaBfW_SHIC*(BfW_SHICC-BfW_SHIC)#Dirhams
CwIC_Wdot=betaCwIC_W*(CwIC_WC-CwIC_W) #Dirhams
CwIC_WC=lambdaCwIC_W*(XBS+MBS)#Dirhams
CwWdot=CwW_SHICdot+CwW_ICdot+Cw_Gdot #Dirhams
#==Actions titres de participation et opcvm non monetaires
EW_Ndot=EW_SHIC_Ndot+EW_IC_Ndot #Dirhams
#EW_SHIC_Ndot et EW_IC_Ndot exogènes
EW_IC_Ndot = betaEWIC*(lambdaEWIC*VSKNA*IPI-EW_IC_N)
EW_SHIC_Ndot = betaEWSHIC*(0.58*E_SHIC_N-EW_SHIC_N)
# EW_SHIC_Ndot = betaEWSHIC*(lambdaEWSHIC*VSKNA*IPI-EW_SHIC_N)
#==Autres
OWdot=((1/(1.0 + exp(-elastOW*(t - initOW))))*(UBOW - LBOW) + LBOW)*PIB
#=================================================
# Societes et entreprises Hors Instituts de credits
#=================================================
#==Numeraire et Depots et opcvm non monetaires
DBSHIC_ICdot=betaDBSHIC_IC*(DBSHIC_ICC-DBSHIC_IC) #Dirhams
DBSHIC_ICC=lambdaDBSHIC_IC*(RSMNAHIC_B+CSEMNAHIC)*(1+lambdaDBSHIC_IC1*DUM08)*(1+lambdaDBSHIC_IC2*DUM1112)+max(0,0.1*(VSKNA*IPI)-LNDSHIC_N) #Dirhams
DgSHICdot=betaDgSHIC*(DgSHICC-DgSHIC) #Dirhams
DgSHICC=lambdaDgSHIC*(RSMNAHIC_B+CSEMNAHIC)*(1+lambdaDgSHIC1*TND0810)*(1+lambdaDgSHIC2*DUM19)+max(0,0.1*(VSKNA*IPI)-LNDSHIC_N) #Dirhams
DwSHIC_ICdot=betaDwSHIC_IC*(DwSHIC_ICC-DwSHIC_IC) #Dirhams
DwSHIC_ICC=(DBSHIC_IC+DgSHIC)*lambdaDwSHIC_IC*(1+lambdaDwSHIC_IC2*TND0711)*(1+lambdaDwSHIC_IC3*DUM19)*(VXBS/VPRODMNA)^lambdaDwSHIC_IC1++max(0,0.1*(VSKNA*IPI)-LNDSHIC_N)#Dirhams
#==Titres autres qu actions (obligations APU)
BgSHICdot=betaBgSHIC*(BgSHICC-BgSHIC) #Dirhams
BgSHICC=lambdaBgSHIC*VAMNAHIC+max(0,0.1*(VSKNA*IPI)-LNDSHIC_N) #Dirhams
#==Credits et obligations priv?es 
LNDSHIC_Ndot=-CAPFSOEIHIC+DBSHIC_ICdot+DgSHICdot+DwSHIC_ICdot+BgSHICdot-BshicWdot-CwIC_SHICdot-CwW_SHICdot+ESHIC_ICdot-E_SHIC_Ndot-OSHICdot #Dirhams
CwIC_SHICdot=betaCwIC_SHIC*(CwIC_SHICC-CwIC_SHIC) #Dirhams
CwIC_SHICC=lambdaCwIC_SHIC*(XBS+MBS)*(1+lambdaCwIC_SHIC1*TND0710)*(1+lambdaCwIC_SHIC2*TND1119)*(1+lambdaCwIC_SHIC3*DUM16) #Dirhams
#CwW_SHICdot exogene
#==Actions titres de participation et opcvm non monetaires
# E_SHIC_Ndot=EIC_SHICdot+EG_SHICdot+EW_SHIC_Ndot #Dirhams
E_SHIC_Ndot=1*(0.22*IPI*VSKNA-E_SHIC_N) #Dirhams
# ESHIC_ICdot exogene
# ESHIC_ICdot = betaESHICIC*(lambdaESHICIC*VSKNA*IPI-ESHIC_IC)
ESHIC_ICdot = betaESHICIC*(0.34*(BgSHIC+DwSHIC_IC+DgSHIC+DBSHIC_IC)-ESHIC_IC)
#==Autres
OSHICdot=((1/(1.0 + exp(-elastOSHIC*(t - initOSHIC))))*(UBOSHIC - LBOSHIC) + LBOSHIC)*PIB
#======================
# Instituts de Credits
#======================
#==Numeraire et Depots et opcvm non monetaires
DBIC_BAMdot=betaDBIC_BAM*(DBIC_BAMC-DBIC_BAM)
DBIC_BAMC=lambdaDBIC_BAM*DB_IC_N*(1+lambdaDBIC_BAM1*TND0913)*(1+lambdaDBIC_BAM2*TND0708)
DgICdot=betaDgIC*(DgICC-DgIC)
DgICC=lambdaDgIC*DB_IC_N*(1+lambdaDgIC1*DUM10)*(1+lambdaDgIC2*DUM08)*(1+lambdaDgIC3*DUM12)*(1+lambdaDgIC4*DUM14)
DB_IC_Ndot=DBSHIC_ICdot+DBM_ICdot+DBG_ICdot
DwIC_Wdot=betaDwIC_W*(DwIC_WC-DwIC_W)
# DwIC_WC=lambdaDwIC_W*(Dw_IC_N+CwIC_N)*(1+lambdaDwIC_W1*DUM0708)*(1+lambdaDwIC_W2*DUM0910)*(1+lambdaDwIC_W3*DUM1112)
DwIC_WC=thetaFX*(Dw_IC_N+CwIC_N)*(lambdaDwIC_W1*((1+TintBw/100)/(1+TintAOR/100))^lambdaDwIC_W2)
Dw_IC_Ndot=DwSHIC_ICdot+DwW_ICdot
#==Titres autres qu actions
# BgICdot=omegaBgIC*Bg_Ndot
# omegaBgIC=omegaBgIC0+omegaBgIC1*(TintBg-TintD)^sigmaBg
BgICdot=Bg_Ndot-BgBAMdot-BgSHICdot
BwICdot=betaBwIC*(BwICC-BwIC)
# BwICC=(1-lambdaDwIC_W)*(1+lambdaDwIC_W4*TND1011)*(1+lambdaDwIC_W5*TND1516)*(Dw_IC_N+CwIC_N)
BwICC=thetaFX*(1-(lambdaDwIC_W1*((1+TintBw/100)/(1+TintAOR/100))^lambdaDwIC_W2))*(Dw_IC_N+CwIC_N)
#==Credits et obligations privees
CBBAM_ICdot=-CAPFIC+DBIC_BAMdot-DB_IC_Ndot+DgICdot+DwIC_Wdot-Dw_IC_Ndot+BgICdot+BwICdot+LNDIC_Ndot+CwIC_Ndot-CwW_ICdot-E_IC_Ndot+EIC_SHICdot-OICdot
LNDIC_Ndot=LNDSHIC_Ndot+CBIC_Mdot+CBIC_Gdot-BfW_SHICdot
CwIC_Ndot=CwIC_SHICdot+CwIC_Wdot
CwW_ICdot =betaCwWIC*(lambdaCwWIC*CwIC_N-CwW_IC)
#==Actions titres de participation et opcvm non monetaires
# EIC_SHICdot exogene
# EIC_SHICdot = betaEICSHIC*(lambdaEICSHIC*VSKNA*IPI-EIC_SHIC)
EIC_SHICdot = E_SHIC_Ndot-EG_SHICdot-EW_SHIC_Ndot
E_IC_Ndot=ESHIC_ICdot+EG_ICdot+EW_IC_Ndot
#==Autres
OICdot=((1/(1.0 + exp(-elastOIC*(t - initOIC))))*(UBOIC - LBOIC) + LBOIC)*PIB
#================
# Banque centrale
#================
#==Avoirs officiels de reserves
AoR=ORDTS+DwBAM_W+BwBAM #Dirhams
BwBAMdot=Bwdot-BwICdot #Dirhams
DwBAM_Wdot=Dw_Wdot-DwIC_Wdot #Dirhams
#==Numeraire et Depots et opcvm non monetaires
DB_BAMdot=DBIC_BAMdot+MFIDMdot #Dirhams
#==Titres autres qu actions
# BgBAMdot=max(0,Bg_Ndot-BgSHICdot-BgICdot)
#BgBAMdot=0 #Dirhams
#==Credits
#CBBAM_ICdot deja modelisee dans IC
#==Actions titres de participation et opcvm non monetaires
#E_BAMdot exogene
#==Autres
# OBAMdot=-CAPFBAM+ORDTSdot-DB_BAMdot+BgBAMdot+CBBAM_ICdot-E_BAMdot #Dirhams ? Verifier que cette egalite est respectee si c'est le cas, le cadre comptable est bon, sinon problème!
OBAMdot = -(OSHICdot+OICdot+OMdot+OGdot+OWdot)
#===========================================================
#================3-18-5/Taux d'interets===========================================
#======Taux d'interets domestiques 
#Mesure les besoins en avances en ratio des actifs des banques. Plus le ratio est élevé plus les banques sont en besoins de liquidité
ratioLiquidite=CBBAM_IC/(BgIC + BwIC + LNDIC_N + CwIC_N)
#Un plut haut besoin de liquidité doit mener à une hause des taux créditeurs
TintCC = TintBgN - rhoTintC1/(1+exp(-rhoTintC2*(1+rhoTintC5*DUM1618)*ratioLiquidite))
#Un plut haut besoin de liquidité doit mener à une hause des taux interbancaire
TintBC = TintBAM + alphaTintB1/(1+exp(-alphaTintB2*ratioLiquidite))


#Mesure le risque de l'emprunteur capturé par le niveau de dette sur l'EBE
ratioRisquePr=(BshicW+LNDSHIC_N+CwIC_SHIC+CwW_SHIC)/EBEMNAHIC
#Plus le risque est élevé, plus la marge sur le taux débiteur sera élevée
TintDC = AFC + xiTintD0*(1+xiTintD3*DUM0712) + xiTintD1/(1+exp(-xiTintD2*ratioRisquePr))

#Mesure le risque emprunteur capturé par le niveau de dette (totale) sur PIB
ratioRisquePub = (Bg_N+BgW+CBIC_G+Cw_G)/PIB
#Plus le risque est élevé plus la marge sera élevée
TintBgC = TintBAM + nuTintBg1/(1+exp(-nuTintBg2*ratioRisquePub))

#Cout moyens de financement des banques
AFC=(TintBAM*CBBAM_IC+TintC*DB_IC_N)/(CBBAM_IC+DB_IC_N)

TintC = TintCN + rhoTintC3*DUM08 + rhoTintC4*DUM10
TintD = TintDN
TintBg = TintBgN + nuTintBg3*DUM10 +nuTintBg4*DUM15
TintB = TintBN + alphaTintB3*DUM08 +alphaTintB4*DUM14

#======Taux d'interets extérieurs 

#DeficitSurPIB=-100*SBUDG/PIB+2
DetteEEPSurPIB=(BshicW)/PIB
DetteSurPIB=(Bg_N+BgW+CBIC_G+Cw_G)/PIB
NIIP=(-Dw_W+DwW_IC+BgW-Bw+BfW_SHIC-CwIC_W+CwW+EW_IC_N+EW_SHIC_N)/PIB

TintBgWC = TintBw + muTintBgW1/(1+exp(-muTintBgW2*NIIP)) + muTintBgW3/(1+exp(muTintBgW4*SBUDGPIB))
TintBgW = TintBgWN + muTintBgW5*DUM09

TintBshicWC = TintBw + muTintBshicW1/(1+exp(-muTintBshicW2*NIIP))+ muTintBshicW3/(1+exp(-muTintBshicW4*DetteEEPSurPIB))
TintBshicW = TintBshicWN + muTintBshicW5*DUM14 + muTintBshicW6*DUM15

TintCDgWC = TintRfw_cd + muTintCDgW1/(1+exp(-muTintCDgW2*NIIP)) + muTintCDgW3/(1+exp(-muTintCDgW4*SBUDGPIB))
TintCDgW = TintCDgWN + muTintCDgW5*DUM0809+ muTintCDgW6*DUM12

TintCDshicWC = TintRfw_cd + muTintCDshicW1/(1+exp(-muTintCDshicW2*NIIP))+ muTintCDshicW3/(1+exp(-muTintCDshicW4*DetteEEPSurPIB))
TintCDshicW = TintCDshicWN + muTintCDshicW5*DUM09

##exogenous variables
DUM07=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM08=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM09=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM10=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM11=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM12=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1217=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM13=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM14=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1415=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM15=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM16=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM17=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM18=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM19=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
DUM0708=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0809=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0810=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0815=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0709=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0710=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0712=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0714=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0715=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0716=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0910=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM0911=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1011=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1012=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1112=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1114=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1214=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1215=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1314=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1417=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1419=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
DUM1516=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1517=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1618=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
DUM1619=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
DUM1719=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
DUM1819=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
DUM1219=c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
TND0708=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
TND0709=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
TND0710=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
TND0711=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
TND0810=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
TND0814=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
TND0913=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
TND1119=c(9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9)
TND1519=c(5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5)
aCONS_SF=c(0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118)
aX_SF=c(0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118)
aBrIDE_R=c(0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555,0.03505555)
aBrIDE_V=c(0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562,0.01190562)
aDABSG_TR=c(0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641,0.7046641)
aDCOMP_TR=c(1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243,1.115243)
aDIG_TR=c(1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049,1.210049)
aDINTG_TR=c(1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664,1.0664)
aDPERG_TR=c(0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106,0.8191106)
aMTCIC=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
aMTCNA=c(-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973,-0.01287973)
aPFATA=c(0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717,0.00902717)
aPFATIC=c(0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905,0.004462905)
aPFATMNA=c(-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521,-0.01551521)
aPFATG=c(0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691,0.09026691)
aRassG=c(0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553,0.01259553)
aRassIC_R=c(0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569,0.02552569)
aRassM=c(0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411,0.7956411)
aRassSHIC_R=c(0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269,0.0005632269)
aRassSHIC_V=c(0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409,0.005337409)
aRassW=c(0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327,0.06071327)
aRDMABMK=c(0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578,0.1163578)
aRDMACS=c(0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692,0.08795692)
aRDMADT=c(0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556,0.06538556)
aRDMAEOD=c(0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333,0.00465333)
aRDMAFM=c(0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205,0.1537205)
aRDMAGON=c(0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517,0.01190517)
aRDMALSH=c(0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132,0.01195132)
aRDMAMS=c(0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506,0.1151506)
aRDMAO=c(0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602,0.0674602)
aRDMARSK=c(0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248,0.1473248)
aRDMASM=c(0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791,0.1287791)
aRDMATTH=c(0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468,0.08935468)
aRDMNABMK=c(0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344,0.04178344)
aRDMNACS=c(0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466,0.3429466)
aRDMNADT=c(0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843,0.02870843)
aRDMNAEOD=c(0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391,0.01333391)
aRDMNAFM=c(0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446,0.0734446)
aRDMNAGON=c(0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295,0.01624295)
aRDMNALSH=c(0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075,0.01972075)
aRDMNAMS=c(0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502,0.08202502)
aRDMNAO=c(0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401,0.04558401)
aRDMNARSK=c(0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264,0.1578264)
aRDMNASM=c(0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713,0.06436713)
aRDMNATTH=c(0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167,0.1140167)
aRIR_TR=c(1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064,1.466064)
aRIS_TR=c(0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941,0.7253941)
aRlyG_V=c(0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506,0.4858506)
aRlyM_R=c(0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494,0.5141494)
aRlyM_V=c(0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142,0.001764142)
aRlyS_V=c(0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863,0.001599863)
aRlyW_V=c(0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227,0.001138227)
ARRTOUR=c(12932260,13416700,13919290,14440710,14981660,15542870,16125110,16729150,17355830,18005980,18680480,19380250,20106240,20859410,21640810,22451470,23292510,24165040,25070270,26009400,26983710,27994520,29043190,30131150,31259860,32430860,33645720,34906090,36213670,37570230,38977620,40437720)
aRSIC=c(0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428,0.07357428)
aRTVA_TR=c(1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934,1.001934)
aSSNF_Flux=c(0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564,0.9782564)
aTRFM_A=c(0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705,0.920705)
aUIAA=c(0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589,0.1333589)
aUIAIA=c(0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429,0.09620429)
aUIANA=c(0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642,0.1584642)
aUIANAHIA=c(0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945,0.09865945)
aUIA_SF=c(0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937,0.1762937)
aUIGA=c(0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567,0.001283567)
aUIGIA=c(0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489,0.005316489)
aUIGNAHIAHIC=c(0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101,0.2312101)
aUIICNAHIAHIC=c(0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103,0.4047103)
aUIG_SF=c(0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706,0.1622706)
aUIMNAA=c(0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517,0.06224517)
aUIMNAIA=c(0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582,0.02624582)
aUIMNANAHIA=c(0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236,0.3916236)
aUIMNAHIC_SF=c(0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149,0.218149)
aUINAA=c(0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413,0.06245413)
aUINANA=c(0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471,0.4661471)
aCG_SF=c(0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118,0.212118)
AURF_TR=c(14948,15545.92,16167.76,16814.47,17487.05,18186.53,18913.99,19670.55,20457.37,21275.66,22126.69,23011.76,23932.23,24889.52,25885.1,26920.5,27997.32,29117.22,30281.91,31493.18,32752.91,34063.03,35425.55,36842.57,38316.27,39848.92,41442.88,43100.59,44824.62,46617.6,48482.31,50421.6)
BgBAMdot=c(205.8004,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
BshicWdot=c(2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151,2928.151)
CGA=c(2100,2186.1,2275.73,2369.035,2466.165,2567.278,2672.537,2782.111,2896.177,3014.92,3138.532,3267.212,3401.168,3540.616,3685.781,3836.898,3994.211,4157.973,4328.45,4505.917,4690.659,4882.976,5083.178,5291.589,5508.544,5734.394,5969.504,6214.254,6469.038,6734.269,7010.374,7297.799)
CGG=c(215317,224140.9,233330.8,242897.6,252856.5,263223.8,274016.2,285251,296946.5,309121.4,321795.6,334989.4,348724.1,363022,377906,393400.3,409529.9,426320.8,443800.1,461996.1,480938.1,500656.8,521183.9,542552.6,564797.4,587954.2,612060.5,637155.2,663278.7,690473.3,718782.9,748253.2)
CGMNA=c(5550,5777.55,6014.43,6261.021,6517.723,6784.95,7063.133,7352.721,7654.183,7968.004,8294.692,8634.775,8988.8,9357.341,9740.992,10140.37,10556.13,10988.93,11439.48,11908.49,12396.74,12905.01,13434.11,13984.91,14558.29,15155.18,15776.55,16423.39,17096.74,17797.71,18527.42,19287.04)
choc=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
CwW_SHICdot=c(-998.4548,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34,6310.34)
DivBAM=c(984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529,984.529)
EXM_BPCN=c(-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5,-5024.5)
HarmBAM=c(530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758,530.9758)
HarmG=c(-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407,-3105.407)
HarmIC=c(6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794,6463.794)
HarmM=c(-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879,-402.2879)
HarmW=c(1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619,1489.619)
IPCUE=c(1.584781,1.608553,1.632681,1.657172,1.682029,1.70726,1.732869,1.758862,1.785244,1.812023,1.839203,1.866792,1.894793,1.923215,1.952064,1.981344,2.011065,2.041231,2.071849,2.102927,2.134471,2.166488,2.198985,2.23197,2.265449,2.299431,2.333923,2.368931,2.404465,2.440532,2.47714,2.514298)
IPCUSA=c(1.535115,1.558142,1.581514,1.605237,1.629315,1.653755,1.678562,1.70374,1.729296,1.755235,1.781564,1.808287,1.835412,1.862943,1.890887,1.91925,1.948039,1.97726,2.006919,2.037022,2.067578,2.098591,2.13007,2.162021,2.194452,2.227368,2.260779,2.294691,2.329111,2.364048,2.399508,2.435501)
IPDBAM=c(18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775,18.775)
IPPRODA=c(1.150732,1.168048,1.185626,1.203468,1.221578,1.239961,1.258621,1.277561,1.296787,1.316301,1.33611,1.356216,1.376625,1.397342,1.41837,1.439714,1.46138,1.483371,1.505694,1.528352,1.551352,1.574697,1.598394,1.622448,1.646863,1.671646,1.696802,1.722336,1.748255,1.774564,1.801268,1.828375)
IPPRODA_Infl=c(0.06791273,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854,0.01504854)
IPS_Flux=c(0.5416289,0.5497533,0.5579996,0.5663696,0.5748651,0.5834881,0.5922404,0.601124,0.6101409,0.619293,0.6285824,0.6380111,0.6475813,0.657295,0.6671545,0.6771618,0.6873192,0.697629,0.7080934,0.7187148,0.7294955,0.740438,0.7515446,0.7628177,0.77426,0.7858739,0.797662,0.8096269,0.8217713,0.8340979,0.8466094,0.8593085)
IPSA_Flux=c(0.5027851,0.5103269,0.5179818,0.5257515,0.5336378,0.5416424,0.549767,0.5580135,0.5663837,0.5748795,0.5835027,0.5922552,0.601139,0.6101561,0.6193085,0.6285981,0.6380271,0.6475975,0.6573114,0.6671711,0.6771787,0.6873364,0.6976464,0.7081111,0.7187328,0.7295138,0.7404565,0.7515633,0.7628368,0.7742793,0.7858935,0.7976819)
IPSNA_Flux=c(0.6538139,0.6636211,0.6735754,0.6836791,0.6939343,0.7043433,0.7149084,0.725632,0.7365165,0.7475643,0.7587777,0.7701594,0.7817118,0.7934375,0.805339,0.8174191,0.8296804,0.8421256,0.8547575,0.8675789,0.8805925,0.8938014,0.9072085,0.9208166,0.9346288,0.9486483,0.962878,0.9773212,0.991981,1.006861,1.021964,1.037293)
lambdaDBG_IC=c(0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509,0.5312509)
lambdaOG=c(0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866,0.007777866)
lambdaOIC=c(0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836,0.03743836)
lambdaOSHIC=c(-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343,-0.01509343)
lambdaOW=c(-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325,-0.02942325)
LG=c(905.571,914.6267,923.773,933.0107,942.3408,951.7642,961.2818,970.8947,980.6036,990.4096,1000.314,1010.317,1020.42,1030.624,1040.93,1051.34,1061.853,1072.472,1083.196,1094.028,1104.969,1116.018,1127.179,1138.45,1149.835,1161.333,1172.947,1184.676,1196.523,1208.488,1220.573,1232.779)
MuA=c(0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672,0.08106672)
MuG=c(0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958,0.2004958)
MuIC=c(0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932,0.237932)
MuMNA=c(0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689,0.2541689)
ORDTSdot=c(1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326,1583.326)
pA_SIF=c(0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717,0.00588717)
PACT=c(12081.82,12266.22,12451.23,12636.79,12822.83,13009.31,13196.16,13383.32,13570.72,13758.3,13946,14133.75,14321.48,14509.13,14696.63,14883.9,15070.88,15257.5,15443.68,15629.35,15814.45,15998.88,16182.59,16365.5,16547.53,16728.6,16908.64,17087.58,17265.34,17441.83,17616.99,17790.74)
Pble_dolr=c(236.9501,240.5043,244.1119,247.7736,251.4902,255.2625,259.0915,262.9779,266.9225,270.9264,274.9903,279.1151,283.3018,287.5514,291.8646,296.2426,300.6862,305.1965,309.7745,314.4211,319.1374,323.9245,328.7834,333.7151,338.7208,343.8016,348.9587,354.193,359.5059,364.8985,370.372,375.9276)
Pble_dolr_2007=c(304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171,304.2171)
pCIG_SIF=c(0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039,0.6564039)
PFATBAM=c(153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697,153.697)
phiAoR=c(0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249,0.4232249)
POP=c(35586.62,35951.66,36313.19,36670.22,37022.39,37369.65,37712.15,38049.58,38381.17,38705.84,39022.52,39329.98,39628,39916.33,40194.74,40463.01,40720.95,40968.34,41205.01,41430.77,41645.45,41848.9,42040.98,42221.55,42390.48,42547.67,42693.02,42826.44,42947.86,43057.22,43154.47,43239.58)
POPR=c(13210.5,13168.22,13124.29,13078.63,13031.34,12982.69,12932.95,12882.19,12830.4,12777.52,12723.54,12668.44,12612.46,12555.83,12498.79,12441.58,12384.43,12327.57,12271.24,12215.67,12161.07,12107.69,12055.73,12005.4,11956.92,11910.48,11866.28,11824.51,11785.34,11748.94,11715.48,11685.11)
POPU=c(22376.12,22783.44,23188.9,23591.58,23991.05,24386.96,24779.2,25167.39,25550.77,25928.31,26298.98,26661.55,27015.54,27360.5,27695.94,28021.43,28336.52,28640.77,28933.77,29215.1,29484.38,29741.22,29985.25,30216.14,30433.56,30637.19,30826.73,31001.93,31162.52,31308.28,31438.99,31554.48)
PPET=c(64.06454,65.02551,66.00089,66.9909,67.99577,69.0157,70.05094,71.1017,72.16823,73.25075,74.34951,75.46476,76.59673,77.74568,78.91186,80.09554,81.29697,82.51643,83.75417,85.01049,86.28564,87.57993,88.89363,90.22703,91.58044,92.95414,94.34846,95.76368,97.20014,98.65814,100.138,101.6401)
PPET_2007=c(72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617,72.69617)
PRODA=c(191976,200710.9,209843.3,219391.1,229373.4,239809.9,250721.3,262129.1,274056,286525.5,299562.4,313192.5,327442.8,342341.4,357917.9,374203.2,391229.4,409030.4,427641.3,447098.9,467441.9,488710.6,510946.9,534195,558500.8,583912.6,610480.7,638257.5,667298.2,697660.3,729403.9,762591.7)
PRODBAM_NV=c(1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861,1464.861)
PSG_STR=c(75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442,75.46442)
PSIC_STR=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
PSSHIC_STR=c(19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827,19.8827)
PSW_STR=c(4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887,4.652887)
PXAPH_dolr=c(675.327,682.0248,692.2551,702.639,713.1785,723.8762,734.7344,745.7554,756.9417,768.2958,779.8203,791.5176,803.3903,815.4412,827.6728,840.0879,852.6892,865.4796,878.4618,891.6387,905.0133,918.5885,932.3673,946.3528,960.5481,974.9563,989.5807,1004.424,1019.491,1034.783,1050.305,1066.059)
PXAPH_dolr_2007=c(488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831,488.5831)
PXENG_dolr=c(322.0348,329.9747,334.9243,339.9481,345.0474,350.2231,355.4764,360.8086,366.2207,371.714,377.2897,382.9491,388.6933,394.5237,400.4416,406.4482,412.5449,418.7331,425.0141,431.3893,437.8601,444.428,451.0944,457.8609,464.7288,471.6997,478.7752,485.9568,493.2462,500.6449,508.1545,515.7769)
PXENG_dolr_2007=c(376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087,376.3087)
PXPH_dolr=c(80.0651,81.26608,82.48507,83.72235,84.97818,86.25285,87.54665,88.85985,90.19274,91.54563,92.91882,94.3126,95.72729,97.1632,98.62065,100.1,101.6015,103.1255,104.6724,106.2424,107.8361,109.4536,111.0954,112.7619,114.4533,116.1701,117.9126,119.6813,121.4765,123.2987,125.1482,127.0254)
PXPH_dolr_2007=c(53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136,53.36136)
RAID_TR=c(4709,4897.36,5093.254,5296.985,5508.864,5729.219,5958.387,6196.723,6444.592,6702.375,6970.47,7249.289,7539.261,7840.831,8154.464,8480.643,8819.869,9172.663,9539.57,9921.153,10318,10730.72,11159.95,11606.35,12070.6,12553.42,13055.56,13577.78,14120.89,14685.73,15273.16,15884.09)
RassIC_V=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
RCST_TR=c(3796,3947.84,4105.754,4269.984,4440.783,4618.414,4803.151,4995.277,5195.088,5402.892,5619.007,5843.768,6077.518,6320.619,6573.444,6836.382,7109.837,7394.23,7689.999,7997.599,8317.503,8650.204,8996.212,9356.06,9730.303,10119.51,10524.3,10945.27,11383.08,11838.4,12311.94,12804.41)
RETR=c(1176114,1214127,1253369,1293879,1335698,1378869,1423435,1469442,1516935,1565964,1616577,1668827,1722764,1778446,1835927,1895265,1956522,2019759,2085039,2152429,2221998,2293815,2367953,2444487,2523495,2605057,2689254,2776173,2865902,2958530,3054152,3152865)
RNF_TR=c(36405.25,37861.46,39375.92,40950.96,42589,44292.56,46064.26,47906.83,49823.1,51816.03,53888.67,56044.22,58285.98,60617.42,63042.12,65563.81,68186.36,70913.81,73750.36,76700.38,79768.39,82959.13,86277.49,89728.59,93317.74,97050.45,100932.5,104969.8,109168.6,113535.3,118076.7,122799.8)
RSBAM_B=c(697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462,697.4462)
RSE_RHMRE=c(17203.7,17575.35,17955.04,18342.92,18739.19,19144.02,19557.59,19980.09,20411.73,20852.69,21303.17,21763.39,22233.55,22713.86,23204.55,23705.85,24217.97,24741.15,25275.64,25821.68,26379.51,26949.39,27531.58,28126.35,28733.97,29354.71,29988.87,30636.72,31298.58,31974.72,32665.48,33371.16)
RSE_V=c(6806.1,7215.003,7648.472,8107.983,8595.102,9111.486,9658.893,10239.19,10854.35,11506.46,12197.76,12930.59,13707.44,14530.97,15403.97,16329.43,17310.48,18350.47,19452.95,20621.66,21860.59,23173.95,24566.21,26042.12,27606.7,29265.28,31023.51,32887.36,34863.2,36957.74,39178.12,41531.89)
RSG=c(155924,162313.9,168968.9,175896.8,183108.6,190616.2,198431.6,206567.4,215036.8,223853.4,233031.6,242586,252532.1,262886.1,273664.5,284884.9,296565.3,308724.6,321382.4,334559.2,348276.2,362555.7,377420.6,392895,409003.8,425773.1,443229.9,461402.4,480320,500013.3,520514,541855.1)
RTIC_TR=c(29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3,29900.3)
S_IDE_E=c(638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2,638438.2)
S_IDE_S=c(61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7,61844.7)
SCST=c(2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536,2287.536)
SRPM=c(-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7,-18423.7)
TAIPTIA=c(0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413,0.009896413)
TAIPTIA0=c(0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516,0.008833516)
TAIPTNAHIA=c(0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858,0.05548858)
TAIPTNAHIA0=c(0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808,0.03881808)
TAIPTIC=c(0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519,0.1518519)
TauH=c(0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168,0.0640168)
tAutreIC=c(-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667,-0.05901667)
tAutreM=c(0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817,0.006008817)
tAutreSHIC=c(0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926,0.00397926)
tAutreW=c(2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06,2.602346e-06)
TCDOLR=c(9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398,9.616398)
TCDOLR_2007=c(8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325,8.195325)
TCEURO=c(10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689,10.7689)
TCEURO_2007=c(11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762,11.21762)
tCSSG=c(0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257,0.4432257)
tCSSIC=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
tCSSSHIC=c(0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178,0.5527178)
tCSSW=c(0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466,0.004056466)
TDDA=c(0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261,0.1366261)
TDDIA=c(0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506,0.05163506)
TDDNAHIA=c(0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683,0.009737683)
TDEC=c(0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769,0.04000769)
tDivG=c(0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146,0.3190146)
tDivIC=c(0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325,0.6541325)
tDivM=c(0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191,0.3855191)
tDivSHIC=c(0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131,0.06335131)
tDivW=c(0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663,0.2954663)
thetaFXCw_G=c(0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917,0.7386917)
thetaFXDetteG=c(0.3365444,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111,0.3390111)
TintAOR=c(1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062,1.238062)
TintBAM=c(2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25,2.25)
TintBw=c(0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498,0.4336498)
TintD_IMB=c(5.1275,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121,5.394121)
TintRfw_cd=c(1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267,1.295267)
TIPDA=c(0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989,0.0005315989)
TIPDG=c(0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755,0.002292755)
TIPDIC=c(0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445,0.02319445)
TIPDMNAHIC=c(0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303,0.00977303)
tLNRMNA=c(0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177,0.03089177)
tLRNSMNA=c(0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575,0.3681575)
tLSAG=c(0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832,0.1772832)
tLSMNA=c(0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507,0.6009507)
tMBIA=c(0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794,0.08344794)
TPSM=c(13739.11,14323.15,14932.02,15566.78,16228.52,16918.39,17637.58,18387.35,19168.99,19983.86,20833.37,21718.99,22642.26,23604.77,24608.2,25654.29,26744.85,27881.76,29067.01,30302.64,31590.79,32933.71,34333.71,35793.22,37314.78,38901.02,40554.69,42278.66,44075.91,45949.57,47902.87,49939.21)
TRF_MRE=c(64779.1,68018.05,71418.96,74989.91,78739.4,82676.37,86810.19,91150.7,95708.23,100493.6,105518.3,110794.2,116334,122150.7,128258.2,134671.1,141404.7,148474.9,155898.6,163693.6,171878.2,180472.1,189495.8,198970.5,208919.1,219365,230333.3,241849.9,253942.4,266639.6,279971.5,293970.1)
TRFBAM=c(1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754,1388.754)
TSG_B=c(143.4266,147.8261,152.3635,157.0401,161.8602,166.8283,171.9489,177.2267,182.6664,188.2731,194.0519,200.008,206.147,212.4744,218.996,225.7177,232.6458,239.7865,247.1464,254.7321,262.5507,270.6093,278.9152,287.4761,296.2997,305.3941,314.7677,324.429,334.3868,344.6502,355.2287,366.1318)
TSPDIC=c(0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403,0.0008559403)
TSPDMNAHIC=c(0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679,0.000841679)
TSPTA=c(0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155,0.002763155)
TSPTIA=c(0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562,0.02425562)
TSPTIA0=c(0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432,0.03994432)
TSPTNAHIA=c(0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219,0.02833219)
TSPTNAHIA0=c(0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244,0.05906244)
TTICIA=c(0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309,0.0840309)
TTICIA0=c(0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409,0.06314409)
TTICNAHIA=c(0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978,0.05933978)
TTICNAHIA0=c(0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333,0.0602333)
tTNSMNA=c(1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159,1.662159)
tTRFIC=c(-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046,-0.02209046)
tTRFSHIC=c(0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769,0.01575769)
TRFW_A=c(4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9,4979.9)
TTRPBAM=c(0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323,0.2621323)
TTRPIC=c(0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573,0.3136573)
TTRPMNA=c(0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907,0.04498907)
TTRPSHIC=c(0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864,0.2614864)
TTVAIA=c(0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395,0.0502395)
TTVAIA0=c(0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583,0.04632583)
TTVANAF=c(0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038,0.1007038)
TTVANAHIAC=c(0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323,0.1476323)
TTVANAHIAC0=c(0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129,0.1466129)
TTVAIC=c(0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038,0.05450038)
UIBAMNAHIA=c(615.836,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689,695.689)
UIICNAHIA=c(17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45,17646.45)
VDMHOCP=c(1.186651,1.222251,1.258919,1.296686,1.335587,1.375654,1.416924,1.459432,1.503215,1.548311,1.59476,1.642603,1.691881,1.742638,1.794917,1.848764,1.904227,1.961354,2.020195,2.0808,2.143225,2.207521,2.273747,2.341959,2.412218,2.484585,2.559122,2.635896,2.714973,2.796422,2.880315,2.966724)
VIG=c(52622.7,54446.26,56333.01,58285.14,60304.92,62394.7,64556.89,66794.01,69108.65,71503.51,73981.35,76545.06,79197.61,81942.09,84781.66,87719.64,90759.43,93904.56,97158.68,100525.6,104009.1,107613.4,111342.6,115201,119193.1,123323.5,127597.1,132018.8,136593.7,141327.2,146224.7,151291.8)
VPRODA=c(166829.5,171834.4,176989.4,182299.1,187768.1,193401.1,199203.2,205179.3,211334.7,217674.7,224204.9,230931.1,237859,244994.8,252344.6,259915,267712.4,275743.8,284016.1,292536.6,301312.7,310352.1,319662.6,329252.5,339130.1,349304,359783.1,370576.6,381693.9,393144.7,404939.1,417087.2)
VSMIG=c(30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05,30879.05)
VTSAG_B=c(13.85031,13.85104,13.85178,13.85251,13.85325,13.85399,13.85472,13.85546,13.8562,13.85693,13.85767,13.85841,13.85914,13.85988,13.86061,13.86135,13.86209,13.86282,13.86356,13.8643,13.86504,13.86577,13.86651,13.86725,13.86798,13.86872,13.86946,13.87019,13.87093,13.87167,13.87241,13.87314)
VXOCP=c(35787.57,36667.01,37568.05,38491.24,39437.11,40406.22,41399.15,42416.48,43458.81,44526.75,45620.94,46742.02,47890.64,49067.49,50273.26,51508.66,52774.42,54071.28,55400.01,56761.39,58156.23,59585.34,61049.57,62549.79,64086.87,65661.72,67275.27,68928.47,70622.3,72357.75,74135.84,75957.63)
XSVpartouriste=c(6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144,6523.144)
lambdaDwW_IC=c(0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826,0.03282826)
lambdaCwIC_W=c(0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946,0.02567946)
lambdaBgSHIC=c(0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504,0.1706504)
VPIB_observe=c(1010147,1039948,1071338,1105375,1140736,1177312,1215138,1254263,1294752,1336695,1380268,1426329,1483383,1539890,1593359,1645244,1697057,1751996,1814267,1875367,1934807,1994945,2061590,2129403,2195753,2264376,2338466,2411698,2485692,2565349,2644634,2725450)
lambdaMFIDM=c(0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425,0.3338425)
VPRODA2=c(170166.1,175271.1,180529.2,185945.1,191523.5,197269.2,203187.2,209282.9,215561.3,222028.2,228689,235549.7,242616.2,249894.7,257391.5,265113.3,273066.7,281258.7,289696.4,298387.3,307338.9,316559.1,326055.9,335837.6,345912.7,356290.1,366978.8,377988.1,389327.8,401007.6,413037.8,425429)
PRODA2=c(195815.5,204725.1,214040.1,223778.9,233960.9,244606.1,255735.7,267371.7,279537.1,292256,305553.7,319456.3,333991.6,349188.2,365076.3,381687.3,399054,417211,436194.1,456040.9,476790.8,498484.8,521165.8,544878.9,569670.9,595590.9,622690.3,651022.7,680644.2,711613.5,743991.9,777843.6)
VISOEI_LEAP=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


##sampling time
DUM07=seq(0,31,1)
DUM08=seq(0,31,1)
DUM09=seq(0,31,1)
DUM10=seq(0,31,1)
DUM11=seq(0,31,1)
DUM12=seq(0,31,1)
DUM1217=seq(0,31,1)
DUM13=seq(0,31,1)
DUM14=seq(0,31,1)
DUM1415=seq(0,31,1)
DUM15=seq(0,31,1)
DUM16=seq(0,31,1)
DUM17=seq(0,31,1)
DUM18=seq(0,31,1)
DUM19=seq(0,31,1)
DUM0708=seq(0,31,1)
DUM0809=seq(0,31,1)
DUM0810=seq(0,31,1)
DUM0815=seq(0,31,1)
DUM0709=seq(0,31,1)
DUM0710=seq(0,31,1)
DUM0712=seq(0,31,1)
DUM0714=seq(0,31,1)
DUM0715=seq(0,31,1)
DUM0716=seq(0,31,1)
DUM0910=seq(0,31,1)
DUM0911=seq(0,31,1)
DUM1011=seq(0,31,1)
DUM1012=seq(0,31,1)
DUM1112=seq(0,31,1)
DUM1114=seq(0,31,1)
DUM1214=seq(0,31,1)
DUM1215=seq(0,31,1)
DUM1314=seq(0,31,1)
DUM1417=seq(0,31,1)
DUM1419=seq(0,31,1)
DUM1516=seq(0,31,1)
DUM1517=seq(0,31,1)
DUM1618=seq(0,31,1)
DUM1619=seq(0,31,1)
DUM1719=seq(0,31,1)
DUM1819=seq(0,31,1)
DUM1219=seq(0,31,1)
TND0708=seq(0,31,1)
TND0709=seq(0,31,1)
TND0710=seq(0,31,1)
TND0711=seq(0,31,1)
TND0810=seq(0,31,1)
TND0814=seq(0,31,1)
TND0913=seq(0,31,1)
TND1119=seq(0,31,1)
TND1519=seq(0,31,1)
aCONS_SF=seq(0,31,1)
aX_SF=seq(0,31,1)
aBrIDE_R=seq(0,31,1)
aBrIDE_V=seq(0,31,1)
aDABSG_TR=seq(0,31,1)
aDCOMP_TR=seq(0,31,1)
aDIG_TR=seq(0,31,1)
aDINTG_TR=seq(0,31,1)
aDPERG_TR=seq(0,31,1)
aMTCIC=seq(0,31,1)
aMTCNA=seq(0,31,1)
aPFATA=seq(0,31,1)
aPFATIC=seq(0,31,1)
aPFATMNA=seq(0,31,1)
aPFATG=seq(0,31,1)
aRassG=seq(0,31,1)
aRassIC_R=seq(0,31,1)
aRassM=seq(0,31,1)
aRassSHIC_R=seq(0,31,1)
aRassSHIC_V=seq(0,31,1)
aRassW=seq(0,31,1)
aRDMABMK=seq(0,31,1)
aRDMACS=seq(0,31,1)
aRDMADT=seq(0,31,1)
aRDMAEOD=seq(0,31,1)
aRDMAFM=seq(0,31,1)
aRDMAGON=seq(0,31,1)
aRDMALSH=seq(0,31,1)
aRDMAMS=seq(0,31,1)
aRDMAO=seq(0,31,1)
aRDMARSK=seq(0,31,1)
aRDMASM=seq(0,31,1)
aRDMATTH=seq(0,31,1)
aRDMNABMK=seq(0,31,1)
aRDMNACS=seq(0,31,1)
aRDMNADT=seq(0,31,1)
aRDMNAEOD=seq(0,31,1)
aRDMNAFM=seq(0,31,1)
aRDMNAGON=seq(0,31,1)
aRDMNALSH=seq(0,31,1)
aRDMNAMS=seq(0,31,1)
aRDMNAO=seq(0,31,1)
aRDMNARSK=seq(0,31,1)
aRDMNASM=seq(0,31,1)
aRDMNATTH=seq(0,31,1)
aRIR_TR=seq(0,31,1)
aRIS_TR=seq(0,31,1)
aRlyG_V=seq(0,31,1)
aRlyM_R=seq(0,31,1)
aRlyM_V=seq(0,31,1)
aRlyS_V=seq(0,31,1)
aRlyW_V=seq(0,31,1)
ARRTOUR=seq(0,31,1)
aRSIC=seq(0,31,1)
aRTVA_TR=seq(0,31,1)
aSSNF_Flux=seq(0,31,1)
aTRFM_A=seq(0,31,1)
aUIAA=seq(0,31,1)
aUIAIA=seq(0,31,1)
aUIANA=seq(0,31,1)
aUIANAHIA=seq(0,31,1)
aUIA_SF=seq(0,31,1)
aUIGA=seq(0,31,1)
aUIGIA=seq(0,31,1)
aUIGNAHIAHIC=seq(0,31,1)
aUIICNAHIAHIC=seq(0,31,1)
aUIG_SF=seq(0,31,1)
aUIMNAA=seq(0,31,1)
aUIMNAIA=seq(0,31,1)
aUIMNANAHIA=seq(0,31,1)
aUIMNAHIC_SF=seq(0,31,1)
aUINAA=seq(0,31,1)
aUINANA=seq(0,31,1)
aCG_SF=seq(0,31,1)
AURF_TR=seq(0,31,1)
BgBAMdot=seq(0,31,1)
BshicWdot=seq(0,31,1)
CGA=seq(0,31,1)
CGG=seq(0,31,1)
CGMNA=seq(0,31,1)
choc=seq(0,31,1)
CwW_SHICdot=seq(0,31,1)
DivBAM=seq(0,31,1)
EXM_BPCN=seq(0,31,1)
HarmBAM=seq(0,31,1)
HarmG=seq(0,31,1)
HarmIC=seq(0,31,1)
HarmM=seq(0,31,1)
HarmW=seq(0,31,1)
IPCUE=seq(0,31,1)
IPCUSA=seq(0,31,1)
IPDBAM=seq(0,31,1)
IPPRODA=seq(0,31,1)
IPPRODA_Infl=seq(0,31,1)
IPS_Flux=seq(0,31,1)
IPSA_Flux=seq(0,31,1)
IPSNA_Flux=seq(0,31,1)
lambdaDBG_IC=seq(0,31,1)
lambdaOG=seq(0,31,1)
lambdaOIC=seq(0,31,1)
lambdaOSHIC=seq(0,31,1)
lambdaOW=seq(0,31,1)
LG=seq(0,31,1)
MuA=seq(0,31,1)
MuG=seq(0,31,1)
MuIC=seq(0,31,1)
MuMNA=seq(0,31,1)
ORDTSdot=seq(0,31,1)
pA_SIF=seq(0,31,1)
PACT=seq(0,31,1)
Pble_dolr=seq(0,31,1)
Pble_dolr_2007=seq(0,31,1)
pCIG_SIF=seq(0,31,1)
PFATBAM=seq(0,31,1)
phiAoR=seq(0,31,1)
POP=seq(0,31,1)
POPR=seq(0,31,1)
POPU=seq(0,31,1)
PPET=seq(0,31,1)
PPET_2007=seq(0,31,1)
PRODA=seq(0,31,1)
PRODBAM_NV=seq(0,31,1)
PSG_STR=seq(0,31,1)
PSIC_STR=seq(0,31,1)
PSSHIC_STR=seq(0,31,1)
PSW_STR=seq(0,31,1)
PXAPH_dolr=seq(0,31,1)
PXAPH_dolr_2007=seq(0,31,1)
PXENG_dolr=seq(0,31,1)
PXENG_dolr_2007=seq(0,31,1)
PXPH_dolr=seq(0,31,1)
PXPH_dolr_2007=seq(0,31,1)
RAID_TR=seq(0,31,1)
RassIC_V=seq(0,31,1)
RCST_TR=seq(0,31,1)
RETR=seq(0,31,1)
RNF_TR=seq(0,31,1)
RSBAM_B=seq(0,31,1)
RSE_RHMRE=seq(0,31,1)
RSE_V=seq(0,31,1)
RSG=seq(0,31,1)
RTIC_TR=seq(0,31,1)
S_IDE_E=seq(0,31,1)
S_IDE_S=seq(0,31,1)
SCST=seq(0,31,1)
SRPM=seq(0,31,1)
TAIPTIA=seq(0,31,1)
TAIPTIA0=seq(0,31,1)
TAIPTNAHIA=seq(0,31,1)
TAIPTNAHIA0=seq(0,31,1)
TAIPTIC=seq(0,31,1)
TauH=seq(0,31,1)
tAutreIC=seq(0,31,1)
tAutreM=seq(0,31,1)
tAutreSHIC=seq(0,31,1)
tAutreW=seq(0,31,1)
TCDOLR=seq(0,31,1)
TCDOLR_2007=seq(0,31,1)
TCEURO=seq(0,31,1)
TCEURO_2007=seq(0,31,1)
tCSSG=seq(0,31,1)
tCSSIC=seq(0,31,1)
tCSSSHIC=seq(0,31,1)
tCSSW=seq(0,31,1)
TDDA=seq(0,31,1)
TDDIA=seq(0,31,1)
TDDNAHIA=seq(0,31,1)
TDEC=seq(0,31,1)
tDivG=seq(0,31,1)
tDivIC=seq(0,31,1)
tDivM=seq(0,31,1)
tDivSHIC=seq(0,31,1)
tDivW=seq(0,31,1)
thetaFXCw_G=seq(0,31,1)
thetaFXDetteG=seq(0,31,1)
TintAOR=seq(0,31,1)
TintBAM=seq(0,31,1)
TintBw=seq(0,31,1)
TintD_IMB=seq(0,31,1)
TintRfw_cd=seq(0,31,1)
TIPDA=seq(0,31,1)
TIPDG=seq(0,31,1)
TIPDIC=seq(0,31,1)
TIPDMNAHIC=seq(0,31,1)
tLNRMNA=seq(0,31,1)
tLRNSMNA=seq(0,31,1)
tLSAG=seq(0,31,1)
tLSMNA=seq(0,31,1)
tMBIA=seq(0,31,1)
TPSM=seq(0,31,1)
TRF_MRE=seq(0,31,1)
TRFBAM=seq(0,31,1)
TSG_B=seq(0,31,1)
TSPDIC=seq(0,31,1)
TSPDMNAHIC=seq(0,31,1)
TSPTA=seq(0,31,1)
TSPTIA=seq(0,31,1)
TSPTIA0=seq(0,31,1)
TSPTNAHIA=seq(0,31,1)
TSPTNAHIA0=seq(0,31,1)
TTICIA=seq(0,31,1)
TTICIA0=seq(0,31,1)
TTICNAHIA=seq(0,31,1)
TTICNAHIA0=seq(0,31,1)
tTNSMNA=seq(0,31,1)
tTRFIC=seq(0,31,1)
tTRFSHIC=seq(0,31,1)
TRFW_A=seq(0,31,1)
TTRPBAM=seq(0,31,1)
TTRPIC=seq(0,31,1)
TTRPMNA=seq(0,31,1)
TTRPSHIC=seq(0,31,1)
TTVAIA=seq(0,31,1)
TTVAIA0=seq(0,31,1)
TTVANAF=seq(0,31,1)
TTVANAHIAC=seq(0,31,1)
TTVANAHIAC0=seq(0,31,1)
TTVAIC=seq(0,31,1)
UIBAMNAHIA=seq(0,31,1)
UIICNAHIA=seq(0,31,1)
VDMHOCP=seq(0,31,1)
VIG=seq(0,31,1)
VPRODA=seq(0,31,1)
VSMIG=seq(0,31,1)
VTSAG_B=seq(0,31,1)
VXOCP=seq(0,31,1)
XSVpartouriste=seq(0,31,1)
lambdaDwW_IC=seq(0,31,1)
lambdaCwIC_W=seq(0,31,1)
lambdaBgSHIC=seq(0,31,1)
VPIB_observe=seq(0,31,1)
lambdaMFIDM=seq(0,31,1)
VPRODA2=seq(0,31,1)
PRODA2=seq(0,31,1)
VISOEI_LEAP=seq(0,31,1)

#===========================================================
#===========================================================
##time derivatives
#===1/- Production et produit ==================
VYDMNAe=betaVYDMNAe*(VYDMNANet-VYDMNAe)+VYDMNAe*grVYDMNAe
grVYDMNAe=betagrVYDMNAe*(VINA/VSKNA-TDEC-grVYDMNAe)

#===2/- PIB en volume==================
#==2-1/ chainage VPIB============
omegaIPC=betaomegaIPC*(IPCF/IPPIB-omegaIPC)
omegaIPI=betaomegaIPI*(IPI/IPPIB-omegaIPI)
omegaIPCG=betaomegaIPCG*(IPCG/IPPIB-omegaIPCG)
omegaIPX=betaomegaIPX*(IPX/IPPIB-omegaIPX)
omegaIPM=betaomegaIPM*(IPM/IPPIB-omegaIPM)
omegaIPS_flux=betaomegaIPS_flux*(IPS_Flux/IPPIB-omegaIPS_flux)
#==2-2/ chainage VVAMNA============
omegaIPVANA=betaomegaIPVANA*(IPVAMNA/IPPRODMNA-omegaIPVANA)
omegaIPUIAMNA=betaomegaIPUIAMNA*(IPUIA/IPPRODMNA-omegaIPUIAMNA)
omegaIPUINAMNA=betaomegaIPUINAMNA*(IPUINA/IPPRODMNA-omegaIPUINAMNA)
#==2-3/ chainage VVAAG============
omegaIPVAAG=betaomegaIPVAAG*(IPVAAG/IPPRODA-omegaIPVAAG)
omegaIPUIAA=betaomegaIPUIAA*(IPUIA/IPPRODA-omegaIPUIAA)
omegaIPUINAA=betaomegaIPUINAA*(IPUINA/IPPRODA-omegaIPUINAA)
#==2-4/ chainage VVAG============
omegaIPVAG=betaomegaIPVAG*(IPVAG/IPPRODG-omegaIPVAG)
omegaIPUIAG=betaomegaIPUIAG*(IPUIA/IPPRODG-omegaIPUIAG)
omegaIPUINAG=betaomegaIPUINAG*(IPUINA/IPPRODG-omegaIPUINAG)
#==2-5/ chainage VVAT============
omegaIPVAAGT=betaomegaIPVAAGT*(IPVAAG/IPVAT-omegaIPVAAGT)
omegaIPVAMNAT=betaomegaIPVAMNAT*(IPVAMNA/IPVAT-omegaIPVAMNAT)
omegaIPVAGT=betaomegaIPVAGT*(IPVAG/IPVAT-omegaIPVAGT)

#===3/- consommation des menages en volume===============
#===3-1/ chainage VCONS============
omegaIPCA=betaomegaIPCA*(IPCFA/IPCF-omegaIPCA)
omegaIPCIA=betaomegaIPCIA*(IPCFIA/IPCF-omegaIPCIA)
omegaIPCNAHIA=betaomegaIPCNAHIA*(IPCFNAHIA/IPCF-omegaIPCNAHIA)
omegaIPCIANA=betaomegaIPCIANA*(IPCFIA/IPCFNA-omegaIPCIANA)
omegaIPCNAHIANA=betaomegaIPCNAHIANA*(IPCFNAHIA/IPCFNA-omegaIPCNAHIANA)

#===4/- consommation APU en volume===============
#===chainage VCG============
omegaIPCGG=betaomegaIPCGG*(IPCGG/IPCG-omegaIPCGG)
omegaIPCGM=betaomegaIPCGM*(IPCGM/IPCG-omegaIPCGM)
#===5/- FBCF en volume===============
VISOEILT=betaVISOEI*(VISOEIC-VISOEI)*(1+kappa10*DUM0715)
VILOGLT=betaVILOG*(VILOGC-VILOG)
ratioDetteRevenu=betaRatioDette*(ratioDetteRevenuC-ratioDetteRevenu)
#==chainage VITOT============
omegaIPINA=betaomegaIPINA*(IPINA/IPI-omegaIPINA)
omegaIPIA=betaomegaIPIA*(IPIA/IPI-omegaIPIA)

#====6/- Stock de Capital====================
VSKNA = VINA - TDEC*VSKNA

#===7/- Stocks en volume==================
VSNA=(PRODMNA-PFATMNA+MTCNA+DDNA+AIPTNA-SPTNA+TICNA+TVANA+MNA-DEMNA)/IPSNA_Flux 

# ===8/- commerce exterieur en volume===========================
#===Import en volume============================
sigmaMBA = betasigmaMBA*(sigmaMBAC-sigmaMBA)
#==Chainage VMBHEA============
sigmaMC=betaVMBHEAMC*(sigmaMCC-sigmaMC)
sigmaMI=betaVMBHEAMI*(sigmaMIC-sigmaMI)
sigmaMX=betaVMBHEAMX*(sigmaMXC-sigmaMX)
sigmaMUI=betaVMBHEAMUI*(sigmaMUIC-sigmaMUI)
sigmaME=betaVMEMC*(sigmaMEC-sigmaME)
#==chainage VMBS============
omegaIPMBA=betaomegaIPMBA*(IPMBA/IPM-omegaIPMBA)
omegaIPMHEA=betaomegaIPMHEA*(IPMHEA/IPM-omegaIPMHEA)
omegaIPME=betaomegaIPME*(IPME/IPM-omegaIPME)
omegaIPMS=betaomegaIPMS*(IPMS/IPM-omegaIPMS)

#===export en volume============================
#==elasticite competitivite prix dans VXBHAOCP
sigmaX=betaVXBHAOCP*(sigmaXC-sigmaX)
#==export agricole et demande en volume
VXBA= betaVXBA*(VXBAC-VXBA)
VDEMAe = betaVDEMAe*(VDEMAeC-VDEMAe)
#==export services en volume
VXSTR=betaVXSTR*(VXSTRC-VXSTR)
VXSHVTR=betaVXHVSTR*(VXSHVTRC-VXSHVTR)
#==chainage VXBS============
omegaIPXBA=betaomegaIPXBA*(IPXBA/IPX-omegaIPXBA)
omegaIPXHAOCP=betaomegaIPXHAOCP*(IPXHAOCP/IPX-omegaIPXHAOCP)
omegaIPXOCP=betaomegaIPXOCP*(IPXOCP/IPX-omegaIPXOCP)
omegaIPXSHVTR=betaomegaIPXSHVTR*(IPXSHVTR/IPX-omegaIPXSHVTR)
omegaIPXSTR=betaomegaIPXSTR*(IPXSTR/IPX-omegaIPXSTR)
omegaIPXSV=betaomegaIPXSV*(IPXSV/IPX-omegaIPXSV)

#===9/- Productivite ===================
# rhoAG=(rhoAG1+rhoAG2*DUM1719+rhoAG3*DUM09)*(LAG/POPR-rhoAG4+rhoAG5*DUM1719+rhoAG6*DUM09)*rhoAG
rhoAG=croissanceExogRhoAG*rhoAG
# rhoMNA=DUM1950*0.02*rhoMNA+(1-DUM1950)*(rhoMNA1+rhoMNA2*(VISOEI/VSKNA-TDEC-rhoMNA3))*rhoMNA*(1-rhoMNA5*DUM0712)
# rhoMNALT=rhoMNA_Cr*rhoMNALT
rhoMNALT=croissanceExogRhoMNA*rhoMNALT

#===10/- Taux de salaire===================================
TSAG_B = betaTSAGB*(TSAG_BC-TSAG_B)
VTSMNA_B=(omega0+omega1*croissanceExogRhoMNA+omega2*LMNA/POP)*VTSMNA_B
#VTSMNA_B=(omega0+omega1*rhoMNA_Cr+omega2*LMNA/POP)*VTSMNA_B

#===11/- Prix production =============================
IPPRODMNALT = IPPRODMNALT_Infl 
muVAMNA = betaMUVANA*(muVAMNAC-muVAMNA)

#===12/- Prix de valeur ajoutee=============================

#===13/- Prix des utilisations intermdiaires=============================
IPUIALT=phiIPUIA0*IPPRODA_Infl*(1+phiIPUIA3*DUM0716)+phiIPUIA1*IPME_Infl+phiIPUIA2*(1-phiIPUIA0-phiIPUIA1)*IPMBA_Infl*(1+phiIPUIA4*DUM0716)
IPUINALT=phiIPUINA0*(1+phiIPUINA3*DUM0715)*IPPRODMNALT_Infl+phiIPUINA1*IPME_Infl+phiIPUINA2*(1-phiIPUINA0-phiIPUINA1)*IPMHEA_Infl
IPUIANALT=phiIPUIANA0*IPPRODMNALT_Infl*(1+phiIPUIANA3*DUM0715)+phiIPUIANA1*IPME_Infl+phiIPUIANA2*(1-phiIPUIANA0-phiIPUIANA1)*IPMHEA_Infl
IPUINANALT=phiIPUINANA0*(1+phiIPUINANA3*DUM0715)*IPPRODMNALT_Infl+phiIPUINANA1*IPME_Infl+phiIPUINANA2*(1-phiIPUINANA0-phiIPUINANA1)*IPMHEA_Infl

#===14/- Prix à la consommation des menages=============================
IPCFALT = betaIPCFA*(IPCFAC-IPCFALT) 
IPCFIALT=phiIPCIA0*IPPRODMNALT_Infl*(1+phiIPCIA3*DUM0715)+phiIPCIA1*IPME_Infl+phiIPCIA2*(1-phiIPCIA0-phiIPCIA1)*IPMHEA_Infl
IPCFNAHIALT=phiIPCNAHIA0*IPPRODMNALT_Infl*(1+phiIPCNAHIA3*DUM0715)*(1+phiIPCNAHIA4*DUM0709)+phiIPCNAHIA1*IPME_Infl+phiIPCNAHIA2*(1-phiIPCNAHIA0-phiIPCNAHIA1)*IPMHEA_Infl
#IPCFLag=betaIPCFLag*(IPCF-IPCFLag)
#inflLT=betaINFLT*(inflLT-infl)

#===15/- Prix à la consommation des APU=============================
IPCGG = betaIPCGG*(IPCGGC-IPCGG) 
# IPCGMLT=phiIPCGM0*(1+phiIPCGM3*DUM0715)*IPPRODMNALT_Infl+phiIPCGM1*IPME_Infl+phiIPCGM2*(1-phiIPCGM0-phiIPCGM1)*IPMHEA_Infl
IPCGMLT=0.015*IPCGM #TODOANTOINE
#===16/- Prix de la FBCF=============================
IPILT=phiIPI0*(1+phiIPI3*DUM0715)*IPPRODMNALT_Infl+phiIPI1*IPME_Infl+phiIPI2*(1-phiIPI0-phiIPI1)*IPMHEA_Infl
IPIALT = betaIPIA*(IPIAC-IPIALT) 
# inflIPILT=betaINFLT*(inflIPILT-inflIPI)
#IPILag=betaIPILag*(IPI-IPILag)

#===17/- Prix commerce exterieur===============================
#===Prix import================================
IPMBALT=IPMBA_Infl
IPMELT=IPME_Infl
IPMHEALT=IPMHEA_Infl
IPMSLT=betaIPMS*(IPMSC-IPMSLT)

#===Prix export================================
IPXBA = betaIPXBA*(IPXBAC-IPXBA) 
IPXHAOCP=betaIPXHAOCP*(IPXHAOCPC-IPXHAOCP)
IPXSHVTR=betaIPXSHVTR*(IPXSHVTRC-IPXSHVTR)
IPXSTR=betaIPXSTR*(IPXSTRC-IPXSTR)
IPXSV=betaIPXSV*(IPXSVC-IPXSV)

#===18/- Prix etrangers et taux de change================================
IPCET=betaIPCET*(IPCETC-IPCET)
TCEN=betaTCEN*(TCENC-TCEN)

#===19/- Interets====================
#CMRBQ=betaCMRBQ*(CMRBQC-CMRBQ)

#===20/- consommation des menages en nominal===============
CONS=betaC*(CONST-CONS)

#===21/- Encours Credit et Depot================================
# Menages
#=========
#Numeraire Depots opcvm non monetaires
MFIDM=MFIDMdot
DBM_IC=DBM_ICdot
#Credits
CBIC_M=CBIC_Mdot
#OM=OMdot
#==============
# APU
#==============
# dette APU: totale, domestique et exterieure
DetteDG=DetteDGdot
FXDetteG=betaFXDetteG*(FXDetteGC-FXDetteG)
#thetaFXDetteG=betathetaFXDetteG*(0.4-thetaFXDetteG)
#phiAoR=betaphiAoR*(7/12-phiAoR)
#Numeraire et Depots et opcvm non monetaires
DBG_IC=DBG_ICdot
Dg_N=Dg_Ndot
#Titres autres qu actions
Bg_N=Bg_Ndot
BgW=BgWdot
#taux d interet sur eurobonds
#Credits
Cw_G=Cw_Gdot
FXCw_G=FXCw_Gdot
CBIC_G=CBIC_Gdot
#Actions titres de participation et opcvm non monetaires
# E_BAM= E_BAMdot
EG_IC= EG_ICdot
EG_SHIC= EG_SHICdot
#================
# Reste du Monde
#================
#==Avoirs de reserves etrangères
AR=ARdot
#==Or monetaire et DTS
ORDTS=ORDTSdot
#==Numeraire et Depots et opcvm non monetaires
Dw_W=Dw_Wdot
DwW_IC = DwW_ICdot
#lambdaDwW_IC=betalambdaDwW_IC*(0.03-lambdaDwW_IC)
#==Titres autres qu actions
BshicgW=BshicgWdot
Bw=Bwdot
BfW_SHIC=BfW_SHICdot
#lambdaBfW_SHIC=betalambdaBfW_SHIC*(0.40-lambdaBfW_SHIC)
#==Credits
CwIC_W=CwIC_Wdot
CwW=CwWdot
#==Actions titres de participation et opcvm non monetaires
#IDE_N=IDE_Ndot
EW_N=EW_Ndot
EW_IC_N=EW_IC_Ndot
EW_SHIC_N=EW_SHIC_Ndot
#=================================================
# Societes et entreprises Hors Instituts de credits
#=================================================
#==Numeraire et Depots et opcvm non monetaires
DBSHIC_IC=DBSHIC_ICdot
DgSHIC=DgSHICdot
DwSHIC_IC=DwSHIC_ICdot
#==Titres autres qu actions
BgSHIC=BgSHICdot
BshicW=BshicWdot
#BfSHIC_IC=BfSHIC_ICdot
#Bf_SHIC_N=Bf_SHIC_Ndot
#==Credits
LNDSHIC_N=LNDSHIC_Ndot
CwIC_SHIC=CwIC_SHICdot
CwW_SHIC=CwW_SHICdot
#==Actions titres de participation et opcvm non monetaires
E_SHIC_N=E_SHIC_Ndot
ESHIC_IC=ESHIC_ICdot
#======================
# Instituts de Credits
#======================
#==Numeraire et Depots et opcvm non monetaires
DBIC_BAM=DBIC_BAMdot
DgIC=DgICdot
DB_IC_N=DB_IC_Ndot
DwIC_W=DwIC_Wdot
Dw_IC_N=Dw_IC_Ndot
#==Titres autres qu actions
BgIC=BgICdot
BwIC=BwICdot
#BfIC_SHIC=BfIC_SHICdot
#==Credits
# CBBAM_IC=CBBAM_ICdot
#CBIC_N=CBIC_Ndot
CwIC_N=CwIC_Ndot
LNDIC_N=LNDIC_Ndot
CwW_IC=CwW_ICdot
#==Actions titres de participation et opcvm non monetaires
EIC_SHIC=EIC_SHICdot
E_IC_N=E_IC_Ndot
#================
# Banque centrale
#================
#==Avoirs officiels de reserves
BwBAM=BwBAMdot
DwBAM_W=DwBAM_Wdot
#==Numeraire et Depots et opcvm non monetaires
DB_BAM=DB_BAMdot
#==Titres autres qu actions
BgBAM=BgBAMdot
#==Credits
CBBAM_IC=CBBAM_ICdot 
#==Actions titres de participation et opcvm non monetaires
E_BAM=E_BAMdot

#================22/-Taux d'interets===========================================
#======Taux d'interets domestiques 
TintCN = betaTintC*(TintCC-TintCN)
TintDN = betaTintD*(TintDC-TintD)
TintBgN = betaTintBg*(TintBgC-TintBgN)
TintBN = betaTintB*(TintBC-TintBN)

#======Taux d'interets extérieurs
TintBgWN = betaTintBgW*(TintBgWC-TintBgWN)
TintBshicWN = betaTintBshicW*(TintBshicWC-TintBshicWN)

TintCDgWN = betaTintCDgW*(TintCDgWC-TintCDgWN)
TintCDshicWN = betaTintCDshicW*(TintCDshicWC-TintCDshicWN)


##initial values
VYDMNAe = 1129067.64186305
grVYDMNAe = 0.00768999245869826
omegaIPC = 0.973832391580934
omegaIPI = 0.991141758941049
omegaIPCG = 1.06928473781957
omegaIPX = 1.0252575440016
omegaIPM = 0.950166882742893
omegaIPS_flux = 0.473024465557037
omegaIPVANA = 1.00709918815487
omegaIPUIAMNA = 1.05052766777703
omegaIPUINAMNA = 0.983609444134001
omegaIPVAAG = 0.939471111837234
omegaIPUIAA = 1.04461768141018
omegaIPUINAA = 0.978075921711465
omegaIPVAG = 1.05029375769585
omegaIPUIAG = 0.978176286440345
omegaIPUINAG = 0.915866819011626
omegaIPVAAGT = 0.924972936495649
omegaIPVAMNAT = 0.985979122579731
omegaIPVAGT = 1.10432428391406
omegaIPCA = 1.06323732485921
omegaIPCIA = 1.09183744359706
omegaIPCNAHIA = 0.946583896393538
omegaIPCIANA = 1.16727358582358
omegaIPCNAHIANA = 0.955871189819095
omegaIPCGG = 1.00369715314859
omegaIPCGM = 0.901883803341154
VISOEILT = 172655.265073895
VILOGLT = 55524.4520319973
ratioDetteRevenu = 0.574521091442765
omegaIPINA = 1.0018010834487
omegaIPIA = 0.963487410285539
VSKNA = 3619783.02861867
VSNA = 143086.947325813
sigmaMBA = 0.090090591318247
sigmaMC = 0.287868376482504
sigmaMI = 0.432942107486661
sigmaMX = 0.279489571800929
sigmaMUI = 0.0265610412992363
sigmaME = 0.999883519411653
omegaIPMBA = 1.24057867167277
omegaIPMHEA = 1.02099850766489
omegaIPME = 0.943735660153048
omegaIPMS = 0.993024223248491
sigmaX = 0.0879638962886691
VXBA = 15504.0184543355
VDEMAe = 195145.077801138
VXSTR = 28690.9630574331
VXSHVTR = 49581.9177895526
omegaIPXBA = 1.50249917684989
omegaIPXHAOCP = 1.01928585356742
omegaIPXOCP = 1.16501141042718
omegaIPXSHVTR = 0.860275971703063
omegaIPXSTR = 0.95733973235151
omegaIPXSV = 1.00014196063114
rhoAG = 47.4040051630268
rhoMNA = 179.802542262605
TSAG_B = 15.4440702707087
VTSMNA_B = 35.7109608332344
IPPRODMNALT = 1.49928564919425
muVAMNA = 0.903
IPUIALT = 1.20207451059883
IPUINALT = 1.12550280915469
IPUIANALT = 1.38403342614422
IPUINANALT = 1.11296849475054
IPCFALT = 1.18558481354524
IPCFIALT = 1.21747596865087
IPCFNAHIALT = 1.05550799061655
IPCGG = 1.22889353101501
IPCGMLT = 1.10423663968396
IPILT = 1.13489052964258
IPIALT = 1.09345273736291
IPMBALT = 0.59450869085899
IPMELT = 0.776989803207688
IPMHEALT = 1.11081873099198
IPMSLT = 1.08038346700036
IPXBA = 1.76386528954064
IPXHAOCP = 1.19659495654219
IPXSHVTR = 1.00992463043757
IPXSTR = 1.12387304446534
IPXSV = 1.17412090212854
IPCET = 1.24080663982488
TCEN = 1.08687877773595
CONS = 660652
MFIDM = 233602.01476481
DBM_IC = 721365.644365402
CBIC_M = 357212.510252179
DetteDG = 439834.628140921
FXDetteG = 136191.709721617
DBG_IC = 89556.6207936653
Dg_N = 47394.32800297
Bg_N = 360348.35585494
BgW = 42512
Cw_G = 105511.879
FXCw_G = 97077.8721246074
CBIC_G = 79486.2722859805
EG_IC = 26459.907
EG_SHIC = 304124.70646428
AR = 246247.52985214
ORDTS = 15936.09111912
Dw_W = 64351.28273302
DwW_IC = 23244.5130224652
BshicgW = 74509.099
Bw = 165960.156
BfW_SHIC = 27229
CwIC_W = 22547.352
CwW = 330511.584701505
EW_N = 512114.039
EW_IC_N = 8309.55344250793
EW_SHIC_N = 503804.485557492
DBSHIC_IC = 238898.525385809
DgSHIC = 41807.9262116057
DwSHIC_IC = 16403.7079228251
BgSHIC = 128921.636904117
BshicW = 31997.099
LNDSHIC_N = 521833.071645949
CwIC_SHIC = 15250.449
CwW_SHIC = 208803.835334205
E_SHIC_N = 857332.165784179
ESHIC_IC = 201685.664953537
DBIC_BAM = 42189.49985677
DgIC = 5586.40179136432
DB_IC_N = 1049820.79054488
DwIC_W = 17439.168
Dw_IC_N = 39648.2209452903
BgIC = 230968.583390003
BwIC = 5662.09999999998
CwIC_N = 37797.801
LNDIC_N = 931302.854184109
CwW_IC = 16195.8703673
EIC_SHIC = 49402.9737624068
E_IC_N = 236455.125396045
BwBAM = 160298.056
DwBAM_W = 46912.11473302
DB_BAM = 275791.51462158
BgBAM = 458.13556082
CBBAM_IC = 73305.23299775
E_BAM = 5528.38653572
TintCN = 1.33
TintDN = 5
TintBgN = 4.34796139538971
TintBN = 2.27769230769231
TintBgWN = 4.32113285660519
TintBshicWN = 5.02545558895824
TintCDgWN = 1.89931220919684
TintCDshicWN = 2.11068441339214
##parameters
alpha10 = 0.8
alpha20 = 0.8
alpha11 = 0.1
alpha21 = 0.1
alpha3 = 0.5
alpha4 = 0.5
alpha5 = 0
betaC = 4
epsilon1A = 0.5
epsilon2A = 0.5
epsilon3A = 0.5
epsilon1IA = 0.5
epsilon2IA = 0.5
epsilon3IA = 0.5
consamin = 5
consiamin = 5
consnahiamin = 5
alpha8 = 0
alpha9 = 0
betaomegaIPCA = 4
betaomegaIPCIA = 4
betaomegaIPCNAHIA = 4
kappa1 = 0.5
kappa2 = 0.05
kappa3 = 0.05
kappa10 = 1
kappa11 = 0
kappa12 = 0
betaVISOEI = 5
betaRatioDette = 5
kappaVILOG1 = 0.5
kappaVILOG2 = 0.05
kappaVILOG3 = 0.05
kappaVILOG4 = 0.5
kappaVILOG5 = 1
kappaVILOG6 = 1
betaVILOG = 5
tVIA = 0.05
tVIA1 = 1
betaomegaIPINA = 5
betaomegaIPIA = 5
rhoAG1 = 10
rhoAG2 = 1
rhoAG3 = 0.5
rhoAG4 = 0.5
rhoAG5 = 10
rhoAG6 = 0
rhoAG7 = 0
rhoAG8 = 0
rhoMNA1 = 1
rhoMNA2 = 0
rhoMNA3 = 0
rhoMNA4 = 0
rhoMNA5 = 0
rhoMNA6 = 0
omega0 = 1
omega1 = 1
omega2 = 1
omega3 = 1
omega4 = 1
omega5 = 1
betaTSAGB = 5
muVANA0 = 1
muVANA1 = 1
muVANA2 = 1
muVANA3 = 0
muVANA4 = 0
muVANA5 = 0
betaIPVANA = 5
betaMUVANA = 5
phiIPVANA0 = 0.5
phiIPVANA1 = 0.5
phiIPVANA2 = 0.5
phiIPVANA3 = 0.5
betaomegaIPVANA = 5
betaomegaIPUIAMNA = 5
betaomegaIPUINAMNA = 5
phiIPCIA0 = 0.5
phiIPCIA1 = 0.5
phiIPCIA2 = 0.5
phiIPCIA3 = 0
phiIPI0 = 0.5
phiIPI1 = 0.5
phiIPI2 = 0.5
phiIPI3 = 0
phiIPUIA0 = 0.5
phiIPUIA1 = 0
phiIPUIA2 = 0
phiIPUIA3 = 0
phiIPUIA4 = 0
phiIPCGM0 = 0.5
phiIPCGM1 = 0.5
phiIPCGM2 = 0.5
phiIPCGM3 = 0
phiIPUIANA0 = 0.5
phiIPUIANA1 = 0.5
phiIPUIANA2 = 0.5
phiIPUIANA3 = 0
phiIPUIANA4 = 0
phiIPUINANA0 = 0.5
phiIPUINANA1 = 0.5
phiIPUINANA2 = 0.5
phiIPUINANA3 = 0
phiIPCNAHIA0 = 0.5
phiIPCNAHIA1 = 0.5
phiIPCNAHIA2 = 0.5
phiIPCNAHIA3 = 0
phiIPCNAHIA4 = 0
phiIPUINA0 = 0.5
phiIPUINA1 = 0.5
phiIPUINA2 = 0.5
phiIPUINA3 = 0
alphaIPMBA0 = 0.5
alphaIPMBA1 = 0.5
alphaIPMBA2 = 0
alphaIPMBA3 = 0
alphaIPMBA4 = 0
alphaIPMBA5 = 0
alphaIPME0 = 0.5
alphaIPME1 = 0.5
alphaIPME2 = 0
alphaIPME3 = 0
alphaIPMHEA0 = 0.5
alphaIPMHEA1 = 0
alphaIPMHEA2 = 0
alphaIPMS = 0.5
alphaIPMS1 = 0
alphaIPMS2 = 0
alphaIPXOCP1 = 0.5
alphaIPXOCP2 = 0.5
alphaIPXOCP3 = 0.5
alphaIPXHAOCP1 = 0.5
alphaIPXHAOCP2 = 0.5
alphaIPXSHVTR1 = 0.5
alphaIPXSHVTR2 = 0.5
alphaIPXSHVTR3 = 0
alphaIPXSHVTR4 = 0
alphaIPXSTR0 = 0.5
alphaIPXSTR1 = 0.5
alphaIPXSTR2 = 0
alphaIPXSV0 = 0.5
alphaIPXSV1 = 0.5
alphaIPXSV2 = 0
alphaIPCET0 = 1
alphaIPCET1 = 0.3
alphaIPCET2 = 0.2
alphaTCEN0 = 0.5
alphaTCEN1 = 0.5
alphaTCEN2 = 0.5
alphaTCEN3 = 0
betaIPMBA = 4
betaIPME = 4
betaIPMHEA = 4
betaIPMS = 4
betaIPXHAOCP = 4
betaIPXSHVTR = 4
betaIPXSTR = 4
betaIPXSV = 4
betaIPCET = 4
betaTCEN = 4
iota1 = 2
iota2 = 0.5
iota3 = 2
iota4 = 0.5
iota5 = 2
iota6 = 0.5
iota7 = 2
iota8 = 0.5
iota9 = 0
iota10 = 0
iota11 = 0
iota12 = 0
betaVMBHEAMC = 5
betaVMBHEAMI = 5
betaVMBHEAMX = 5
betaVMBHEAMUI = 5
sigmaMS1 = 0.5
sigmaMS2 = 0.5
sigmaMS3 = 0.5
sigmaMS4 = 0
sigmaMS5 = 0
iotaME1 = 0.5
iotaME2 = 0.5
iotaME3 = 0.5
iotaME4 = 0
iotaME5 = 0
iotaME6 = 0
betaVMEMC = 5
betaMFIDM = 1
betaCBIC_M = 1
rhoCBIC_M = 2
rhoCBIC_M2 = 1
rhoCBIC_M3 = 1
betaFXDetteG = 1
betaDBG_IC = 1
betaBgW = 1
betaCw_G = 1
betaFXCw_G = 1
betaCBIC_G = 1
thetaCBIC_G = 1
lambdaCBIC_G1 = 0
lambdaCBIC_G2 = 0.5
lambdaCBIC_G3 = 1
lambdaEBAM = 0.02
betaEBAM = 2
lambdaEGIC = 0.02
betaEGIC = 2
lambdaEGSHIC = 0.02
betaEGSHIC = 2
betaDwW_IC = 2
betaBw = 2
gammaBw1 = 0.5
gammaBw2 = 10
gammaBw3 = 0
gammaBw4 = 0
gammaBw5 = 0
betaBfW_SHIC = 2
lambdaBfW_SHIC = 0.5
lambdaBfW_SHIC1 = 0
lambdaBfW_SHIC2 = 0
lambdaBfW_SHIC3 = 0
betaCwIC_W = 2
lambdaEWIC = 0.02
betaEWIC = 2
lambdaEWSHIC = 0.02
betaEWSHIC = 2
betaDBSHIC_IC = 2
lambdaDBSHIC_IC = 0.5
lambdaDBSHIC_IC1 = 0
lambdaDBSHIC_IC2 = 0
betaDgSHIC = 2
lambdaDgSHIC = 1
lambdaDgSHIC1 = 1
lambdaDgSHIC2 = 1.5
betaDwSHIC_IC = 2
lambdaDwSHIC_IC = 1
lambdaDwSHIC_IC1 = 2.5
lambdaDwSHIC_IC2 = 0
lambdaDwSHIC_IC3 = 0
betaBgSHIC = 2
betaCwIC_SHIC = 2
lambdaCwIC_SHIC = 0.5
lambdaCwIC_SHIC1 = 2
lambdaCwIC_SHIC2 = 0.5
lambdaCwIC_SHIC3 = 0.5
lambdaESHICIC = 0.02
betaESHICIC = 2
betaDBIC_BAM = 1
lambdaDBIC_BAM = 0.1
lambdaDBIC_BAM1 = 0.1
lambdaDBIC_BAM2 = 1
betaDgIC = 2
lambdaDgIC = 0.01
lambdaDgIC1 = 1
lambdaDgIC2 = 0
lambdaDgIC3 = 0
lambdaDgIC4 = 0
betaDwIC_W = 2
thetaFX = 0.5
lambdaDwIC_W1 = 0
lambdaDwIC_W2 = 0.7
betaBwIC = 2
lambdaCwWIC = 0.7
betaCwWIC = 2
lambdaEICSHIC = 0.02
betaEICSHIC = 2
tetamba0 = -5
tetamba1 = 5
tetamba2 = 1
tetamba3 = 1
tetamba4 = 1
tetamba6 = 5
tetamba7 = 0
betasigmaMBA = 10
betaVDEMAe = 10
betaVXBA = 10
nuXBHAOCP0 = 4
nuXBHAOCP1 = 4
nuXBHAOCP2 = 50500000
nuXBHAOCP3 = 0.5
nuXBHAOCP4 = 0.5
nuXBHAOCP5 = 0.5
betaVXBHAOCP = 4
nuVXSTR0 = 4
nuVXSTR1 = 4
betaVXSTR = 4
nuVXSHVTR0 = 4
nuVXSHVTR1 = 4
betaVXHVSTR = 4
betaVYDMNAe = 4
betagrVYDMNAe = 4
upsilon2NA = 1
upsilon1NA = 0.1
betaomegaIPC = 4
betaomegaIPI = 4
betaomegaIPCG = 4
betaomegaIPM = 4
betaomegaIPX = 4
betaomegaIPS_flux = 4
betaomegaIPMBA = 4
betaomegaIPME = 4
betaomegaIPMS = 4
betaomegaIPMHEA = 4
betaomegaIPXBA = 4
betaomegaIPXHAOCP = 4
betaomegaIPXOCP = 4
betaomegaIPXSHVTR = 4
betaomegaIPXSTR = 4
betaomegaIPXSV = 4
betaomegaIPCGG = 4
betaomegaIPCGM = 4
betaomegaIPCIANA = 4
betaomegaIPCNAHIANA = 4
betaomegaIPVAAG = 4
betaomegaIPUIAA = 4
betaomegaIPUINAA = 4
betaomegaIPVAG = 4
betaomegaIPUIAG = 4
betaomegaIPUINAG = 4
betaomegaIPVAAGT = 4
betaomegaIPVAMNAT = 4
betaomegaIPVAGT = 4
rhoTintC1 = 0.1
rhoTintC2 = 0.8
rhoTintC3 = 0
rhoTintC4 = 0
rhoTintC5 = 0
betaTintC = 2
xiTintD0 = 2
xiTintD1 = 0.3
xiTintD2 = 1
xiTintD3 = 0
betaTintD = 2
nuTintBg1 = 0.5
nuTintBg2 = 0.6
nuTintBg3 = 0
nuTintBg4 = 0
betaTintBg = 2
alphaTintB1 = 0.01
alphaTintB2 = 0.2
alphaTintB3 = 0
alphaTintB4 = 0
betaTintB = 2
muTintBgW1 = 1
muTintBgW2 = 1
muTintBgW3 = 1
muTintBgW4 = 1
muTintBgW5 = 0
betaTintBgW = 2
muTintBshicW1 = 1
muTintBshicW2 = 1
muTintBshicW3 = 1
muTintBshicW4 = 1
muTintBshicW5 = 0
muTintBshicW6 = 0
betaTintBshicW = 2
muTintCDgW1 = 1
muTintCDgW2 = 1
muTintCDgW3 = 1
muTintCDgW4 = 1
muTintCDgW5 = 1
muTintCDgW6 = 1
betaTintCDgW = 2
muTintCDshicW1 = 1
muTintCDshicW2 = 1
muTintCDshicW3 = 1
muTintCDshicW4 = 1
muTintCDshicW5 = 1
betaTintCDshicW = 2
phiIPCFA0 = 0.5
phiIPCFA1 = 0.5
phiIPCFA2 = 0.1
phiIPCFA3 = 0.5
phiIPCFA4 = 0.5
phiIPCFA5 = 0.5
phiIPCFA6 = 0.5
phiIPCFA7 = 0.5
phiIPCFA8 = 0.5
phiIPCFA9 = 0
betaIPCFA = 10
phiIPIA0 = 0.5
phiIPIA1 = 0.5
phiIPIA2 = 0.1
phiIPIA3 = 0.5
phiIPIA4 = 0.5
phiIPIA5 = 0.5
phiIPIA6 = 0.5
phiIPIA7 = 0.5
phiIPIA8 = 0.5
phiIPIA9 = 0
phiIPIA10 = 0
betaIPIA = 0.5
phiIPXBA0 = 0.5
phiIPXBA1 = 0.5
phiIPXBA2 = 0.1
phiIPXBA3 = 0.5
betaIPXBA = 10
muG0 = 10
muG1 = 10
betaIPCGG = 10
varianteVIG = 0
VARPROD_VIG = 0
varianteTCEN = 0
elastOSHIC = 0
initOSHIC = 0
LBOSHIC = 0
UBOSHIC = 0
elastOIC = 0
initOIC = 0
LBOIC = 0
UBOIC = 0
elastOM = 0
initOM = 0
LBOM = 0
UBOM = 0
elastOG = 0
initOG = 0
LBOG = 0
UBOG = 0
elastOW = 0
initOW = 0
LBOW = 0
UBOW = 0
elastExp = 0
initExp = 0
LBExp = 0
UBExp = 0
croissanceExogRhoMNA = 0
croissanceExogRhoAG = 0
##time
begin= 0 
end= 31 
by= 0.01
