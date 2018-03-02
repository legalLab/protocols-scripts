#! /usr/bin/env python3

#######################################
# A simple script to estimate number of fragmens within a particular 
# size range flanked on either side by two different restiction sites
#
# for ddRAD estimates
#
# Author:  Tomas Hrbek
# Email: hrbek@evoamazon.net
# Date: 23.11.2014
# Version: 1.1
#######################################

__author__ = 'legal'

import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio.Restriction import *
from Bio.Alphabet.IUPAC import *

re_key = ['AanI', 'AarI', 'AasI', 'AatII', 'AbaSI', 'AbsI', 'Acc16I', 'Acc36I', 'Acc65I', 'AccB1I', 'AccB7I', 'AccBSI', 'AccI', 'AccII', 'AccIII', 'AceIII', 'AciI', 'AclI', 'AclWI', 'AcoI', 'AcsI', 'AcuI', 'AcvI', 'AcyI', 'AdeI', 'AfaI', 'AfeI', 'AfiI', 'AflII', 'AflIII', 'AgeI', 'AgsI', 'AhaIII', 'AhdI', 'AhlI', 'AjiI', 'AjnI', 'AjuI', 'AleI', 'AlfI', 'AllEnzymes', 'AloI', 'AluBI', 'AluI', 'Alw21I', 'Alw26I', 'Alw44I', 'AlwFI', 'AlwI', 'AlwNI', 'Ama87I', 'Analysis', 'Aor13HI', 'Aor51HI', 'AoxI', 'ApaBI', 'ApaI', 'ApaLI', 'ApeKI', 'ApoI', 'ApyPI', 'AquII', 'AquIII', 'AquIV', 'ArsI', 'AscI', 'AseI', 'Asi256I', 'AsiGI', 'AsiSI', 'Asp700I', 'Asp718I', 'AspA2I', 'AspBHI', 'AspLEI', 'AspS9I', 'AssI', 'AsuC2I', 'AsuHPI', 'AsuI', 'AsuII', 'AsuNHI', 'AvaI', 'AvaII', 'AvaIII', 'AvrII', 'AxyI', 'BaeGI', 'BaeI', 'BalI', 'BamHI', 'BanI', 'BanII', 'BarI', 'BasI', 'BauI', 'Bbr7I', 'BbrPI', 'BbsI', 'Bbv12I', 'BbvCI', 'BbvI', 'BbvII', 'BccI', 'Bce83I', 'BceAI', 'BcefI', 'BcgI', 'BciT130I', 'BciVI', 'BclI', 'BcnI', 'BcoDI', 'BcuI', 'BdaI', 'BetI', 'BfaI', 'BfiI', 'BfmI', 'BfoI', 'BfrI', 'BfuAI', 'BfuCI', 'BfuI', 'BglI', 'BglII', 'BinI', 'BisI', 'BlnI', 'BlpI', 'BlsI', 'BmcAI', 'Bme1390I', 'Bme18I', 'BmeDI', 'BmeRI', 'BmeT110I', 'BmgBI', 'BmgI', 'BmgT120I', 'BmiI', 'BmrFI', 'BmrI', 'BmsI', 'BmtI', 'BmuI', 'BoxI', 'BpiI', 'BplI', 'BpmI', 'Bpu10I', 'Bpu1102I', 'Bpu14I', 'BpuEI', 'BpuMI', 'BpvUI', 'Bsa29I', 'BsaAI', 'BsaBI', 'BsaHI', 'BsaI', 'BsaJI', 'BsaWI', 'BsaXI', 'BsbI', 'Bsc4I', 'BscAI', 'BscGI', 'Bse118I', 'Bse1I', 'Bse21I', 'Bse3DI', 'Bse8I', 'BseAI', 'BseBI', 'BseCI', 'BseDI', 'BseGI', 'BseJI', 'BseLI', 'BseMI', 'BseMII', 'BseNI', 'BsePI', 'BseRI', 'BseSI', 'BseX3I', 'BseXI', 'BseYI', 'BsgI', 'Bsh1236I', 'Bsh1285I', 'BshFI', 'BshNI', 'BshTI', 'BshVI', 'BsiEI', 'BsiHKAI', 'BsiHKCI', 'BsiI', 'BsiSI', 'BsiWI', 'BsiYI', 'BslFI', 'BslI', 'BsmAI', 'BsmBI', 'BsmFI', 'BsmI', 'BsnI', 'Bso31I', 'BsoBI', 'Bsp119I', 'Bsp120I', 'Bsp1286I', 'Bsp13I', 'Bsp1407I', 'Bsp143I', 'Bsp1720I', 'Bsp19I', 'Bsp24I', 'Bsp68I', 'BspACI', 'BspCNI', 'BspD6I', 'BspDI', 'BspEI', 'BspFNI', 'BspGI', 'BspHI', 'BspLI', 'BspLU11I', 'BspMI', 'BspMII', 'BspNCI', 'BspOI', 'BspPI', 'BspQI', 'BspT104I', 'BspT107I', 'BspTI', 'BsrBI', 'BsrDI', 'BsrFI', 'BsrGI', 'BsrI', 'BsrSI', 'BssAI', 'BssECI', 'BssHII', 'BssKI', 'BssMI', 'BssNAI', 'BssNI', 'BssSI', 'BssT1I', 'Bst1107I', 'Bst2BI', 'Bst2UI', 'Bst4CI', 'Bst6I', 'BstACI', 'BstAFI', 'BstAPI', 'BstAUI', 'BstBAI', 'BstBI', 'BstC8I', 'BstDEI', 'BstDSI', 'BstEII', 'BstENI', 'BstF5I', 'BstFNI', 'BstH2I', 'BstHHI', 'BstKTI', 'BstMAI', 'BstMBI', 'BstMCI', 'BstMWI', 'BstNI', 'BstNSI', 'BstOI', 'BstPAI', 'BstPI', 'BstSCI', 'BstSFI', 'BstSLI', 'BstSNI', 'BstUI', 'BstV1I', 'BstV2I', 'BstX2I', 'BstXI', 'BstYI', 'BstZ17I', 'BstZI', 'Bsu15I', 'Bsu36I', 'BsuI', 'BsuRI', 'BtgI', 'BtgZI', 'BthCI', 'BtrI', 'BtsCI', 'BtsI', 'BtsIMutI', 'BtuMI', 'BveI', 'Cac8I', 'CaiI', 'CauII', 'CchII', 'CchIII', 'CciI', 'CciNI', 'Cdi630V', 'CdiI', 'CdpI', 'CfoI', 'Cfr10I', 'Cfr13I', 'Cfr42I', 'Cfr9I', 'CfrI', 'Cgl13032I', 'Cgl13032II', 'ChaI', 'CjeFIII', 'CjeFV', 'CjeI', 'CjeNII', 'CjeNIII', 'CjeP659IV', 'CjePI', 'CjuI', 'CjuII', 'ClaI', 'CommOnly', 'CpoI', 'CseI', 'CsiI', 'Csp6I', 'CspAI', 'CspCI', 'CspI', 'CstMI', 'CviAII', 'CviJI', 'CviKI_1', 'CviQI', 'CviRI', 'DdeI', 'DinI', 'DpnI', 'DpnII', 'DraI', 'DraII', 'DraIII', 'DraRI', 'DrdI', 'DrdII', 'DriI', 'DsaI', 'DseDI', 'EaeI', 'EagI', 'Eam1104I', 'Eam1105I', 'EarI', 'EciI', 'Ecl136II', 'EclXI', 'Eco105I', 'Eco130I', 'Eco147I', 'Eco24I', 'Eco31I', 'Eco32I', 'Eco47I', 'Eco47III', 'Eco52I', 'Eco53kI', 'Eco57I', 'Eco57MI', 'Eco72I', 'Eco81I', 'Eco88I', 'Eco91I', 'EcoHI', 'EcoICRI', 'EcoNI', 'EcoO109I', 'EcoO65I', 'EcoRI', 'EcoRII', 'EcoRV', 'EcoT14I', 'EcoT22I', 'EcoT38I', 'EgeI', 'EheI', 'ErhI', 'EsaBC3I', 'EsaSSI', 'Esp3I', 'EspI', 'FaeI', 'FaiI', 'FalI', 'FaqI', 'FatI', 'FauI', 'FauNDI', 'FbaI', 'FblI', 'FinI', 'FmuI', 'Fnu4HI', 'FnuDII', 'FokI', 'FormattedSeq', 'FriOI', 'FseI', 'Fsp4HI', 'FspAI', 'FspBI', 'FspEI', 'FspI', 'GauT27I', 'GdiII', 'GlaI', 'GluI', 'GsaI', 'GsuI', 'HaeI', 'HaeII', 'HaeIII', 'HapII', 'HauII', 'HgaI', 'HgiAI', 'HgiCI', 'HgiEII', 'HgiJII', 'HhaI', 'Hin1I', 'Hin1II', 'Hin4I', 'Hin4II', 'Hin6I', 'HinP1I', 'HincII', 'HindII', 'HindIII', 'HinfI', 'HpaI', 'HpaII', 'HphI', 'Hpy166II', 'Hpy178III', 'Hpy188I', 'Hpy188III', 'Hpy8I', 'Hpy99I', 'Hpy99XIII', 'Hpy99XIV', 'HpyAV', 'HpyCH4III', 'HpyCH4IV', 'HpyCH4V', 'HpyF10VI', 'HpyF3I', 'HpySE526I', 'Hsp92I', 'Hsp92II', 'HspAI', 'Jma19592I', 'KasI', 'KflI', 'Kpn2I', 'KpnI', 'KroI', 'Ksp22I', 'Ksp632I', 'KspAI', 'KspI', 'Kzo9I', 'LguI', 'LpnI', 'LpnPI', 'Lsp1109I', 'LweI', 'MabI', 'MaeI', 'MaeII', 'MaeIII', 'MalI', 'MaqI', 'MauBI', 'MbiI', 'MboI', 'MboII', 'McaTI', 'McrI', 'MfeI', 'MflI', 'MhlI', 'MjaIV', 'MkaDII', 'MlsI', 'MluCI', 'MluI', 'MluNI', 'Mly113I', 'MlyI', 'MmeI', 'MnlI', 'Mph1103I', 'MreI', 'MroI', 'MroNI', 'MroXI', 'MscI', 'MseI', 'MslI', 'Msp20I', 'MspA1I', 'MspCI', 'MspI', 'MspJI', 'MspR9I', 'MssI', 'MstI', 'MunI', 'Mva1269I', 'MvaI', 'MvnI', 'MvrI', 'MwoI', 'NaeI', 'NarI', 'NciI', 'NcoI', 'NdeI', 'NdeII', 'NgoAVIII', 'NgoMIV', 'NhaXI', 'NheI', 'NlaCI', 'NlaIII', 'NlaIV', 'Nli3877I', 'NmeAIII', 'NmeDI', 'NmuCI', 'NonComm', 'NotI', 'NruI', 'NsbI', 'NsiI', 'NspBII', 'NspI', 'NspV', 'OliI', 'PabI', 'PacI', 'PaeI', 'PaeR7I', 'PagI', 'PalAI', 'PasI', 'PauI', 'PceI', 'PciI', 'PciSI', 'PcsI', 'PctI', 'PdiI', 'PdmI', 'PenI', 'PfeI', 'Pfl1108I', 'Pfl23II', 'PflFI', 'PflMI', 'PfoI', 'PinAI', 'PlaDI', 'Ple19I', 'PleI', 'PluTI', 'PmaCI', 'PmeI', 'PmlI', 'PpiI', 'PpsI', 'Ppu10I', 'Ppu21I', 'PpuMI', 'PrintFormat', 'PscI', 'PshAI', 'PshBI', 'PsiI', 'Psp03I', 'Psp124BI', 'Psp1406I', 'Psp5II', 'Psp6I', 'PspCI', 'PspEI', 'PspGI', 'PspLI', 'PspN4I', 'PspOMI', 'PspOMII', 'PspPI', 'PspPPI', 'PspPRI', 'PspXI', 'PsrI', 'PssI', 'PstI', 'PstNI', 'PsuI', 'PsyI', 'PteI', 'PvuI', 'PvuII', 'R2_BceSIV', 'RanaConfig', 'RceI', 'RdeGBI', 'RdeGBII', 'RdeGBIII', 'Restriction', 'RestrictionBatch', 'Restriction_Dictionary', 'RflFIII', 'RgaI', 'RigI', 'RlaI', 'RleAI', 'RpaB5I', 'RpaBI', 'RpaI', 'RpaTI', 'RruI', 'RsaI', 'RsaNI', 'RseI', 'Rsr2I', 'RsrII', 'SacI', 'SacII', 'SalI', 'SanDI', 'SapI', 'SaqAI', 'SatI', 'Sau3AI', 'Sau96I', 'SauI', 'SbfI', 'ScaI', 'SchI', 'SciI', 'ScrFI', 'SdaI', 'SdeAI', 'SdeOSI', 'SduI', 'SecI', 'SelI', 'SetI', 'SexAI', 'SfaAI', 'SfaNI', 'SfcI', 'SfeI', 'SfiI', 'SfoI', 'Sfr274I', 'Sfr303I', 'SfuI', 'SgeI', 'SgfI', 'SgrAI', 'SgrBI', 'SgrDI', 'SgrTI', 'SgsI', 'SimI', 'SlaI', 'SmaI', 'SmiI', 'SmiMI', 'SmlI', 'SmoI', 'SnaBI', 'SnaI', 'Sno506I', 'SpeI', 'SphI', 'SplI', 'SpoDI', 'SrfI', 'Sse232I', 'Sse8387I', 'Sse8647I', 'Sse9I', 'SseBI', 'SsiI', 'SspD5I', 'SspDI', 'SspI', 'SstE37I', 'SstI', 'Sth132I', 'Sth302II', 'StrI', 'StsI', 'StuI', 'StyD4I', 'StyI', 'SwaI', 'TaaI', 'TaiI', 'TaqI', 'TaqII', 'TasI', 'TatI', 'TauI', 'TfiI', 'Tru1I', 'Tru9I', 'TscAI', 'TseFI', 'TseI', 'TsoI', 'Tsp45I', 'Tsp4CI', 'TspDTI', 'TspEI', 'TspGWI', 'TspMI', 'TspRI', 'TssI', 'TstI', 'TsuI', 'Tth111I', 'Tth111II', 'UbaF11I', 'UbaF12I', 'UbaF13I', 'UbaF14I', 'UbaF9I', 'UbaPI', 'UcoMSI', 'UnbI', 'Van91I', 'Vha464I', 'VneI', 'VpaK11AI', 'VpaK11BI', 'VspI', 'WviI', 'XagI', 'XapI', 'XbaI', 'XceI', 'XcmI', 'XhoI', 'XhoII', 'XmaI', 'XmaIII', 'XmaJI', 'XmiI', 'XmnI', 'XspI', 'YkrI', 'ZraI', 'ZrmI', 'Zsp2I']
re_value = [AanI, AarI, AasI, AatII, AbaSI, AbsI, Acc16I, Acc36I, Acc65I, AccB1I, AccB7I, AccBSI, AccI, AccII, AccIII, AceIII, AciI, AclI, AclWI, AcoI, AcsI, AcuI, AcvI, AcyI, AdeI, AfaI, AfeI, AfiI, AflII, AflIII, AgeI, AgsI, AhaIII, AhdI, AhlI, AjiI, AjnI, AjuI, AleI, AlfI, AllEnzymes, AloI, AluBI, AluI, Alw21I, Alw26I, Alw44I, AlwFI, AlwI, AlwNI, Ama87I, Analysis, Aor13HI, Aor51HI, AoxI, ApaBI, ApaI, ApaLI, ApeKI, ApoI, ApyPI, AquII, AquIII, AquIV, ArsI, AscI, AseI, Asi256I, AsiGI, AsiSI, Asp700I, Asp718I, AspA2I, AspBHI, AspLEI, AspS9I, AssI, AsuC2I, AsuHPI, AsuI, AsuII, AsuNHI, AvaI, AvaII, AvaIII, AvrII, AxyI, BaeGI, BaeI, BalI, BamHI, BanI, BanII, BarI, BasI, BauI, Bbr7I, BbrPI, BbsI, Bbv12I, BbvCI, BbvI, BbvII, BccI, Bce83I, BceAI, BcefI, BcgI, BciT130I, BciVI, BclI, BcnI, BcoDI, BcuI, BdaI, BetI, BfaI, BfiI, BfmI, BfoI, BfrI, BfuAI, BfuCI, BfuI, BglI, BglII, BinI, BisI, BlnI, BlpI, BlsI, BmcAI, Bme1390I, Bme18I, BmeDI, BmeRI, BmeT110I, BmgBI, BmgI, BmgT120I, BmiI, BmrFI, BmrI, BmsI, BmtI, BmuI, BoxI, BpiI, BplI, BpmI, Bpu10I, Bpu1102I, Bpu14I, BpuEI, BpuMI, BpvUI, Bsa29I, BsaAI, BsaBI, BsaHI, BsaI, BsaJI, BsaWI, BsaXI, BsbI, Bsc4I, BscAI, BscGI, Bse118I, Bse1I, Bse21I, Bse3DI, Bse8I, BseAI, BseBI, BseCI, BseDI, BseGI, BseJI, BseLI, BseMI, BseMII, BseNI, BsePI, BseRI, BseSI, BseX3I, BseXI, BseYI, BsgI, Bsh1236I, Bsh1285I, BshFI, BshNI, BshTI, BshVI, BsiEI, BsiHKAI, BsiHKCI, BsiI, BsiSI, BsiWI, BsiYI, BslFI, BslI, BsmAI, BsmBI, BsmFI, BsmI, BsnI, Bso31I, BsoBI, Bsp119I, Bsp120I, Bsp1286I, Bsp13I, Bsp1407I, Bsp143I, Bsp1720I, Bsp19I, Bsp24I, Bsp68I, BspACI, BspCNI, BspD6I, BspDI, BspEI, BspFNI, BspGI, BspHI, BspLI, BspLU11I, BspMI, BspMII, BspNCI, BspOI, BspPI, BspQI, BspT104I, BspT107I, BspTI, BsrBI, BsrDI, BsrFI, BsrGI, BsrI, BsrSI, BssAI, BssECI, BssHII, BssKI, BssMI, BssNAI, BssNI, BssSI, BssT1I, Bst1107I, Bst2BI, Bst2UI, Bst4CI, Bst6I, BstACI, BstAFI, BstAPI, BstAUI, BstBAI, BstBI, BstC8I, BstDEI, BstDSI, BstEII, BstENI, BstF5I, BstFNI, BstH2I, BstHHI, BstKTI, BstMAI, BstMBI, BstMCI, BstMWI, BstNI, BstNSI, BstOI, BstPAI, BstPI, BstSCI, BstSFI, BstSLI, BstSNI, BstUI, BstV1I, BstV2I, BstX2I, BstXI, BstYI, BstZ17I, BstZI, Bsu15I, Bsu36I, BsuI, BsuRI, BtgI, BtgZI, BthCI, BtrI, BtsCI, BtsI, BtsIMutI, BtuMI, BveI, Cac8I, CaiI, CauII, CchII, CchIII, CciI, CciNI, Cdi630V, CdiI, CdpI, CfoI, Cfr10I, Cfr13I, Cfr42I, Cfr9I, CfrI, Cgl13032I, Cgl13032II, ChaI, CjeFIII, CjeFV, CjeI, CjeNII, CjeNIII, CjeP659IV, CjePI, CjuI, CjuII, ClaI, CommOnly, CpoI, CseI, CsiI, Csp6I, CspAI, CspCI, CspI, CstMI, CviAII, CviJI, CviKI_1, CviQI, CviRI, DdeI, DinI, DpnI, DpnII, DraI, DraII, DraIII, DraRI, DrdI, DrdII, DriI, DsaI, DseDI, EaeI, EagI, Eam1104I, Eam1105I, EarI, EciI, Ecl136II, EclXI, Eco105I, Eco130I, Eco147I, Eco24I, Eco31I, Eco32I, Eco47I, Eco47III, Eco52I, Eco53kI, Eco57I, Eco57MI, Eco72I, Eco81I, Eco88I, Eco91I, EcoHI, EcoICRI, EcoNI, EcoO109I, EcoO65I, EcoRI, EcoRII, EcoRV, EcoT14I, EcoT22I, EcoT38I, EgeI, EheI, ErhI, EsaBC3I, EsaSSI, Esp3I, EspI, FaeI, FaiI, FalI, FaqI, FatI, FauI, FauNDI, FbaI, FblI, FinI, FmuI, Fnu4HI, FnuDII, FokI, FormattedSeq, FriOI, FseI, Fsp4HI, FspAI, FspBI, FspEI, FspI, GauT27I, GdiII, GlaI, GluI, GsaI, GsuI, HaeI, HaeII, HaeIII, HapII, HauII, HgaI, HgiAI, HgiCI, HgiEII, HgiJII, HhaI, Hin1I, Hin1II, Hin4I, Hin4II, Hin6I, HinP1I, HincII, HindII, HindIII, HinfI, HpaI, HpaII, HphI, Hpy166II, Hpy178III, Hpy188I, Hpy188III, Hpy8I, Hpy99I, Hpy99XIII, Hpy99XIV, HpyAV, HpyCH4III, HpyCH4IV, HpyCH4V, HpyF10VI, HpyF3I, HpySE526I, Hsp92I, Hsp92II, HspAI, Jma19592I, KasI, KflI, Kpn2I, KpnI, KroI, Ksp22I, Ksp632I, KspAI, KspI, Kzo9I, LguI, LpnI, LpnPI, Lsp1109I, LweI, MabI, MaeI, MaeII, MaeIII, MalI, MaqI, MauBI, MbiI, MboI, MboII, McaTI, McrI, MfeI, MflI, MhlI, MjaIV, MkaDII, MlsI, MluCI, MluI, MluNI, Mly113I, MlyI, MmeI, MnlI, Mph1103I, MreI, MroI, MroNI, MroXI, MscI, MseI, MslI, Msp20I, MspA1I, MspCI, MspI, MspJI, MspR9I, MssI, MstI, MunI, Mva1269I, MvaI, MvnI, MvrI, MwoI, NaeI, NarI, NciI, NcoI, NdeI, NdeII, NgoAVIII, NgoMIV, NhaXI, NheI, NlaCI, NlaIII, NlaIV, Nli3877I, NmeAIII, NmeDI, NmuCI, NonComm, NotI, NruI, NsbI, NsiI, NspBII, NspI, NspV, OliI, PabI, PacI, PaeI, PaeR7I, PagI, PalAI, PasI, PauI, PceI, PciI, PciSI, PcsI, PctI, PdiI, PdmI, PenI, PfeI, Pfl1108I, Pfl23II, PflFI, PflMI, PfoI, PinAI, PlaDI, Ple19I, PleI, PluTI, PmaCI, PmeI, PmlI, PpiI, PpsI, Ppu10I, Ppu21I, PpuMI, PrintFormat, PscI, PshAI, PshBI, PsiI, Psp03I, Psp124BI, Psp1406I, Psp5II, Psp6I, PspCI, PspEI, PspGI, PspLI, PspN4I, PspOMI, PspOMII, PspPI, PspPPI, PspPRI, PspXI, PsrI, PssI, PstI, PstNI, PsuI, PsyI, PteI, PvuI, PvuII, R2_BceSIV, RanaConfig, RceI, RdeGBI, RdeGBII, RdeGBIII, Restriction, RestrictionBatch, Restriction_Dictionary, RflFIII, RgaI, RigI, RlaI, RleAI, RpaB5I, RpaBI, RpaI, RpaTI, RruI, RsaI, RsaNI, RseI, Rsr2I, RsrII, SacI, SacII, SalI, SanDI, SapI, SaqAI, SatI, Sau3AI, Sau96I, SauI, SbfI, ScaI, SchI, SciI, ScrFI, SdaI, SdeAI, SdeOSI, SduI, SecI, SelI, SetI, SexAI, SfaAI, SfaNI, SfcI, SfeI, SfiI, SfoI, Sfr274I, Sfr303I, SfuI, SgeI, SgfI, SgrAI, SgrBI, SgrDI, SgrTI, SgsI, SimI, SlaI, SmaI, SmiI, SmiMI, SmlI, SmoI, SnaBI, SnaI, Sno506I, SpeI, SphI, SplI, SpoDI, SrfI, Sse232I, Sse8387I, Sse8647I, Sse9I, SseBI, SsiI, SspD5I, SspDI, SspI, SstE37I, SstI, Sth132I, Sth302II, StrI, StsI, StuI, StyD4I, StyI, SwaI, TaaI, TaiI, TaqI, TaqII, TasI, TatI, TauI, TfiI, Tru1I, Tru9I, TscAI, TseFI, TseI, TsoI, Tsp45I, Tsp4CI, TspDTI, TspEI, TspGWI, TspMI, TspRI, TssI, TstI, TsuI, Tth111I, Tth111II, UbaF11I, UbaF12I, UbaF13I, UbaF14I, UbaF9I, UbaPI, UcoMSI, UnbI, Van91I, Vha464I, VneI, VpaK11AI, VpaK11BI, VspI, WviI, XagI, XapI, XbaI, XceI, XcmI, XhoI, XhoII, XmaI, XmaIII, XmaJI, XmiI, XmnI, XspI, YkrI, ZraI, ZrmI, Zsp2I]
res = dict(zip(re_key, re_value))

parser = argparse.ArgumentParser(description='Script to estimate number of RADs and ddRADs given a genome and size range')
parser.add_argument('input', help='input genome file name')
parser.add_argument('min', help='minimum fragment length, default 350 bp', type=int, nargs='?', default=350)
parser.add_argument('max', help='maximum fragment length, default 450 bp', type=int, nargs='?', default=450)
parser.add_argument('extract', help='extract fragments from genome, default TRUE', type=int, nargs='?', default=1)
parser.add_argument('-re1','--restriction_enzyme_1', help='rare restriction enzyme', required=False, default='SbfI')
parser.add_argument('-re2','--restriction_enzyme_2', help='common restriction enzyme', required=False, default='Csp6I')
parser.add_argument('-v', '--verbose', help='increase output verbosity', action='store_true')
args = parser.parse_args()

if args.verbose:
    print("verbosity turned on")


count_cuts = 0 #Initialize rad counter
count_contigs = 0  #Initialize contig counter
no_cut = 0 #Initialize no cut counter
range_seqs = [] #Setup an empty list of sequences
count_seq = 0 #Initialize counter
len_seqs = [] #Setup an empty list of sequence lengths
gc_seqs = [] #Setup an empty list of sequence GC%

input_handle = open(args.input, "rU")
min = args.min
max = args.max
re = args.restriction_enzyme_1
re1 = res.get(re)
re = args.restriction_enzyme_2
re2 = res.get(re)

if args.restriction_enzyme_1 not in re_key:
	print ("The restriction enzyme {0} does not exist.".format(args.restriction_enzyme_1))
if args.restriction_enzyme_2 not in re_key:
	print ("The restriction enzyme {0} does not exist.".format(args.restriction_enzyme_2))
	
if args.restriction_enzyme_1 not in re_key or args.restriction_enzyme_2 not in re_key:
	print ("Execute program again with valid enzymes.")
	input_handle.close()
	sys.exit(1)


print ("Analysing {0} genome with {1} and {2} enzymes".format(args.input, re1, re2))
print ("Please wait till finished".format())


for record in SeqIO.parse(input_handle, "fasta") :
	count_contigs += 1
	cut1 = re1.catalyse(record.seq) #replace first enzyme here
	if len(cut1) > 1 :
		count_cuts += len(cut1)
		for cut in cut1 :
			cut2 = re2.catalyse(cut) #replace second enzyme here
			if len(cut2) > 1 :
				if len(cut2[0]) < max and len(cut2[0]) > min : #pass min and max size as argument
					count_seq += 1
					id = '>RAD' + str(count_seq)
					range_seqs.append(id)
					range_seqs.append(cut2[0].reverse_complement())
					#range_seqs.append(SeqRecord(Seq(cut2[0], IUPACAmbiguousDNA()), id = args.input))
					len_seqs.append(len(cut2[0]))
				if len(cut2[-1]) < max and len(cut2[-1]) > min :#pass min and max size as argument
					count_seq += 1
					id = '>RAD' + str(count_seq)
					range_seqs.append(id)
					range_seqs.append(cut2[-1])
					#range_seqs.append(SeqRecord(Seq(cut2[-1].reverse_complement(), IUPACAmbiguousDNA()), id = args.input))
					len_seqs.append(len(cut2[-1]))
	else :
		no_cut += 1

input_handle.close()

if args.extract:
	print ("********".format())
	print ("Extracting ddRAD fragments to a file".format())
	output_handle = open("extracted_seq.fas", "w")
	output_handle.write("\n".join(str(x) for x in range_seqs))
	output_handle.close()
	output_handle = open("extracted_length", "w")
	output_handle.write("\n".join(str(x) for x in len_seqs))
	output_handle.close()
	

print ("********".format())
print ("Genome analyzed: {0}".format(args.input))
print ("Minimum fragment size: {0}".format(min))
print ("Maximum fragment size: {0}".format(max))
print ("Rare restriction enzyme: {0}".format(re1))
print ("Common restriction enzyme: {0}".format(re2))
print ("Number of contigs: {0}".format(count_contigs))
print ("Number of contigs without {0} cuts: {1}".format(re1, no_cut))
print ("Number of RAD {0} cuts: {1}".format(re1, count_cuts))
print ("Number of RADseq fragments: {0}".format((count_cuts-(count_contigs-no_cut))*2))
print ("Number of ddRAD fragments within range: {0}".format(count_seq))
print ("********".format())

    
#output_handle = open("range_sequences.fas", "w")
#SeqIO.write(range_seqs, output_handle, "fasta")
#output_handle.close()
