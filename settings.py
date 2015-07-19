# This is a file with all settings and options that we have for the AnnotateCilitate program

import multiprocessing

Options = {
'Telomeres' : True, 
'MDS' : True,
'Pointers' : True, 
'TelomereEndLimit' : 100,
'TelomereLength' : 10,
'BlastMaskLowercase' : True,
'RoughBlastTask' : 'megablast',
'RoughBlastWordSize' : 28,
'RoughBlastDust' : False,
'RoughBlastUngapped' : True,
'RoughCoverageLimit' : 5,
'FineBlastTask' : 'blastn-short',
'FineBlastWordSize' : 12,
'FineBlastDust' : False,
'FineBlastUngapped' : True,
'ThreadCount' : 2*multiprocessing.cpu_count(),
'DatabaseUpdate' : True,
'MIC_Coverage_Threshold' : 10,
'TelomericErrorTolerance' : 5,
'MDS_id_start' : 1,
#'Min_hsp_to_mds_ratio' : 0.8,
#'Min_mds_num_percentage' : 0.5,
#'Tel_Reg_Exp_5' : "((AA){0,1}(CCCCAA)+(CCCC){0,1})",
#'Tel_Reg_Exp_3' : "((TT){0,1}(GGGGTT)+(GGGG){0,1})"
'Tel_Reg_Exp_5' : "(A{0,4}(C{4}A{4})+C{0,4})|(C{0,4}(A{4}C{4})+A{0,4})|(A{1,4}C{4}A{1,4})|(C{1,4}A{4}C{1,4})",
'Tel_Reg_Exp_3' : "(T{0,4}(G{4}T{4})+G{0,4})|(G{0,4}(T{4}G{4})+T{0,4})|(T{1,4}G{4}T{1,4})|(G{1,4}T{4}G{1,4})"
}
