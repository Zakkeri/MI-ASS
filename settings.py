# This is a file with all settings and options that we have for the AnnotateCilitate program

import multiprocessing

Options = {
'Telomeres' : True, 
'MDS' : True,
'Pointers' : True, 
'TelomereEndLimit' : 100,
'TelomereOffset' : 5,
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
'MIC_Annotation_MDS_Overlap_Threshold' : .65,
#'Tel_Reg_Exp' : "AAAACCCCAAAACCCC"
'Tel_Reg_Exp' : "(A{1,4}(C{4}A{4})*C{1,4})|(C{1,4}(A{4}C{4})*A{1,4})|(T{1,4}(G{4}T{4})*G{1,4})|(G{1,4}(T{4}G{4})*T{1,4})"
}