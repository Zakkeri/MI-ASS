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
'RoughBlastDust' : True,
'RoughBlastUngapped' : True,
'RoughCoverageLimit' : 5,
'FineBlastTask' : 'blastn-short',
'FineBlastWordSize' : 12,
'FineBlastDust' : False,
'FineBlastUngapped' : True,
'ThreadCount' : 2*multiprocessing.cpu_count(),
'DatabaseUpdate' : True
}