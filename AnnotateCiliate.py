# Used libraries
import argparse
import os.path
import re
from pyfasta import Fasta

# Debugging
DEBUGGING = True

# Define argument parser
parser = argparse.ArgumentParser()
parser.add_argument('-mic', '--mic', dest='MIC')
parser.add_argument('-re', '--re', dest='RE')
parser.add_argument('-mac', '--mac', dest='MAC')

# Get arguments
if DEBUGGING:
	args = parser.parse_args('-re AAAACCCC -mic Test_Files/mic.fasta -mac Test_Files/mac.fasta'.split())
else:
	args = parser.parse_args()
	
regExp = args.RE
MICfile = args.MIC
MACfile = args.MAC

while not os.path.isfile(MICfile):
	print(MICfile, ' is not a file, please input a valid MIC file:')
	MICfile = input()

while not os.path.isfile(MACfile):
	print(MACfile, ' is not a file, please input a valid MAC file:')
	MACfile = input()
	
# Test whether a regulat expression is a valid one
re_comp = 0
is_Valid_Re = False
while not is_Valid_Re:
	try:
		re_comp = re.compile(regExp)
		is_Valid_Re = True
	except re.error:
		print('The regular expression is invalid, please type a valid regular expression:')
		regExp = input()

# Print parsed arguments
print(regExp, MICfile, MACfile)

# Read MIC fasta file
mic_fasta = Fasta(MICfile)
print(list(mic_fasta.keys()))

#Read MAC fasta file
mac_fasta = Fasta(MACfile)
print(list(mac_fasta.keys()))
