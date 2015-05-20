# Used libraries
import argparse
import os.path
import re
import errno
import datetime
import subprocess
import sys
from io import StringIO
from subprocess import call
from pyfasta import Fasta

# Debugging
DEBUGGING = True

# Define argument parser
parser = argparse.ArgumentParser()
parser.add_argument('-mic', '--mic', dest='MIC')
parser.add_argument('-re', '--re', dest='RE')
parser.add_argument('-mac', '--mac', dest='MAC')
parser.add_argument('-o', '--o', dest='OUT')

# Get arguments
if DEBUGGING:
	args = parser.parse_args('-re AAAACCCC -mic Test_Files/mic.fasta -mac Test_Files/mac.fasta -o ../Output'.split())
else:
	args = parser.parse_args()
	
regExp = args.RE
MICfile = args.MIC
MACfile = args.MAC
Output_dir = args.OUT

# Check if files are real files
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

LogFile = None
# Create output directory if needed
try:
	os.makedirs(Output_dir)
	LogFile = open(Output_dir + '/log.txt', 'w')
	LogFile.write(datetime.datetime.now().strftime("%I:%M%p %B %d %Y") + ' - Directory ' + Output_dir + ' created\n')
except OSError as exception:
	if exception.errno == errno.EEXIST:
		LogFile = open(Output_dir + '/log.txt', 'w')
		LogFile.write(datetime.datetime.now().strftime("%I:%M%p %B %d %Y") + ' - Directory ' + Output_dir + ' found\n')
	else:
		raise

#-----------------------------------------------------------------------------------------------------
# Define log comment function for the ease of use
def logComment(comment):
	global LogFile
	LogFile.write(datetime.datetime.now().strftime("%I:%M%p %B %d %Y") + ' - ' + comment + '\n')
#-----------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------
# Define create directory function
def safeCreateDirectory(dir):
	try:
		os.makedirs(dir)
		logComment('Directory ' + dir + ' created')
	except OSError as exception:
		if exception.errno != errno.EEXIST:
			raise
#-----------------------------------------------------------------------------------------------------
		
# Print parsed arguments
if DEBUGGING:
	print(regExp, MICfile, MACfile)

# Check if BLAST is installed
output = None
try:
	#call('blastn -version')
	output = subprocess.check_output(["blastn", "-version"])
except:
	print('BLAST is not installed on the computer. Please install BLAST to run the program')
	exit()
	
logComment("BLAST version: " + output.decode(sys.stdout.encoding).split('\n')[0])

# Read MIC fasta file (We might need to surround with try/except block for error checking!!!!!!!!!!!!!)
mic_fasta = Fasta(MICfile)

if DEBUGGING:
	print(list(mic_fasta.keys()))

#Read MAC fasta file
mac_fasta = Fasta(MACfile)

if DEBUGGING:
	print(list(mac_fasta.keys()))


# Close all files
LogFile.close()












