#!/usr/bin/env python
# -*- coding: utf-8	-*-

import os, sys , argparse, glob
import subprocess
from utils import listFiles, run_shell_command, associationMRISession, bcolors
import bids.layout
import scipy.stats
import numpy as	np
import nibabel as nb
import json
import gzip
import shutil


### NEED TO BE DONE: Frequency and StartTime

def	get_arguments():
	parser	= argparse.ArgumentParser(
			formatter_class=argparse.RawDescriptionHelpFormatter,
			description="",
			epilog="""
			Convert Physiologic	data to	BIDS format	(PNM-FSL compatible)

			Input: .cfg

			Although no arguments is	mandatory there	is an order, if	case 1 then	it won't do	case 2 or 3	and	so on
			1. dataset Folder
			2. subject Folder
			3. single File
			""")

	parser.add_argument(
			"-t", "--tar",
			required=True, nargs="+",
			help="Tar files (output from Siemens MRI) single archive files can be used as input: physio2bids.py -t file.tar.gz -p S02 -o bidsFolder",
			)

	parser.add_argument(
			"-s", "--subject",
			required=True, nargs="+",
			help="Subject name (i.e BIDS format)",
			)

	parser.add_argument(
			"-o", "--bidsfolder",
			required=True, nargs="+",
			help="BIDS folder",
			)

	parser.add_argument(
			"-v", "--verbose",
			required=False, nargs="+",
			help="Verbose",
			)

	args =	parser.parse_args()
	if	len(sys.argv) == 1:
		parser.print_help()
		sys.exit()
	else:
		return args

class physio2bids(object):
	"""
	"""
	def __init__(
			self, tar, subject,
		bidsfolder,	verbose=False, log_level="INFO"):
		self.tar = tar
		self.subject = subject[0]
		self.bidsFolder = bidsfolder[0]
		self.verbose = verbose

		self.iCardiac	= ''
		self.iRespiratory	= ''
		self.iTrigger	= ''
		self.iNIFTI = ''

		# Additional information
		self.iNIFTIdim = ''
		self.frequency = ''
		self.startTime = '\"None\"'
		self.TR = ''

		#	For	each functional	MRI	we'll output a json	and	TSV	file with "cardiac", "respiratory",	"trigger" recorded
		self.oJSON = ''
		self.oTSV = ''

	def getPhysiologicalFiles(self):

		physioFiles = []

		if self.tar:
			for nTar in self.tar:
				cmd = 'tar -xf ' + nTar	+ '	-C '+ os.path.dirname(nTar) + ' -v'
				output = run_shell_command(cmd)

				tmpPhysio = []

				'''
				Get all base name of each file
				'''
				for x in output.split():  # Starting at 1 to remove folder name
					tmpPhysio.append(os.path.dirname(nTar) + os.path.sep + os.path.splitext(x.decode("utf-8"))[0])

				'''
				Remove duplicates
				'''
				tmpPhysio = list(set(tmpPhysio))

				'''
				Get all files if .ext .puls and .resp exist
				'''
				for nPhysio in tmpPhysio:
					if os.path.exists(nPhysio+'.ext') and os.path.exists(nPhysio+'.puls') and os.path.exists(nPhysio+'.resp'):
						physioFiles.append(nPhysio)

		elif self.sourcedata:
			pass

		elif self.subjectFolder:
			pass

		physioFiles.sort()

		if self.verbose:
			print(physioFiles)

		print(physioFiles)
		return physioFiles

	'''
	1. Associate each physiological data to a MRI session
	2. Get JSON file name -> self.oJSON
	3. Get TSV file name -> self.oTSV
	'''
	def associateWithMRSession(self):

		try:
			self.validateBIDSFolder()

		except OSError as exception:
			print('Not able to associate MRI acquisition')
			logger.error("Exception: {}".format(exception))
			logger.info("subprocess failed")

	def writeJSON(self):
		text = '{\n\"SamplingFrequency\": '+str(self.frequency)+',\n\"StartTime\": '+self.startTime+',\n\"Columns\": [\"cardiac\", \"respiratory\", \"trigger\"]\n}'
		oJSON = open(self.oJSON,"w")
		oJSON.write(text)
		oJSON.close()

	def run(self):
		'''
		Create a list	of triplet with	cardiac, respiratory and trigger input
		'''
		physioFiles = self.getPhysiologicalFiles()

		for nPhysioFile in physioFiles:

			print('#########')
			print(nPhysioFile)

			'''
			Initialization
			'''
			self.oJSON = ''
			self.oTSV = ''
			self.iNIFTI = ''

			self.baseName = nPhysioFile
			self.iCardiac = nPhysioFile + '.puls'
			self.iRespiratory = nPhysioFile + '.resp'
			self.iTrigger = nPhysioFile + '.ext'

			self.associateWithMRSession()

			if self.verbose:
				print('Output phyisio json: ' + self.oJSON)
				print('Output phyisio data: ' + self.oTSV)
				print('Associated nifti file: ' + self.iNIFTI)


			'''
			Extract TRIGGER, CARDIAC and RESPIRATORY DATA
			'''
			triggerIndexes, dataTRIGGER = self.createTRIGGER()

			dataCARDIAC = self.createCARDIAC(triggerIndexes)

			dataRESPIRATORY = self.createRespiratory(triggerIndexes)

			self.writeJSON()

			'''
			Combine everything
			'''
			self.combinePhysio(dataCARDIAC, dataRESPIRATORY, dataTRIGGER)
		#subprocess.call("/cerebro/cerebro1/dataset/bmpd/derivatives/2021-02_last_dicom_AP_acq/extract_first_last_dicom_AP_acqs.sh")


	def combinePhysio(self, dataCARDIAC, dataRESPIRATORY, dataTRIGGER):
		'''
		Concatenate and write covariates (Resp and then Puls)
		'''
		covariateData = np.vstack((dataCARDIAC, dataRESPIRATORY))
		covariateData = np.vstack((covariateData, dataTRIGGER))

		np.savetxt(self.oTSV, covariateData.transpose(), fmt="%g")

		with open(self.oTSV, 'rb') as f_in, gzip.open(self.oTSV + '.gz', 'wb') as f_out:
			shutil.copyfileobj(f_in, f_out)
		os.remove(self.oTSV)

		if self.verbose:
			print('Save combination all physio data: {}'.format(self.oTSV))

	def createCARDIAC(self, triggerIndexes):
		'''
		Creation of the BIDS version of the cardiac data (puls file)
		'''
		if self.verbose:
			print('Load CARDIAC File: {}'.format(self.iCardiac))

		cardiacFile  = open(self.iCardiac, 'r')

		for line in cardiacFile:
			values = line
			break

		# Convert string to int
		origCARDIAC = np.fromstring(values[values.rfind(':')+1:],sep=' ')

		# Reverse array and select data
		origReverseCARDIAC = origCARDIAC[::-1]
		selectedReverseCARDIAC = origReverseCARDIAC[triggerIndexes[0]:triggerIndexes[1]+1]

		# Reverse again
		selectedCARDIAC = selectedReverseCARDIAC[::-1]

		# Detection
		selectedCARDIAC[selectedCARDIAC!=5000] = 0
		selectedCARDIAC[selectedCARDIAC==5000] = 1

		# Save selected CARDIAC DATA
		with open(self.iCardiac.replace('.puls', '_corrected.puls'), "w") as text_file:
			selectedCARDIAC.tofile(text_file, sep="\n", format="%g")

		if self.verbose:
			print('Save CARDIAC file: {}'.format(self.iCardiac.replace('.puls', '_corrected.puls')))

		return selectedCARDIAC

	def createTRIGGER(self):
		'''
		Creation of the BIDS version of the trigger data (EXT file)
		'''
		print('Load TRIGGER file: {}'.format(self.iTrigger))
		triggerFile  = open(self.iTrigger, 'r')

		for line in triggerFile:
			values = line
			break

		# Convert string to int
		origTRIGGER = np.fromstring(values, sep=' ')

		# Get first and last occurance of 5000
		indexes5000 = np.where(origTRIGGER[::-1]==5000)[0]

		# Find Mode
		mode5000 = mode5000 = int(scipy.stats.mode(np.diff(indexes5000), keepdims=True).mode) #int(scipy.stats.mode(np.diff(indexes5000)).mode)

		# Boolean of values close to Mode
		boolMode5000 = np.abs(np.diff(indexes5000) - mode5000)<10

		firstIdx = -1
		groupIdx = []
		groupLength = []
		for idx, val in enumerate(boolMode5000):
			if val:
				if firstIdx<0:
					firstIdx = idx
			elif firstIdx>=0:
				lastIdx = idx + 1
				groupIdx.append([firstIdx, lastIdx])
				groupLength.append(lastIdx-firstIdx)
				firstIdx = -1

		if firstIdx>=0:
			groupIdx.append([firstIdx, idx + 2])
			groupLength.append(idx + 1-firstIdx)

		# Correct group
		print(groupIdx)
		print(groupLength)
		print(self.iNIFTIdim)
		groupValid = [i for i,x in enumerate(groupLength) if np.abs(x-self.iNIFTIdim[3])<10]

		if len(groupValid) != 1:
			print(groupIdx)
			print(groupLength)
			print(self.iNIFTIdim[3])
			raise Exception('Trigger file need to be cropped because it contains multiple acquisition')
		else:
			groupValid = groupValid[0]


		# Crop to correct indexes
		indexes5000 = indexes5000[groupIdx[groupValid][0]:groupIdx[groupValid][1]]

		if self.verbose:
			if np.abs(len(indexes5000)-self.iNIFTIdim[3])<5: # Correct
				print('{}Number of TR: {}{}'.format('\033[92m',str(len(indexes5000)),'\033[0m'))
				print('{}Number of fMRI volumes: {}{}'.format('\033[92m',str(self.iNIFTIdim[3]),'\033[0m'))
			else:
				print('{}Number of TR: {}{}'.format('\033[91m',str(len(indexes5000)),'\033[0m'))
				print('{}Number of fMRI volumes: {}{}'.format('\033[91m',str(self.iNIFTIdim[3]),'\033[0m'))

		self.frequency = int(mode5000/float(self.TR))

		if self.verbose:
			print('Frequency: {} Hz'.format(self.frequency))

		first5000ReverseIndex = indexes5000[0]
		last5000ReverseIndex = indexes5000[-1]

		# Reverse array and select data
		origReversTRIGGER = origTRIGGER[::-1]
		selectedReverseTRIGGER = origReversTRIGGER[first5000ReverseIndex:last5000ReverseIndex+1]

		# Reverse again
		selectedTRIGGER = selectedReverseTRIGGER[::-1]

		# Detection
		selectedTRIGGER[selectedTRIGGER!=5000] = 0
		selectedTRIGGER[selectedTRIGGER==5000] = 1

		# Save selected TRIGGER DATA
		with open(self.iTrigger.replace('.ext', '_corrected.ext'), "w") as text_file:
			selectedTRIGGER.tofile(text_file, sep="\n", format="%g")

		if self.verbose:
			print('Save TRIGGER file: {}'.format(self.iTrigger.replace('.ext', '_corrected.ext')))

		return [first5000ReverseIndex, last5000ReverseIndex] , selectedTRIGGER

	def createRespiratory(self, triggerIndexes):
		'''
		Creation of the BIDS version of the Respiratory data (RESP file)
		'''
		print('Load RESPIRATORY file: {}'.format(self.iRespiratory))
		respiratoryFile  = open(self.iRespiratory, 'r')

		for line in respiratoryFile:
			values = line
			break

		# Convert string to int
		origRespiratory = np.fromstring(values[values.rfind(': ')+1:],sep=' ')

		indexFirstValue = np.where(origRespiratory==6002)[0]
		origRespiratory = origRespiratory[indexFirstValue[0]:]

		# Reverse array and select data
		origReverseRespiratory = origRespiratory[::-1]
		selectedReverseRespiratory = origReverseRespiratory[triggerIndexes[0]:triggerIndexes[1]+1]

		# Reverse again
		selectedRespiratory = selectedReverseRespiratory[::-1]

		# Detection already done
		# Remove peaks and put mean(t-1,t+1)
		selectedRespiratoryPeaks = np.where(selectedRespiratory==5000)[0]

		for currentPeak in selectedRespiratoryPeaks:
			if int(currentPeak)==0:
				selectedRespiratory[currentPeak] = selectedRespiratory[currentPeak+1]
			elif int(currentPeak)==selectedRespiratory.shape[0]-1:
				selectedRespiratory[currentPeak] = selectedRespiratory[currentPeak-1]
			else:
				selectedRespiratory[currentPeak] = np.mean([selectedRespiratory[currentPeak-1], selectedRespiratory[currentPeak+1]])

		'''
		Remove mean
		'''
		selectedRespiratory = selectedRespiratory-np.mean(selectedRespiratory)

		# Save selected RESP DATA
		with open(self.iRespiratory.replace('.resp', '_corrected.resp'), "w") as text_file:
			selectedRespiratory.tofile(text_file, sep="\n", format="%g")

		if self.verbose:
			print('Save RESPIRATORY file: {}'.format(self.iRespiratory.replace('.resp', '_corrected.resp')))

		return selectedRespiratory

	def validateBIDSFolder(self):
		'''
		Load BIDS layout
		'''
		layout = bids.layout.BIDSLayout(self.bidsFolder)

		print(layout)

		try:
			'''
			Because so many files were not named correctly
			'''
			if len(os.path.basename(self.baseName).split('_')) == 1:
				taskName = os.path.basename(self.baseName)[0]
			else:
				taskName = "_".join(os.path.basename(self.baseName).split('_')[1:])

			print('Task: {}'.format(taskName))

			currKey = [i for i in associationMRISession.keys() if taskName in i][0]

			print('Keys: {} : {}'.format(currKey, ' , '.join(associationMRISession[currKey])))

			fMRI_func_nifti = layout.get(subject=self.subject, datatype='func', task=associationMRISession[currKey][0], extension='.nii.gz', return_type='file')

			print(fMRI_func_nifti)
			
			if len(fMRI_func_nifti) == 1:
				self.oJSON = fMRI_func_nifti[0].replace('bold.nii.gz', 'physio.json')
				self.oTSV = fMRI_func_nifti[0].replace('bold.nii.gz', 'physio.tsv')
				self.iNIFTI = fMRI_func_nifti[0]
			else:
				'''
				Find Session
				'''
				print(associationMRISession[currKey])
				fMRI_func_nifti_run = [i for i in fMRI_func_nifti if associationMRISession[currKey][1] in i]
				fMRI_func_nifti_run
				if len(fMRI_func_nifti_run) == 1:
					self.oJSON = fMRI_func_nifti_run[0].replace('bold.nii.gz', 'physio.json')
					self.oTSV = fMRI_func_nifti_run[0].replace('bold.nii.gz', 'physio.tsv')
					self.iNIFTI = fMRI_func_nifti_run[0]
				else:
					print('HOUSTON we are not able to find a correct MRI session')

			'''
			Get dimension of the nifti file
			'''
			nifti_data = nb.load(self.iNIFTI)
			self.iNIFTIdim = nifti_data.shape
			val = True

			'''
			Get TR
			'''
			with open(self.iNIFTI.replace('.nii.gz', '.json')) as jsonFile:
				jsonData = json.load(jsonFile)

			self.TR = jsonData['RepetitionTime']

			if self.verbose:
				print('TR: {}'.format(self.TR))

		except Exception as e:
			print(e)
			val = False

		return val

def	main():
	"""Let's go"""
	args =	get_arguments()
	app = physio2bids(**vars(args))
	return	app.run()

if __name__	== '__main__':
	sys.exit(main())
