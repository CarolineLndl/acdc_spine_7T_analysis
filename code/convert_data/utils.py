import os , sys, glob

associationMRISession = {
		"task-rest_acq-shimBase+3mm_bold" : ["rest", "shimBase+3mm"],
		"task-rest_acq-shimSlice+3mm_bold" : ["rest", "shimSlice+3mm"],
		"task-motor_acq-shimBase+3mm_bold" : ["motor", "shimBase+3mm"],
		"task-motor_acq-shimBase+3mm_bold_run-01" : ["motor", "shimBase+3mm_run-01"],
		"task-motor_acq-shimBase+3mm_bold_run-02" : ["motor", "shimBase+3mm_run-02"],
		"task-motor_acq-shimSlice+3mm_bold" : ["motor", "shimSlice+3mm"],
		"task-motor_acq-shimSlice+3mm_bold_run-01" : ["motor", "shimSlice+3mm_run-01"],
		"task-motor_acq-shimSlice+3mm_bold_run-02" : ["motor", "shimSlice+3mm_run-02"],
		"task-rest_acq-shimBase+1mm+sms2_bold" : ["rest", "shimBase+1mm+sms2"],
		"task-rest_acq-shimSlice+1mm+sms2_bold" : ["rest", "shimSlice+1mm+sms2"],
		"task-motor_acq-shimSlice+1mm+sms2_bold" : ["motor", "shimSlice+1mm+sms2"],

}


def	listFiles(rootFolder, iFilter=None,	allFiles=[]):
	for path, subdirs,	files in os.walk(rootFolder):
		for name in files:
			try:
				if not iFilter:
					if not 'edited'	in name:
						allFiles.append(os.path.join(path,	name))
						
				elif iFilter in name and not 'edited'	in name:
						allFiles.append(os.path.join(path,	name))

			except:
				continue

	return allFiles

def	run_shell_command(commandLine):
	import	logging
	import	shlex
	import	subprocess

	logger	= logging.getLogger("physio2bids")

	cmd = shlex.split(commandLine)
	logger.info("subprocess: {}".format(commandLine))

	output	= None
	try:
		process =	subprocess.Popen(
				cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		output, _	= process.communicate()

		try:
			logger.info("\n"	+ output.decode("utf-8"))
		except:
			logger.info(output)

	except OSError as exception:
		logger.error("Exception: {}".format(exception))
		logger.info("subprocess failed")

	return	output

class bcolors(object):
    normal = '\033[0m'
    red = '\033[91m'
    green = '\033[92m'
    yellow = '\033[93m'
    blue = '\033[94m'
    magenta = '\033[95m'
    cyan = '\033[96m'
    bold = '\033[1m'
    underline = '\033[4m'
	
    @classmethod
    def colors(cls):
        return [v for k, v in cls.__dict__.items() if not k.startswith("_") and k != "colors"]
	