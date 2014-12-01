import sys

class ColorText:
	INFOBLUE = '\033[1;94m'
	INFOGREEN = '\033[1;92m'
	WARNING = '\033[93m'
	ERROR = '\033[7;91m'
	ENDC = '\033[0m'

	def disable(self):
		self.HEADER = ''
		self.OKBLUE = ''
		self.OKGREEN = ''
		self.WARNING = ''
		self.FAIL = ''
		self.ENDC = ''

	def info(self, str, ostream="stdout"):
		if ostream == "stderr":
			sys.stderr.write(self.INFOGREEN + str + self.ENDC)
		elif ostream == "stdout":
			sys.stdout.write(self.INFOGREEN + str + self.ENDC)

	def error(self, str):
		sys.stderr.write(self.ERROR + str + self.ENDC)

	def warning(self, str, ostream):
		if ostream == "stderr":
			sys.stderr.write(self.WARNING + str + self.ENDC)
		elif ostream == "stdout":
			sys.stdout.write(self.WARNING + str + self.ENDC)

	def __debug(self):
		print ColorText.FAIL + "Error: No active frommets remain. Continue?" + ColorText.ENDC
		print "\033[1;34mRED TEXT\033[0m"
		print "\033[4;34mRED TEXT\033[0m"
		print "\033[5;34mRED TEXT\033[0m"
		print "\033[7;34mRED TEXT\033[0m"
