import os,sys

class options:

	#constructor
	def __init__(self,fname):
		#option defaults
		self.cellsNumber = 64
		self.photoBinsNumber = 100
		self.network = ""
		self.energyMin = None
		self.energyMax = None
		self.plotRates = True
		self.useEntropyProduction = False
		self.species = []
		self.usePhotochemistry = True

		#required casting for options
		integerType = ["cellsNumber","photoBinsNumber"]
		floatType = ["energyMin","energyMax"]
		boolType = ["plotRates","useEntropyProduction","usePhotochemistry"]
		listType = ["species"]

		#check if file exists
		if(not(os.path.isfile(fname))):
			print ("ERROR: option file not found: "+fname)
			sys.exit()

		print ("reading option file "+fname)
		#open option file
		fh = open(fname,"r")
		for row in fh:
			srow = row.strip()
			if(srow==""): continue
			if(srow.startswith("#")): continue
			(option,value) = [x.strip() for x in srow.split("=") if x!=""]
			if(option in integerType): value = int(value)
			if(option in floatType): value = float(value)
			if(option in boolType): value = (value.upper()=="T")
			if(option in listType): value = [x.strip() for x in value.split(",")]
			#convert option string into variable
			setattr(self, option, value)
		fh.close()

