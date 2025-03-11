import patmo_string
import patmo_options
import patmo_commons

def buildPhotoRates(network,options):

	allRates = ""
	for reaction in network.photoReactions:
		allRates +=  "!"+reaction.getVerbatim()+"\n"
		allRates +=  "krate(:,"+str(reaction.index)+") = "+reaction.rate+"\n\n"
	
	dE = 0.0
	dE += (options.wavelengMax-options.wavelengMin)/options.photoBinsNumber
	res =""
	res += "dE = " +str(dE) 
	
	#replace commons pragma
	pragmaList = ["#PATMO_photoRates","#PATMO_resolution"]
	replaceList = [allRates,res]

	patmo_string.fileReplaceBuild("src_f90/patmo_photoRates.f90", "build/patmo_photoRates.f90", \
		pragmaList, replaceList)
