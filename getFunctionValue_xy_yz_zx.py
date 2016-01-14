#!/usr/bin/env python
"""
Calculate the deformation under the given condition and return the deviation angle as a function value.

Fixed parameters: widening ration
Fitting parameters: xy shear ratio, xz shear ratio, zx shear ratio

Basic parameters should be prepared in the parameter file with a key-value pair format:
	"reductionRate": initial reduction rate
	"wideningRatio": widening per reduction
	"finalTime": finale time (s)
	"normalSteps": steps per time in normal, stable calculation region after yielding
	"yieldingTime": initial time for careful calculation region around yielding to avoid the conversion error
	"yieldingSteps": steps per time for careful calculation region around yielding to avoid the conversion error
	"dest_phi1": phi1 in euler angle for destination orientation
	"dest_Phi": Phi in euler angle for destination orientation
	"dest_phi2": phi2 in euler angle for destination orientation
	"initialShear1": [xy,yz,zx]
	"initialShear2": [xy,yz,zx]
	"initialShear3": [xy,yz,zx]
	"initialShear4": [xy,yz,zx]
"""
import numpy,os,damask,string,sys,subprocess,re,glob,ConfigParser,time
import deviationAngle

# initialize global parameters
# First set a kind of default, for the case that no initial file is given,
basename = "test"
reductionRate = 1.0e3
wideningRatio = 0.0
finalTime	  = 500
normalSteps	  = 2
yieldingTime  = 10
yieldingSteps = 100
dest_phi1 = 0.0
dest_Phi = 0.0
dest_phi2 = 0.0
initialShear1 = [0.0,0.0,0.0]
initialShear2 = [1.0,0.0,0.0]
initialShear3 = [0.0,1.0,0.0]
initialShear4 = [0.0,0.0,1.0]
initialFunctionValues = []
geomFile = "4x4x4.geom"

# Initialize with the initial file
initFile = "initial"	
conf = ConfigParser.SafeConfigParser()
conf.read(initFile)
basename = conf.get('sample', 'basename')
reductionRate = conf.getfloat('load', 'reductionRate')
wideningRatio = conf.getfloat('load', 'wideningRatio')
finalTime	  = conf.getfloat('load', 'finalTime')
normalSteps	  = conf.getint('load', 'normalSteps')
yieldingTime  = conf.getfloat('load', 'yieldingTime')
yieldingSteps = conf.getint('load', 'yieldingSteps')
dest_phi1 = conf.getfloat('destination', 'dest_phi1')
dest_Phi = conf.getfloat('destination', 'dest_Phi')
dest_phi2 = conf.getfloat('destination', 'dest_phi2')
initialShear1 = map(float,conf.get('initial', 'initialShear1').split(","))
initialShear2 = map(float,conf.get('initial', 'initialShear2').split(","))
initialShear3 = map(float,conf.get('initial', 'initialShear3').split(","))
initialShear4 = map(float,conf.get('initial', 'initialShear4').split(","))
initialXlist = [initialShear1,initialShear2,initialShear3,initialShear4]
flistStr = conf.get('initial', 'initialFunctionValues')
if not flistStr == "": initialFunctionValues = map(float,flistStr.split(","))
geomFile = conf.get('file', 'geomFile')

dest_widening = 1.0+wideningRatio*finalTime*reductionRate
#print initialShear

# Prepare the summary file
summaryFile=basename+".summary"
summaryIndex="file\txy\tyz\tzx\tdeviation_angle\twidening\tphi1_eulerangles\tPhi_eulerangles\tphi2_eulerangles\n"
f = open(summaryFile,"w")
f.write(summaryIndex)
f.close

# parameters for load statement
## common functions
def addSpace(text):
	tmp = " " + text + " "
	return tmp

def shearStr(shearValue):
	shearValue *= reductionRate
	return addSpace("{0:e}".format(shearValue))
	
## fdot
asterisk = addSpace("*")
zero = addSpace("0")
yyStr = addSpace("{0:e}".format(wideningRatio*reductionRate))
zzStr = addSpace("{0:e}".format(-reductionRate))

## yielding calculation
yieldingTimeStr = "time" + addSpace(str(yieldingTime))
yieldingIncs = int(yieldingTime*yieldingSteps)
yieldingIncsStr = "incs" + addSpace(str(yieldingIncs))
yieldingFreqStr =  "freq" + addSpace(str(yieldingIncs))

## normal calculation
normalTime = finalTime-yieldingTime
normalTimeStr = "time" + addSpace(str(normalTime))
normalIncs = int(normalTime*normalSteps)
normalIncsStr = "incs" + addSpace(str(normalIncs))
normalFreqStr =  "freq" + addSpace(str(normalIncs/10)) # assuming normalIncs is a multiple of 10

def MakeNewLoadFile(basename,loadStatement):
	""" Make new load file of the current iteration series.
		The name of load file is sequentially pointed, as follows:
			test0001, test0002, and so on.
		The last 4 digits show the number of the load files
		 having the same 'basename', i.e. 'test' in this example.
		The number is expressed in the hex form.
	"""
	key = "./"+basename + "*.load"
	loadFileList = glob.glob(key)
	#print loadFileList
	# remove extension, path, and basename
	loadFileList = [ x.replace(".load","") for x			# remove extension
					 in [ x.replace("./","") for x			# remove path
					   in [ x. replace(basename,"") for x	# remove basename
					     in loadFileList]]] 
	if loadFileList == []: loadFileList=['_0000']			# case for the first load
	SuffixNumberList = [ int(x.replace("_",""),16) for x in loadFileList ] # number is in hex form.
	newSuffix = "_{0:0>4}".format("{0:x}".format(max(SuffixNumberList) + 1)) # express with four digits
	newFile = basename + newSuffix + ".load"
	print "load file: {0}".format(newFile)
 	f = open(newFile,'w')
 	f.write(loadStatement)
 	f.close
	return newFile

def printProceeding():
	print ".",

def execDamask_spectral(loadFile,geomFile):
	debugMode = "off"			# "on" makes the output of damask printed
	exitCode = 2
	#
	while exitCode == 2:
		proc=subprocess.Popen(executable='/opt/petsc/gnu-4.9/bin/mpiexec^n 1 DAMASK_spectral',\
								args=['--load', loadFile, '--geom', geomFile],\
								stderr=subprocess.PIPE,stdout=subprocess.PIPE, bufsize=1)
		while proc.poll() is None:               # while process is running
			myLine = proc.stdout.readline()
			if debugMode == "on":
				if len(myLine)>1: print myLine[0:-1]   # print output without extra newline
			exitCode = proc.returncode
		#print exitCode
		err = proc.stderr.readlines()
		if not err == ['0\n']:
			print '-------------------------------------------------------'
			print 'error messages', err
			print '-------------------------------------------------------'
		if exitCode==2:
			print "The calculation is terminated with exit code of",exitCode
			break
	return exitCode

def execPostProcessing(outputFile):
	proc=subprocess.Popen(['postResults','--cr','f,eulerangles',outputFile])
	print "waiting for finishing postResults from",outputFile
	while proc.poll() is None:               # while process is running
		printProceeding
	return proc.returncode
	
def calcDeviationAngle(postFile):
	print
	print "#" * 10
	# read the result file
	devAng = -1
	widening = -1
	phi1 = -1
	Phi = -1
	phi2 = -1
	if os.path.exists(postFile):
		table=damask.ASCIItable(open(postFile))
		table.head_read()
		wideningIndex=table.labels.index('5_f')
		phi1Index=table.labels.index('1_eulerangles')
		PhiIndex=table.labels.index('2_eulerangles')
		phi2Index=table.labels.index('3_eulerangles')
		while table.data_read():
			widening=table.data[wideningIndex]
			phi1=table.data[phi1Index]
			Phi=table.data[PhiIndex]
			phi2=table.data[phi2Index]
		if numpy.testing.assert_approx_equal(widening,dest_widening) is None:
			devAng=deviationAngle.DeviationAngle(phi1,Phi,phi2,dest_phi1,dest_Phi,dest_phi2)
	else: print "No post file"
	result = {"devAng": devAng, "widening": widening, "eulerangles": [phi1,Phi,phi2]}
	return result


def f(xy,yz,zx):
	# Make load file
	xyStr = shearStr(float(xy))
	yzStr = shearStr(float(yz))
	zxStr = shearStr(float(zx))
	loadStatementBase = "fdot" + \
							asterisk + xyStr + zxStr + \
 							zero + yyStr + yzStr + \
 							zero + zero + zzStr + \
 						"stress" + \
 							zero + asterisk + asterisk + \
 							asterisk + asterisk + asterisk + \
 							asterisk + asterisk + asterisk
	yieldingStatement = loadStatementBase \
							+ yieldingTimeStr + yieldingIncsStr + yieldingFreqStr
	normalStatement = loadStatementBase \
							+ normalTimeStr + normalIncsStr + normalFreqStr
	loadStatement = yieldingStatement + "\n" + normalStatement
	loadFile = MakeNewLoadFile(basename,loadStatement)
 	proc=subprocess.Popen(['cat',loadFile],\
 				stderr=subprocess.PIPE,stdout=subprocess.PIPE)
	print 
	for loadCondition in iter(proc.stdout.readline,''):
		print loadCondition[0:-1] # without extra newline
	print 

	# Execute DAMASK_spectral
	exitCode = execDamask_spectral(loadFile,geomFile)
	
	# Post processing
	outputFileBase = geomFile.replace(".geom","") + \
					"_" + loadFile.replace(".load","")
	outputFile = outputFileBase + ".spectralOut"
	exitCode=execPostProcessing(outputFile)
	
	#  Calc deviation from the destination
	postFile = "postProc/"+outputFileBase+".txt"
	print "The post processed file is ",postFile
	result = calcDeviationAngle(postFile) # result = {devAng:,widening:,eulerangles:}
	
	# Delete spectral output files
	for outputFile in glob.glob(outputFileBase+"*"):
		os.remove(outputFile)

	# add to the summary file
	summary = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(loadFile,xy,yz,zx,\
		result["devAng"],result["widening"],result["eulerangles"][0],result["eulerangles"][1],result["eulerangles"][2])
	f = open(summaryFile,"a")
	f.write(summary)
	f.close
	print
	print "Summary of this calculation"
	print summary
	print "#" * 10
	print
	return result["devAng"]
	
# reload(deviationAngle)
# print dest_phi1,dest_Phi,dest_phi2
# print deviationAngle.DeviationAngle(120,0,0,dest_phi1,dest_Phi,dest_phi2)
#print f(0.1,0.1,0.1)
# postFile="postProc/4x4x4_test001b.txt"
# table=damask.ASCIItable(open(postFile))
# table.head_read()
# print table.labels
# print table.labels.index('5_f')
# while table.data_read():
# 	widening = table.data[12]
# print widening
