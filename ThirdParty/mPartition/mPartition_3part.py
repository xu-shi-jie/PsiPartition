#Partition an alignment in Phylip format to 3 subsets 
#Author: Thulek@gmail.com
#using these tools: IQ-TREE, TIGER
 
from os import listdir
from os.path import isfile, join
import os
import sys
import platform
import subprocess
import time
import array as arr
import config
import random
import numpy as np
import string
#import runCommand

import argparse

text = "===================\nSplit partitions from an alignments\n===================="

parser = argparse.ArgumentParser(description = text)  
#parser.parse_args()  

parser.add_argument("-f", "--filex", help="The alignment file")
parser.add_argument("-t", "--treefile", help="The tree file")
parser.add_argument("-tiger", "--tiger", help="The tiger file")
parser.add_argument("-m", "--maxlength", help="Max length of sequences")
parser.add_argument("-p", "--parfile", help="Partition file")
parser.add_argument("-tper", "--tper", help="Percent - minimum length")
parser.add_argument("-mset", "--mset", help="mset model, separated by commas")
parser.add_argument("-o", "--output", help="The output directory")
parser.add_argument("-s", "--strategy", help="The strategy to choose partition")
parser.add_argument("-inv", "--inv", help="The invariant site list file")
parser.add_argument("-prot", "--prot", help="protein: 1, DNA: 0")

# read arguments from the command line
args = parser.parse_args()

filex = args.filex
output = args.output
if args.treefile:
	treefile = args.treefile
else:
	treefile = "NOTTREE"
	
strategy = 2
if args.strategy:
	strategy = args.strategy

prot = 1
if args.prot:
	prot = int(args.prot)

tig_rate = "0"
if args.tiger:
	tig_rate = args.tiger
	
invF = "0"
if args.inv:
	invF = args.inv

#iqtree_path = "$IQTREE/" #IQ-TREE path including splash
#tiger_path = "/home/lkthu/Partitions/tiger_original/"

iqtree_path = config.iqtree_path
tiger_path = config.tiger_path

RANPROB = 0.1 #probability that a site is assigned to a group randomly

#mset = "LG,WAG"
mset = "JC69,HKY,GTR"
if args.mset:
	mset = args.mset

tper = 10.0
if args.tper:
	tper = float(args.tper)

if tper >= 50: 
	sys.exit("The minimum percent parameter must be less than 50.")

f = open(filex)
line = f.readline()
f.close()
currLength = 0
if " " in line:
	currLength = int(line.split(" ",1)[1].strip())
elif "\t" in line:
	currLength = int(line.split("\t",1)[1].strip())

if args.maxlength:
	maxlength = int(args.maxlength)
else:
	maxlength = currLength

if not os.path.isdir(output):
	cmd = 'mkdir '+output
	os.system(cmd)
else:
	print(output+" is exists")
	
	
def countLines(files):
	if(os.path.isfile(files)):
		with open(files) as f:
			return len(f.readlines())
	else:
		return 0

def convertPhylip2Fasta(filex1,tofile):
	infile = open(filex1,"r")
	outfile = open(tofile,"w")
	i = 1
	for line in infile:
		if i > 1:
			if len(line.strip()) >= 1:
				if "\t" in line:
					tx = line.split("\t",1)
					outfile.write(">"+tx[0].strip()+"\n")
					outfile.write(tx[1].strip()+"\n")
				elif " " in line:
					tx = line.split(" ",1)
					outfile.write(">"+tx[0].strip()+"\n")
					outfile.write(tx[1].strip()+"\n")
		i = i + 1
	outfile.close()
	infile.close()
	if(countLines(tofile)>0):
		print "Finish to create the fasta file."
	else:
		print "Create fasta file error."

def line_prepender(filename, line):
	with open(filename, 'r+') as f:
		content = f.read()
		f.seek(0, 0)
		f.write(line.rstrip('\r\n') + '\n' + content)

def removeGapsPhylip(filex):
	cfile = open(filex,"r")
	i = 0
	t = 0
	treefn = ""
	if "/" in filex:
		treef = filex.split("/")
		treefn = "ck"+treef[len(treef)-1].strip()
	else:
		treefn = "ck"+filex
	if os.path.isfile(output+"/"+treefn):
		os.system("rm "+output+"/"+treefn)
		
	rfile = open(output+"/"+treefn,"w")
	for line in cfile:
		if i == 0:
			if " " in line:
				tax = int(line.strip().split(" ",1)[0].strip())
				sitex = int(line.strip().split(" ",1)[1].strip())
			elif "\t" in line:
				tax = int(line.strip().split("\t",1)[0].strip())
				sitex = int(line.strip().split("\t",1)[1].strip())
		elif i >= 1:
			if " " in line:
				checkLine = line.split(" ",1)
				if checkLine[1].count("-") == len(checkLine[1].strip()) or checkLine[1].count("N") == len(checkLine[1].strip()):
					
					t += 1
				else:
					rfile.write(line)
					
			elif "\t" in line:
				checkLine = line.split("\t",1)
				if checkLine[1].count("-") == len(checkLine[1].strip()) or checkLine[1].count("N") == len(checkLine[1].strip()):
					
					t += 1
				else:
					rfile.write(line)
		i += 1
	cfile.close()
	rfile.close()
	if t == 0:
		print("This file is ok.")
		if os.path.isfile(output+"/"+treefn):
			os.system("rm "+output+"/"+treefn)
		#newfile = filex
	else:
		print("This file has "+str(t)+" seqs only gaps or missing data.")
		tax = tax - t
		line_prepender(output+"/"+treefn,str(tax)+" "+str(sitex))
		
def getStringIndx(string1,idx):
	returnStr = ""
	i = 1
	for line in string1.split("\n"):
		#print("line "+ str(i) + " "+str(len(line))+" "+str(idx))
		if(len(line.strip()) >= int(idx)):
			returnStr += line[int(idx)-1]+"\n"
		elif(len(line.strip())>0):
			returnStr += "error"+str(len(line.strip())) + " "+str(idx)
			#print(str(idx) + " "+str(len(line.strip()))) 
		i+=1
	return returnStr

def strIntersection(s1, s2):
	out = ""
	for c in s1:
		if c in s2 and not c in out:
			out += c
	return out

def isConst(string1,prot):
	if prot ==0:
		stringa = list(set(string1.upper()))
		bases = 'ACGT'
		ot = strIntersection(stringa,bases.upper())
		if len(ot) <= 1:
			return 1
		else:
			return 0
	else:
		stringa = list(set(string1.upper()))
		bases = string.ascii_uppercase[:26]
		result = ""
		result = strIntersection(stringa,bases.upper())
		if len(result) <= 1:
			return 1
		else:
			return 0

def getParIdx(logfile,lh1,lh2,lh3,stg):
	v_return = 1
	openLog = open(logfile,"a+")
	if stg == 1: #assign to best likelihood subset
		if lh1 > lh2 and lh1 > lh3:
			v_return = 1
		elif lh2 > lh3:
			v_return = 2
		else:
			v_return = 3
		openLog.write("\n"+str(stg)+"\t"+str(v_return)+"\t"+str(lh1)+"\t"+str(lh2)+"\t"+str(lh3))
	elif stg==2: #assign to subset based on probability distribution
		minLH = min(lh1,lh2,lh3)
		eps1 = lh1-minLH
		eps2 = lh2-minLH
		eps3 = lh3-minLH
		total = np.exp(eps1) + np.exp(eps2) + np.exp(eps3)
		
		prob1 = float(np.exp(eps1)/total)
		prob2 = float(np.exp(eps2)/total)
		randPrb = random.uniform(0, 1)
		if randPrb <= prob1:
			v_return = 1
		elif randPrb <= (prob1+prob2):
			v_return = 2
		else:
			v_return = 3
		openLog.write("\n"+str(stg)+"\t"+str(v_return)+"\t"+str(lh1)+"\t"+str(lh2)+"\t"+str(lh3)+"\t"+str(randPrb)+"\t"+str(prob1)+"\t"+str(prob2)+"\t"+str(1.0-prob1-prob2))
	else: #assign to best likelihood part with prob 1-RANPROB, random with RANPROB
		random.seed()
		#v_return = random.randint(1,3)
		
		randPrb = random.uniform(0,1)
		if randPrb <= RANPROB:
			random.seed()
			v_return = random.randint(1,3)
			openLog.write("\n"+str(stg)+"\t"+str(v_return)+"\t"+str(randPrb)+"\t"+str(RANPROB))
		else:
			if lh1 > lh2 and lh1 > lh3:
				v_return = 1
			elif lh2 > lh3:
				v_return = 2
			else:
				v_return = 3
			openLog.write("\n"+str(stg)+"\t"+str(v_return)+"\t"+str(randPrb)+"\t"+str(lh1)+"\t"+str(lh2)+"\t"+str(lh3))
	openLog.close()
	return v_return

min_rate = 1.0
max_rate = 0.0
#tg_rate = []

treefn = ""

if "/" in filex:
	treef = filex.split("/")
	treefn = treef[len(treef)-1].strip()
else:
	treefn = filex
path = filex.strip(treefn)

#print treefn

parfile = ""
if args.parfile:
	parfile = args.parfile
else:
	parfile = output+"/par."+treefn+""
	if(os.path.isfile(parfile)):
		os.system("rm "+parfile)

removeGapsPhylip(filex) #remove gaps from phylip
newfile = filex
if os.path.isfile(output+"/ck"+treefn):
	newfile = output+"/ck"+treefn
print(newfile)
checkstate = 0
rvalue = []
ivalue = []
min_rate = 0.0
max_rate = 0.0
if tig_rate == "0":
	convertPhylip2Fasta(newfile,output+"/"+treefn+".FASTA")
	if not os.path.isfile(output+"/rate_"+treefn):
		os.system(tiger_path+"tiger -in "+output+"/"+treefn+".FASTA -f s,r -rl "+output+"/rate_"+treefn) 

	os.system("rm "+output+"/"+treefn+".FASTA")
	
	#check inv site in original alignment
	#seq = list of sequence content in file (not include taxon name)
	inFile = open(filex,"r")
	i = 1
	numSite = 0
	seq = ""
	for line in inFile:
		if i ==1:
			if(len(line.strip())==0):
				#print "File format's error."
				sys.exit("File format's error.")
			else:
				if " " in line:
					numSite = int(line.split(" ",1)[1].strip())
				elif "\t" in line:
					numSite = int(line.split("\t")[1].strip())
		else:
			if(len(line.strip()) > 0):
				if " " in line:
					if len(line.split(" ",1)[1].strip()) <> numSite:
						#print "File format's error."
						sys.exit("File format's error.")
					else:
						seq += line.split(" ",1)[1].strip()+"\n"
				elif "\t" in line:
					if len(line.split("\t")[1].strip()) <> numSite:
						#saprint "File format's error"
						sys.exit("File format's error")
					else:
						seq += line.split("\t")[1].strip()+"\n"
		i += 1
	inFile.close()
	os.system("touch "+output+"/inv_"+treefn)
	#lines = seq.splitlines()
	invF = open(output+"/inv_"+treefn,"r+")
	for i in range (0,numSite):
		siteCi = getStringIndx(seq,i+1)
		ivalue.append(isConst(siteCi,prot))
		if isConst(siteCi,prot):
			invF.write("1\n")
		else:
			invF.write("0\n")
	invF.close()
	
	
	if os.path.isfile(output+"/rate_"+treefn):
		ratefile = open(output+"/rate_"+treefn,"r")
		for line in ratefile:
			rvalue.append(float(line.strip()))
		ratefile.close()
		checkstate = 1
	else:
		print "Rate file's not exists, please check the input file. \n"
else:
	#if tiger_rate file's not exists.
	if not os.path.isfile(output+"/"+tig_rate) or not os.path.isfile(parfile):
		convertPhylip2Fasta(newfile,output+"/"+treefn+".FASTA")
		os.system(tiger_path+"tiger -in "+output+"/"+treefn+".FASTA -f s,r -rl "+output+"/rate_"+treefn) 
		os.system("rm "+output+"/"+treefn+".FASTA")
		if os.path.isfile(output+"/rate_"+treefn):
			ratefile = open(output+"/rate_"+treefn,"r")
			for line in ratefile:
				rvalue.append(float(line.strip()))
			ratefile.close()
		else:
			print "Rate file's not exists, please check the input file. \n"
		if os.path.isfile(output+"/inv_"+treefn):
			invfile = open(output+"/inv_"+treefn,"r")
			for line in invfile:
				ivalue.append(int(line.strip()))
			invfile.close()
		else:
			print "Inv site list file's not exists, please check the input file. \n"
	else:
		ftiger = open(output+"/"+tig_rate,"r")
		lines=ftiger.readlines()
		ftiger.close()	
		finv = open(output+"/"+invF,"r")
		ilines = finv.readlines()
		finv.close()
		opar = open(parfile,"r")
		i = 0
		for ln in opar:
			if ln.split(";")[1].strip() in treefn:
				rvalue.append(float(lines[i].strip()))
				ivalue.append(int(ilines[i].strip()))
			i+=1	
		opar.close()
	checkstate = 1
	
par1 = []
par2 = []
par3 = []
if checkstate == 1:
	
	newList = sorted(rvalue)
	if int(round(len(rvalue)*1/100))-1 >= 50:
		min_rate = float(newList[int(round(len(rvalue)*1/100))-1])
		max_rate = float(newList[int(round(len(rvalue)*99/100))-1])
	else:
		min_rate = float(newList[0])
		max_rate = float(newList[int(round(len(rvalue)))-1])
	avg_rate = min_rate+(max_rate-min_rate)/3
	a2vg_rate = min_rate+(max_rate-min_rate)*2/3 
	
	print "avg rate: "+str(avg_rate)
	print("2 check max rate: "+str(max_rate) + "; min rate: "+ str(min_rate)+"; 1/3 rate: "+str(avg_rate)+" 2/3 rate: "+str(a2vg_rate)+"\n")
	par1 = []
	par2 = []
	par3 = []
	i=1
	#ratefile = open(output+"/rate_"+treefn,"r")
	for line in rvalue:
		if float(line) > a2vg_rate:
			par3.append(i)
			#print("par2 added seq: "+str(i))
		elif float(line) < avg_rate:
			par1.append(i)
			#print("par1 added seq: "+str(i))
		else:
			par2.append(i)
		i = i + 1
	#ratefile.close()
	if (len(par1)==0 or len(par2)==0 or len(par3)==0):
		par1 = []
		par2 = []
		par3 = []
		i=1
		#ratefile = open(output+"/rate_"+treefn,"r")
		for line in rvalue:
			if float(line) >= a2vg_rate and len(par3) <= len(rvalue)/3:
				par3.append(i)
				#print("par2 added seq: "+str(i))
			elif float(line) <= avg_rate and len(par1) <= len(rvalue)/3:
				par1.append(i)
				#print("par1 added seq: "+str(i))
			else:
				par2.append(i)
			i = i + 1
		#ratefile.close()
	
	parFile = open(output+"/S1_Par_"+treefn,"w")
	parFile.write("#nexus\nbegin sets;\n")
	line = "charset Par1 ="
	for p in par1:
		line += " "+str(p)
	line += ";"
	parFile.write(line+"\n")
	line = "charset Par2 ="
	for p in par2:
		line += " "+str(p)
	line += ";"
	parFile.write(line+"\n")
	line = "charset Par3 ="
	for p in par3:
		line += " "+str(p)
	line += ";"
	parFile.write(line+"\nend;")
	parFile.close()
	command = iqtree_path+"iqtree -s "+newfile+" -m MFP -fast -mset "+mset+" -spp "+output+"/S1_Par_"+treefn+" -safe -pre "+output+"/"+treefn+" -seed 0 -redo \n"
	os.system(command)
	if not os.path.isfile(output+"/"+treefn+"_AICc.iqtree"):
		command = iqtree_path+"iqtree -s "+newfile+" -m MFP -fast -mset "+mset+" -safe -pre "+output+"/"+treefn+"_AICc  -seed 0  -redo \n"
		os.system(command)
	
model =""
first_model = ""
second_model = ""
third_model = ""
if(os.path.isfile(output+"/"+treefn+".iqtree")):
	gtModel = open(output+"/"+treefn+".iqtree","r")
	for l in gtModel:
		if "Best-fit model according to BIC: " in l:
			model = l.split("Best-fit model according to BIC: ")[1].strip()
			fmodel = model.split(":")
			first_model = fmodel[0].strip()
			second_model = fmodel[1].split(",")[1].strip()
			third_model = fmodel[2].split(",")[1].strip()
	gtModel.close()
first_model = first_model.replace("+ASC","")
second_model = second_model.replace("+ASC","")
third_model = third_model.replace("+ASC","")
print("First Best Model: "+first_model+" | Treefile: "+output+"/"+treefn+".treefile\n")
print("Second Best Model: "+second_model+" | Treefile: "+output+"/"+treefn+".treefile\n")
print("Third Best Model: "+third_model+" | Treefile: "+output+"/"+treefn+".treefile\n")

if treefile == "NOTTREE":
	os.system("cp "+output+"/"+treefn+"_AICc.treefile "+output+"/"+treefn+".tree")
	treefile = output+"/"+treefn+".tree"


def getAIC(fromFile):
	text = ""
	with open(fromFile) as infile:
		for line in infile:
			if "Akaike information criterion (AIC) score" in line.strip():
				text += line.split(" ")[5]
				return text

def getAICc(fromFile):
	text = ""
	with open(fromFile) as infile:
		for line in infile:
			if "Corrected Akaike information criterion (AICc) score:" in line.strip():
				text += line.split(" ")[6]
				return text

def getBIC(fromFile):
	text = ""
	with open(fromFile) as infile:
		for line in infile:
			if "Bayesian information criterion (BIC) score:" in line.strip():
				text += line.split(" ")[5]
				return text

#get n
def getN(fromFile):
	text = ""
	with open(fromFile) as infile:
		for line in infile:
			if "Input data" in line.strip():
				if line.count(" ")<8:
					text += line.split(" ")[5]
				else:
					text += line.split(" ")[8]
				#print(text)
				return text

#get free para
def getK(fromFile):
	text = ""
	with open(fromFile) as infile:
		for line in infile:
			if "free parameters" in line.strip():
				text += line.split(" ")[8]
				return text

command = iqtree_path+"iqtree -s "+newfile+" -m "+first_model+" -te "+treefile+" -wsl -safe -pre "+output+"/"+treefn+"_G1  -seed 0 -redo \n"
os.system(command)
command = iqtree_path+"iqtree -s "+newfile+" -m "+second_model+" -te "+treefile+" -wsl -safe -pre "+output+"/"+treefn+"_G2  -seed 0 -redo \n"
os.system(command)
command = iqtree_path+"iqtree -s "+newfile+" -m "+third_model+" -te "+treefile+" -wsl -safe -pre "+output+"/"+treefn+"_G3  -seed 0 -redo \n"
os.system(command)

g1 = []
g2 = []
g3 = []
if(os.path.isfile(output+"/"+treefn+"_G1.sitelh") and os.path.isfile(output+"/"+treefn+"_G2.sitelh") and os.path.isfile(output+"/"+treefn+"_G3.sitelh")):
	if(os.path.isfile(output+"/siteLH/"+treefn+".sitelh")):
		os.system("rm "+output+"/siteLH/"+treefn+".sitelh")
	rwSite = open(output+"/siteLH/"+treefn+".sitelh","w+")
	getG1 = open(output+"/"+treefn+"_G1.sitelh","r")
	i = 1
	for l in getG1:
		if i > 1:
			rwSite.write(l.replace("Site_Lh","G1"))
			lh = l.split()
			z = 1
			for x in lh:
				if z > 1:
					g1.append(float(x))
					#print(x+"\n")
				z+=1
		i = i + 1
	getG1.close()
	getG2 = open(output+"/"+treefn+"_G2.sitelh","r")
	i = 1
	for l in getG2:
		if i > 1:
			rwSite.write(l.replace("Site_Lh","G2"))
			lh = l.split()
			z = 1
			for x in lh:
				if z > 1:
					g2.append(float(x))
				z+=1
		i = i + 1
	getG2.close()
	getG3 = open(output+"/"+treefn+"_G3.sitelh","r")
	i = 1
	for l in getG3:
		if i > 1:
			rwSite.write(l.replace("Site_Lh","G3"))
			lh = l.split()
			z = 1
			for x in lh:
				if z > 1:
					g3.append(float(x))
				z+=1
		i = i + 1
	getG3.close()
	rwSite.close()
	
	if(os.path.isfile(output+"/siteLH/"+treefn+".parsitelh")):
		os.system("rm "+output+"/siteLH/"+treefn+".parsitelh")
	
	par1 = []
	par2 = []
	par3 = []
	i = 0
	while i < len(g1):
		#print(ivalue)
		if (ivalue[i]==1):
			pIdx = getParIdx(output+"/siteLH/"+treefn+".parsitelh",g1[i],g2[i],g3[i],2)
		else:
			pIdx = getParIdx(output+"/siteLH/"+treefn+".parsitelh",g1[i],g2[i],g3[i],1)
		if pIdx == 1:
			par1.append(i+1)
		elif pIdx == 2:
			par2.append(i+1)
		else: 
			par3.append(i+1)
		i += 1
	print("===============PARTITIONS LENGTH AFTER USING STRATEGY===================\n")
	print("par1: "+str(len(par1))+" | par2: "+str(len(par2))+" | par3: "+str(len(par3)))
	
	isOk = 1
	zlen = 0
	if len(par1) < int((maxlength)*tper/100):
		zlen += 1
	
	if len(par2) < int((maxlength)*tper/100):
		zlen += 1
	
	if len(par3) < int((maxlength)*tper/100):
		zlen += 1
	
	#Not Ok if 2 partitions are very small.
	if zlen >= 2:
		isOk = 0
	else:
		#if 1 partition is very small, assign it's elements to other partitions
		if len(par1) > 0 and len(par1) < int((maxlength)*tper/100):
			for pIx in par1:
				if g2[pIx-1] < g3[pIx-1]:
					par2.append(pIx)
				else:
					par3.append(pIx)
			par1 = []
			par2.sort()
			par3.sort() 
			print("Moving the elements of partition 1 to partition 2 and partition 3")
		elif len(par2) > 0 and len(par2) < int((maxlength)*tper/100):
			for pIx in par2:
				if g1[pIx-1] < g3[pIx-1]:
					par1.append(pIx)
				else:
					par3.append(pIx)
			par2 = []
			par1.sort()
			par3.sort()
			print("Moving the elements of partition 2 to partition 1 and partition 3")
		elif len(par3) > 0 and len(par3) < int((maxlength)*tper/100):
			for pIx in par3:
				if g1[pIx-1] < g2[pIx-1]:
					par1.append(pIx)
				else:
					par2.append(pIx)
			par3 = []
			par1.sort()
			par2.sort()
			print("Moving the elements of partition 3 to partition 2 and partition 1")
		print("===============PARTITIONS LENGTH AFTER USING RE-ASIGNED===================\n")
		print("par1: "+str(len(par1))+" | par2: "+str(len(par2))+" | par3: "+str(len(par3)))
	if isOk == 0:
		print("Can not split. Length of 3 partitions: "+str(len(par1))+"|"+str(len(par2))+"|"+str(len(par3)))
	#if (len(par1) < int((maxlength)*tper/100) and len(par1)>0) or (len(par2)>0 and len(par2) < int((maxlength)*tper/100)) or (len(par3)>0 and len(par3) < int((maxlength)*tper/100)):
	#	print "Not good. Partition's length is very small."
	#elif (len(par1) + len(par2) ==0) or (len(par1) + len(par3) ==0) or (len(par3) + len(par2) ==0):
	#	print "Don't need split."
	else:
		parFile = open(output+"/F1_Par_"+treefn,"w")
		if len(par1)>0:
			line = "P1 ="
			for p in par1:
				line += " "+str(p)+","
			line += ""
			parFile.write(line.strip(",")+"\n")
		if len(par2)>0:
			line = "P2 ="
			for p in par2:
				line += " "+str(p)+","
			line += ""
			parFile.write(line.strip(",")+"\n")
		if len(par3)>0:
			line = "P3 ="
			for p in par3:
				line += " "+str(p)+","
			line += ""
			parFile.write(line.strip(",")+"\n")
		parFile.close()
		#os.system("python extractPartitions.py -f "+filex+" -p "+output+"/F1_Par_"+treefn+" -o "+output)
		os.system("python splitPartition.py -f "+newfile+" -p "+output+"/F1_Par_"+treefn+"")
		
		if(os.path.isfile(path+""+treefn+"P1") and not os.path.isfile(output+"/"+treefn+"P1")):
			os.system("cp "+path+""+treefn+"P1 "+output+"/"+treefn+"P1")
			os.system("rm "+path+""+treefn+"P1")
		if(os.path.isfile(path+""+treefn+"P2") and not os.path.isfile(output+"/"+treefn+"P2")):
			os.system("cp "+path+""+treefn+"P2 "+output+"/"+treefn+"P2")
			os.system("rm "+path+""+treefn+"P2")
		if(os.path.isfile(path+""+treefn+"P3") and not os.path.isfile(output+"/"+treefn+"P3")):
			os.system("cp "+path+""+treefn+"P3 "+output+"/"+treefn+"P3")
			os.system("rm "+path+""+treefn+"P3")
		
		if(os.path.isfile(output+"/ck"+treefn+"P1") and not os.path.isfile(output+"/"+treefn+"P1")):
			os.system("cp "+output+"/ck"+treefn+"P1 "+output+"/"+treefn+"P1")
			os.system("rm "+output+"/ck"+treefn+"P1")
		if(os.path.isfile(output+"/ck"+treefn+"P2") and not os.path.isfile(output+"/"+treefn+"P2")):
			os.system("cp "+output+"/ck"+treefn+"P2 "+output+"/"+treefn+"P2")
			os.system("rm "+output+"/ck"+treefn+"P2")
		if(os.path.isfile(output+"/ck"+treefn+"P3") and not os.path.isfile(output+"/"+treefn+"P3")):
			os.system("cp "+output+"/ck"+treefn+"P3 "+output+"/"+treefn+"P3")
			os.system("rm "+output+"/ck"+treefn+"P3")
		
		
		#os.system("mv "+output+"/"+treefn+"_Par1 "+output+"/"+treefn+"P1")
		#os.system("mv "+output+"/"+treefn+"_Par2 "+output+"/"+treefn+"P2")
		command = ""
		if(os.path.isfile(output+"/"+treefn+"P1")):
			command = iqtree_path+"iqtree -s "+output+"/"+treefn+"P1 -m MFP -fast -mset "+mset+" -t "+treefile+" -safe  -seed 0  -pre "+output+"/"+treefn+"P1_AICc\n"
			os.system(command)
		if(os.path.isfile(output+"/"+treefn+"P2")):
			command = iqtree_path+"iqtree -s "+output+"/"+treefn+"P2 -m MFP -fast -mset "+mset+" -t "+treefile+" -safe  -seed 0 -pre "+output+"/"+treefn+"P2_AICc\n"
			os.system(command)
		if(os.path.isfile(output+"/"+treefn+"P3")):
			command = iqtree_path+"iqtree -s "+output+"/"+treefn+"P3 -m MFP -fast -mset "+mset+" -t "+treefile+" -safe  -seed 0 -pre "+output+"/"+treefn+"P3_AICc\n"
			os.system(command)
		
		# aic = 0.0
		# k = 0.0
		# n = 0.0
		# files = [f for f in listdir(output)]
		# for f in files:
			# if (treefn+"P1_AICc.iqtree" in f) or (treefn+"P2_AICc.iqtree" in f)  or (treefn+"P3_AICc.iqtree" in f):
				# n += float(getN(output+"/"+f))
				# k += float(getK(output+"/"+f))
				# aic += float(getAIC(output+"/"+f))
		# AICc = aic + 2*k*(k+1)/(n-k-1)
		# FirstAICc = 0
		# if(os.path.isfile(output+"/"+treefn+"_AICc.iqtree")):
			# FirstAICc = float(getAICc(output+"/"+treefn+"_AICc.iqtree"))
		# print "FirstAICc: "+str(FirstAICc)+" | New AICc: "+str(AICc)
		
		bic = 0.0
		files = [f for f in listdir(output)]
		for f in files:
			if (treefn+"P1_AICc.iqtree" in f) or (treefn+"P2_AICc.iqtree" in f)  or (treefn+"P3_AICc.iqtree" in f):
				#n += float(getN(output+"/"+f))
				#k += float(getK(output+"/"+f))
				bic += float(getBIC(output+"/"+f))
		#BIC = aic + 2*k*(k+1)/(n-k-1)
		FirstBIC = 0
		if(os.path.isfile(output+"/"+treefn+"_AICc.iqtree")):
			FirstBIC = float(getBIC(output+"/"+treefn+"_AICc.iqtree"))
		print "FirstAICc: "+str(FirstBIC)+" | New AICc: "+str(bic)
		
		checkAICcFile = 1
		if ((not os.path.isfile(output+"/"+treefn+"P1_AICc.iqtree")) and os.path.isfile(output+"/"+treefn+"P1")):
			checkAICcFile = 0
			print "IQ-TREE core dump ERROR - "+treefn+"P1_AICc.iqtree"
		
		if ((not os.path.isfile(output+"/"+treefn+"P2_AICc.iqtree")) and os.path.isfile(output+"/"+treefn+"P2")):
			checkAICcFile = 0
			print "IQ-TREE core dump ERROR - "+treefn+"P2_AICc.iqtree"
		
		if ((not os.path.isfile(output+"/"+treefn+"P3_AICc.iqtree")) and os.path.isfile(output+"/"+treefn+"P3")):
			checkAICcFile = 0
			print "IQ-TREE core dump ERROR - "+treefn+"P3_AICc.iqtree"
				
		if (not os.path.isfile(output+"/"+treefn+"_AICc.iqtree")):
			checkAICcFile = 0
			print "IQ-TREE core dump ERROR - "+treefn+"_AICc.iqtree"
			
		if checkAICcFile == 0:
			print "IQ-TREE core dump ERROR."
			os.system("cp "+output+"/"+treefn+"P*.log logs/")
			os.system("rm "+output+"/"+treefn+"P*")
			os.system("rm core.*")
			os.system("rm "+output+"/core.*")
		elif bic > FirstBIC:
			print "Worser."
			os.system("cp "+output+"/"+treefn+"P*.log logs/")			
			os.system("rm "+output+"/"+treefn+"P*")
			
			#os.system("rm "+output+"/"+treefn+"P2*")
		else:
			if(os.path.isfile(output+"/"+treefn+"P3") and os.path.isfile(output+"/"+treefn+"P2") and os.path.isfile(output+"/"+treefn+"P1")):
				print "Better. Extract "+ filex+ " to 3 partitions: "+output+"/"+treefn+"P[1,2,3]."
			else:
				print "Better. Extract "+ filex+ " to 2 partitions: "+output+"/"+treefn+"P[1,2]."
			if(os.path.isfile(parfile)):
				opar = open(parfile,"r")
				ipar = open(parfile+"_1","w")
				i = 1
				z = 1
				for line in opar:
					if line.split(";")[1].strip() in treefn:
						if z in par1:
							ipar.write(str(i)+";"+treefn+"P1\n")
						elif z in par2:
							ipar.write(str(i)+";"+treefn+"P2\n")
						elif z in par3:
							ipar.write(str(i)+";"+treefn+"P3\n")
						#if g1[z-1] > g2[z-1] and g1[z-1]>g3[z-1]:
						##if g1[z-1] > g2[z-1]:
						#	ipar.write(str(i)+";"+treefn+"P1\n")
						##else:
						#elif g2[z-1]>g3[z-1]:
						#	ipar.write(str(i)+";"+treefn+"P2\n")
						#else:
						#	ipar.write(str(i)+";"+treefn+"P3\n")
						z = z+1
					else:
						ipar.write(line.strip()+"\n")
					i += 1
				ipar.close()
				opar.close()
				os.system("mv "+parfile+"_1 "+parfile)
			else:
				opar = open(parfile,"w")
				i = 0
				while i < len(g1):
					if i+1 in par1:
						opar.write(str(i+1)+";"+treefn+"P1\n")
					elif i+1 in par2:
						opar.write(str(i+1)+";"+treefn+"P2\n")
					elif i+1 in par3:
						opar.write(str(i+1)+";"+treefn+"P3\n")
					#if g1[i] > g2[i] and g1[i]>g3[i]:
					#	opar.write(str(i+1)+";"+treefn+"P1\n")
					#elif g2[i]>g3[i]:
					#	opar.write(str(i+1)+";"+treefn+"P2\n")
					#else:
					#	opar.write(str(i+1)+";"+treefn+"P3\n")
					i += 1
				opar.close()
			#os.system("rm "+output+"/"+treefn+"P*AICc*")

		
