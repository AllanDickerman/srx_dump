from pymongo import MongoClient
from pymongo import errors
import glob
import sys
import math
#import cPickle
ref_genome = "GRCm38.p5"		#Name of reference genome
program = "Platypus"			#Name of program that output VCF file 

#startup MongoDB: Assumes the server is already running
#try port 27020
client = MongoClient('localhost', 27018) #change this for running through slurm
db = client.mytestdb
vcf_collection = db.test_vcf_collection

#read info from file
v = "variants"

vcf_drc = "/home/clark96/MongoDB/test_VCF/*"
#vcf_drc = "/home/clark96/MongoDB/test_VCF/SRR3746168.sangerGeneSnps.vcf"
numPrints = 0
large_files = []
recordThreshold = 475000.0
for f in glob.glob(vcf_drc):
    if numPrints < 20:
        sys.stdout.write(".")
        numPrints = numPrints + 1
    else:
        sys.stdout.write(".\n")
        numPrints = 0
    sys.stdout.flush()
    expData = []
    recordCount = 0
    with open(f,"r") as fopen:
        for line in fopen:
            if line[0] == "#":
                continue			
            parsed = line.rstrip().split("\t")
            if parsed[6] == "PASS":	
                recordCount = recordCount+1
                data = []
                data.append(parsed[0])
                data.append(parsed[1])
                data.append(parsed[2])
                data.append(parsed[3])
                data.append(parsed[4])
                data.append(parsed[5])
                alleleIndices = parsed[9].split(":")[0]
                data.append(alleleIndices)	
                expData.append(data)
    numFiles_needed = math.ceil(recordCount/recordThreshold)	
    if numFiles_needed == 0:
        #input these files
        f_data = {}
        f_data["run_id"] = run_id
        f_data["referenceGenome"] = ref_genome
        f_data["program"] = program
        f_data["fileNum"] = 1
        f_data["totalFiles"] = numFiles_needed
        f_data[v] = []
        vcf_collection.insert_one(f_data)
        continue
    parts = f.split("/")
    
    currFile = 1 
    startIndex = 0
    endIndex = int(recordThreshold)
    run_id = parts[len(parts)-1].split(".")[0]
    #loop through and add necessary number of files	
    while currFile <= numFiles_needed:
        file_data = {}
        file_data["run_id"] = run_id 
        file_data["referenceGenome"] = ref_genome
        file_data["program"] = program
        file_data["fileNum"] = currFile
        file_data["totalFiles"] = numFiles_needed	
        
        #add data to file_data entry
        		
        file_data[v] = data[startIndex:endIndex]
        startIndex = endIndex
        endIndex = endIndex + int(recordThreshold)
        
        #add file data to mongo
        #print("inserting document")
        currFile = currFile + 1
        try:
            vcf_collection.insert_one(file_data)
        except errors.DocumentTooLarge as e:
            large_files.append(f)
