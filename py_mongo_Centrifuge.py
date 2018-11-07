from pymongo import MongoClient
import glob
import sys

ref_genome = "GRCm38.p5"
program = "Centrifuge"
num_reads_threshold = 100

s = "species"

#startup MongoDB: Assumes the server is already running
client = MongoClient('localhost', 27017) #change this for running through slurm
db = client.mytestdb
cent_collection = db.test_centrifuge_collection

#read info from file
cent_drc = "/home/clark96/MongoDB/test_Centrifuge/*"
for f in glob.glob(cent_drc):
    file_data = {}
    parts = f.split("/")
    #get run id from file name
    cent_file_name = parts[len(parts)-1]
    run_id = cent_file_name.split("_")[0]
    file_data["run_id"] = run_id
    file_data["referenceGenome"] = ref_genome
    file_data["program"] = program
    file_data[s] = {}	
    #add relevant info to file_data 
    f_open = open(f,"r")
    next(f_open)
    name = ""
    for line in f_open:
        parsed = line.rstrip().split("\t")	
        name = parsed[0].replace(".","")
        taxon_id = parsed[1]
        taxon_rank = parsed[2]
        num_reads = parsed[4] 
        if int(num_reads) >= num_reads_threshold:
            file_data[s][name] = {}
            file_data[s][name]["taxon_id"] = taxon_id
            file_data[s][name]["taxon_rank"] = taxon_rank
            file_data[s][name]["num_reads"] = num_reads
    f_open.close()
    cent_collection.insert_one(file_data)
