from pymongo import MongoClient
import glob
import sys
import collections
import hashlib

def insertStarGeneOrder (collection,gene_order,genome_name,program,hashKey):
        new_record = {}
        new_record['program'] = program
        new_record['ref_genome'] = genome_name
        new_record['gene_order'] = gene_order
        new_record['hashKey'] = hashKey
        collection.insert_one(new_record)

ref_genome = "GRCm38.p5"		#Name of the reference genome
program = "STAR"			#Name of the program that output the ReadsPerGene info

#startup MongoDB: Assumes the server is already running 
client = MongoClient('localhost', 27018) #change this for running through slurm
db = client.mytestdb
order_collection = db.test_genomeOrder_collection
rpg_collection = db.test_rpg_collection

#read info from file
rpg_drc = "/home/clark96/MongoDB/test_RPG/*"
for f in glob.glob(rpg_drc):
    #print(f)
    file_data = {}
    gene_order = []
    parts = f.split("/")
    run_id = parts[len(parts)-1]
    file_data["run_id"] = run_id
    file_data["referenceGenome"] = ref_genome
    file_data["program"] = program
    file_data["geneCounts"] = []
    f_open = open(f,"r")
    for line in f_open:
        parsed = line.rstrip().split()
        gene_order.append(parsed[0])
        count = int(parsed[1])
        file_data["geneCounts"].append(count)
			
    #check gene order
    order_exists = False
    geneOrder_hash = hashlib.md5(",".join(gene_order)).hexdigest()
    for record in order_collection.find({"ref_genome" : ref_genome,"program" : program}):
        record_order = record["gene_order"]
        recOrder_hash = record["hashKey"]
        if geneOrder_hash == recOrder_hash:
            order_exists = True
            break	
    if not order_exists: #insert gene order
        insertStarGeneOrder(order_collection,gene_order,ref_genome,program,geneOrder_hash)

    rpg_collection.insert_one(file_data)
    f_open.close()
