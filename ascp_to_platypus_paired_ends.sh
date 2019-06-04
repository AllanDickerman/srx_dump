#Processes a mouse SRA file name through the programs ascp, parallel fastq-dump, fastq-mc || seqtk, STAR, velvet, opossum, and platypus
#outputs a log file containing information on program running times and file sizes
#Assumptions: all necessary programs to run are in $HOME 

centrifugeReads=1000000 #number of reads for centrifuge
pairedEnds=true #false == doesn't process the $base_2 file if it exists, false == does process the $base_2 file if it exists  
directoryBase="SRA_RNA_SEQ_Data"
seqtk=true #false == program runs fastq-mc, true == program runs seqtk
while getopts ":d:s"  opt; do
  case $opt in 
    s)
      pairedEnds=false
      ;;
    d)
      directoryBase=$OPTARG
      ;;
    \?)
      echo "invalid options" 2>&1
      echo "<optional arguments> -v [# of velvet reads] : -d [output directory] : -s (single end only)" 2>&1
      exit 1
      ;; 
  esac
done
shift $((OPTIND -1))
base=$1

#print information on run to stdout
echo "SRA ID: $base"
echo "output directory = $directoryBase" 
if [ "$pairedEnds" = true ] #conditional not working properly, check down the pipeline if others are 
  then 
    echo "Paired ends processing if applicable"
  else
    echo "Single ends processing only"
fi

#create working directory
mkdir $base
cd $base 

#creates output directory if it doesn't exist
if [ ! -d $HOME/$directoryBase ] 
  then
    mkdir $HOME/$directoryBase
    echo "created directory"
fi
if [ ! -d $HOME/$directoryBase/VCF_Files ]
  then 
    mkdir $HOME/$directoryBase/VCF_Files
    echo "created VCF_Files"
fi  
if [ ! -d $HOME/$directoryBase/ReadsPerGene_Files ] 
  then
    mkdir $HOME/$directoryBase/ReadsPerGene_Files
    echo "created ReadsPerGene_Files"
fi
if [ ! -d $HOME/$directoryBase/Centrifuge_Files ] 
  then
    mkdir $HOME/$directoryBase/Centrifuge_Files 
    echo "created Centrifuge_Files"
fi
if [ ! -d $HOME/$directoryBase/LogFiles ]
  then
    mkdir $HOME/$directoryBase/LogFiles
    echo "created LogFiles"
fi

#check if ID has already been processed by checking VCF file existence
if [ -f $HOME/$directoryBase/VCF_Files/$base.platypus.vcf ]
  then
    echo "$base VCF file exists"
    cd ..
    rm -r $base
    exit 1
fi 

#log file variables
log=./Log.$base.out #set log file name to variable "log"
single=true #assume single end read unless it goes through the paired end read option, outputs at bottom of log file
echo "Stage	Time_Start	Time_End	Size" > $log

#being pipeline, start with ascp download
if [ ! -f $base.sra ]
  then
    printf "ascp\t" >> $log
    date +"%s" | xargs echo -n >> $log
    printf "\t" >> $log
    $HOME/.aspera/cli/bin/ascp -i $HOME/.aspera/cli/etc/asperaweb_id_dsa.openssh -T anonftp@ftp.wip.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/${base:0:3}/${base:0:6}/$base/$base.sra $base.sra
    date +"%s" | xargs echo -n >> $log 
    printf "\t" >> $log 
    du -k $base.sra >> $log
fi
if [ ! -f $base.sra ]
  then
    echo "failed to download $base.sra"
    echo $base >> ../failedRuns.txt
    cd ..
    exit 1
fi

#parallel-fastq-dump
echo " " >> $log
printf "fastq-dump\t" >> $log
date +"%s" | xargs echo -n >> $log
printf "\t" >> $log
if [ ! -f ${base}_1.fastq ]
  then 
    time python $HOME/parallel-fastq-dump.py -s $base.sra -t 8 -X 30000000 --split-files
fi
if [ -f $base.fastq ]
  then 
    mv $base.fastq ${base}_1.fastq
fi
date +"%s" | xargs echo -n >> $log
printf "\t" >> $log
#only getting the size of $base_1.fastq
du -k ${base}_1.fastq >> $log

#trimm (fastq-mc || seqtk) && STAR
echo " " >> $log
printf "trim\t" >> $log
if [ "$pairedEnds" = true ] && [ -f ${base}_2.fastq ] #checks and modified
  then
    single=false #will print paired ends at bottom of log file
    echo "${base}_2.fastq exists: paired ends mode"
    date +"%s" | xargs echo -n >> $log
    printf "\t" >> $log
    if [ "$seqtk"  == true ]
      then 
        echo "running seqtk"
        seqtk mergepe ${base}_1.fastq ${base}_2.fastq | seqtk trimfq - | seqtk dropse > ${base}_12.trim.fastq
        seqtk seq -1 ${base}_12.trim.fastq > ${base}_1.trim.fastq
        seqtk seq -2 ${base}_12.trim.fastq > ${base}_2.trim.fastq	
      else
        echo "running fastq-mcf"
        time fastq-mcf -o ${base}_1.trim.fastq -o ${base}_2.trim.fastq $HOME/Trimmomatic/all_TruSeq_adapters.fa ${base}_1.fastq ${base}_2.fastq
    fi
    date +"%s" | xargs echo -n >> $log
    printf "\t"  >> $log
    du -k ${base}_1.trim.fastq >> $log
     
    echo " " >> $log
    printf "STAR\t" >> $log
    date +"%s" | xargs echo -n >> $log
    printf "\t" >> $log
    time STAR --runThreadN 10 --genomeDir $HOME/genome/index_for_STAR/ --readFilesIn ${base}_1.trim.fastq ${base}_2.trim.fastq --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --genomeLoad NoSharedMemory --outReadsUnmapped Fastx --outSAMattributes MD
    date +"%s" | xargs echo -n >> $log
    printf "\t" >> $log
    du -k ReadsPerGene.out.tab >> $log
else     
  echo "Single end mode"
  date +"%s" | xargs echo -n >> $log
  printf "\t" >> $log
  if [ "$seqtk" == true ]
    then
      echo "running seqtk"
      seqtk trimfq ${base}_1.fastq > ${base}_1.trim.fastq
    else
      echo "running fastq-mcf"
      time fastq-mcf -o ${base}_1.trim.fastq $HOME/Trimmomatic/all_TruSeq_adapters.fa ${base}_1.fastq 
  fi
  date +"%s" | xargs echo -n >> $log
  printf "\t" >> $log
  du -k ${base}_1.trim.fastq >> $log

  echo " " >> $log
  printf "STAR\t" >> $log
  date +"%s" | xargs echo -n >> $log
  printf "\t" >> $log
  time STAR --runThreadN 10 --genomeDir $HOME/genome/index_for_STAR/ --readFilesIn ${base}_1.trim.fastq --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --genomeLoad NoSharedMemory --outReadsUnmapped Fastx --outSAMattributes MD
  date +"%s" | xargs echo -n >> $log
  printf "\t" >> $log
  du -k ReadsPerGene.out.tab >> $log
fi
 
mv ReadsPerGene.out.tab $HOME/$directoryBase/ReadsPerGene_Files/$base.STAR.ReadsPerGene.out.tab 
mv Aligned.sortedByCoord.out.bam $base.Aligned.sortedByCoord.md.bam 

#centrifuge
echo " " >> $log
printf "centrifuge\t" >> $log
date +"%s" | xargs echo -n >> $log
printf "\t" >> $log
time $HOME/centrifuge/centrifuge --qupto $centrifugeReads -x ~/centrifuge/centrifuge_index/p_compressed+h+v -U Unmapped.out.mate1 -S /dev/null
date +"%s" | xargs echo -n >> $log
printf "\t" >> $log
du -k centrifuge_report.tsv >> $log
mv centrifuge_report.tsv $HOME/$directoryBase/Centrifuge_Files/${base}_centrifuge_out.txt
rm Unmapped.out.mate1

#opossum
echo " " >> $log
printf "opossum\t" >> $log
date +"%s" | xargs echo -n >> $log
printf "\t" >> $log
time python $HOME/Opossum/Opossum.py --BamFile=./$base.Aligned.sortedByCoord.md.bam --OutFile=./$base.opossum.md.bam --ProperlyPaired False
date +"%s" | xargs echo -n >> $log
printf "\t" >> $log
du -k $base.opossum.md.bam >> $log

#platypus
echo " " >> $log
printf "platypus\t" >> $log
date +"%s" | xargs echo -n >> $log
printf "\t" >> $log
time python $HOME/Platypus/Platypus.py callVariants --bamFiles=./$base.opossum.md.bam  --refFile=$HOME/genome/GRCm38.primary_assembly.genome.fa --output=$HOME/$directoryBase/VCF_Files/$base.platypus.vcf
date +"%s" | xargs echo -n >> $log
printf "\t" >> $log
du -k $HOME/$directoryBase/VCF_Files/$base.platypus.vcf >> $log 

#finish log file
echo " " >> $log
if [ "$single" == true ] 
  then 
    echo "Single End Run" >> $log
  else
    echo "Paired End Run" >> $log
fi
mv $log $HOME/$directoryBase/LogFiles

echo "Finished Processing $base"

cd ..
rm -r $base   
echo "done"
