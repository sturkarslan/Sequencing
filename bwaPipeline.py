import glob, sys, os, string

# Process fastQ Files
def runPipeline():
    fastqFolders = glob.glob('%s/*/' %(dataDir))
    print
    print 'FASTQ Folders: %s' %fastqFolders
    
    for folder in fastqFolders:
        sampleFolder = folder.split("/")[1] # Folder containing reads
        #print sampleFolder
        fastqFilesFirst = glob.glob('%s/*R1*.fastq.gz' %(folder)) # 1st file the folder
        #print 'FASTQ Files: %s' %fastqFilesFirst
        
        for file in fastqFilesFirst:
            print
            print "\033[33m Processing %s \033[0m" %(file)
            fileName = file.split("_R1")[0]
            print fileName
            sampleResultsDir = resultsDir+ '/'+organism+'/'+fileName.split("/")[1]
            sampleTitle = fileName.split("/")[2] 
            lane = fileName.split("/")[2].split("_")[2]
            sampleName = fileName.split("/")[2].split("_")[1]
            #sampleId = fileName.split("/")[2].split("_")[0].split("-")[3] # uncomment for actual run
            sampleId = fileName.split("/")[2].split("_")[0] # uncomment for test run
            
            print "FileName: %s Lane: %s sample: %s ID: %s" %(fileName, lane, sampleName, sampleId)
            firstPair = fileName + "_R1_001.fastq.gz"
            secondPair = fileName + "_R2_001.fastq.gz"
            
            print "First Pair: %s, Second Pair: %s" %(firstPair, secondPair)
            
            # Run Fastqc
            runQC(firstPair, secondPair)
            # Run trimmomatic
            firstPairTrimmedPaired, secondPairTrimmedPaired = runTrim(firstPair, secondPair)
            # Run bwa    
            runBWA(firstPairTrimmedPaired,  secondPairTrimmedPaired, lane, sampleName, sampleId, sampleResultsDir, sampleTitle)
            # Run samtools fixmate
            runSamtoolsFixmate(sampleResultsDir,sampleTitle)
            # Run samtools sort
            runSamtoolsSort(sampleResultsDir,sampleTitle)
            # Run GATK Indel realignment
            runGATK(sampleResultsDir,sampleTitle)
        # Run markduplicates
        libraryName = runMarkDuplicates(sampleResultsDir,sampleTitle,sampleFolder)
        # Samtools flagstat
        flagstats(sampleResultsDir, libraryName)
        # Get BAM Statistics
        #bamStats(sampleResultsDir)
        # Run Samtools Variant calling
        variantCalling(sampleResultsDir,sampleTitle,libraryName)
        # Run variant calling with Varscan
        varscan(sampleResultsDir,sampleTitle,libraryName)
        # Variant calling with GATK HaploTypeCaller
        haplotypeCaller(sampleResultsDir,sampleTitle,libraryName)
        # Run SnpEff annotation
        runSnpEff(sampleResultsDir, libraryName)
        # run SnpSift Annotation
        runSnpSift(sampleResultsDir, libraryName)
            
        
    return firstPair, secondPair, libraryName, sampleResultsDir





# Quality control
def runQC(firstPair, secondPair):
    print
    print "\033[34m Running FastQC Quality Control \033[0m"
    
    # create results folder
    if not os.path.exists('%s' %(fastqcDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(fastqcDir)  
        os.makedirs('%s' %(fastqcDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(fastqcDir)
    # run Command and write output into both screen and logfile with 2>&1 | tee -a %s
    cmd = 'unbuffer %s -t 4 -o %s %s %s 2>&1 | tee -a %s' %(fastqc, fastqcDir, firstPair, secondPair, pipelineLog)
    print 'FastQC Command:', cmd
    #os.system(cmd)
    print
    print
  
  
  
      

# Trim reads
def runTrim(firstPair, secondPair):
    print
    print "\033[34m Running Read Trimming... \033[0m"
    # Program Parameters
    illuminaClip = "ILLUMINACLIP:/users/sturkars/Trimmomatic-0.35/adapters/NexteraPE-PE.fa:2:30:10" #Remove adapters 
    leading = "LEADING:3" #Remove leading low quality or N bases
    trailing = "TRAILING:3" #Remove trailing low quality or N bases
    slidingWindow = "SLIDINGWINDOW:4:20" #Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 20
    minLen = "MINLEN:36"
    paired = "PE" # or SE for single end
    threads = "-phred33 -threads 8"
      
    # define result files
    filesFolder = firstPair.split('/')[0] + "/" + firstPair.split('/')[1]
    firstPairTrimmedPaired = filesFolder+"/trimmed/"+firstPair.split('.fastq.gz')[0].split("/")[2] + "_paired_trimmed.fastq.gz"
    secondPairTrimmedPaired = filesFolder+"/trimmed/"+secondPair.split('.fastq.gz')[0].split("/")[2] + "_paired_trimmed.fastq.gz"
    firstPairTrimmedUnpaired = filesFolder+"/trimmed/"+firstPair.split('.fastq.gz')[0].split("/")[2] + "_unpaired_trimmed.fastq.gz"
    secondPairTrimmedUnpaired = filesFolder+"/trimmed/"+secondPair.split('.fastq.gz')[0].split("/")[2] + "_unpaired_trimmed.fastq.gz"
    trimDir = filesFolder+"/trimmed/"
    
    # create trim folder
    if not os.path.exists('%s' %(trimDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(trimDir)  
        os.makedirs('%s' %(trimDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(trimDir)
    
    # define command
    cmd = '/users/sturkars/java/bin/java -Xmx128m -jar %s %s %s %s %s %s %s %s %s %s %s %s %s %s 2>&1 | tee -a %s' %(trimmomaticPath, paired, threads, firstPair, secondPair, firstPairTrimmedPaired, firstPairTrimmedUnpaired, secondPairTrimmedPaired, secondPairTrimmedUnpaired, illuminaClip, leading, trailing, slidingWindow, minLen, pipelineLog)
    print "Trimmomatic Command: ", cmd
    #os.system(cmd)
    return firstPairTrimmedPaired, secondPairTrimmedPaired
    print
    print



    
def runBWA(firstPairTrimmedPaired, secondPairTrimmedPaired,lane,sampleName,sampleId,sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running BWA alignment... \033[0m"
    
    # create results folder
    if not os.path.exists('%s' %(sampleResultsDir)):
        print '\033[31m %s directory does NOT exists.Creating one! \033[0m' %(sampleResultsDir)  
        os.makedirs('%s' %(sampleResultsDir))
    else:
        print '\033[31m %s directory exists. Not creating. \033[0m' %(sampleResultsDir)
    
    # modify read group information
    readGroup = "@RG\\tID:%s\\tPL:ILLUMINA\\tSM:%s\\tLB:%s" %(sampleId, sampleName, lane)
    # bwa run command
    cmd = "%s mem -R '%s' %s %s %s > %s/%s.sam" %(bwaPath, readGroup, genomeFasta, firstPairTrimmedPaired, secondPairTrimmedPaired, sampleResultsDir, sampleTitle)
    print "Run BWA Command", cmd
    #os.system(cmd)
    print
 
  
  
      
def runSamtoolsFixmate(sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running SAMtools fixmate... \033[0m"
    # fixmate run command
    cmd = '%s fixmate -O bam %s/%s.sam %s/%s_fixmate.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, pipelineLog)
    print "Samtools Fixmate Command: ", cmd
    #os.system(cmd)



    
def runSamtoolsSort(sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running SAMtools sort.. \033[0m"
    # run command
    cmd = '%s sort -O bam -o %s/%s_sorted.bam -T tmp/%s_temp %s/%s_fixmate.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, sampleTitle, sampleTitle, sampleResultsDir, sampleTitle, pipelineLog)
    print "Samtools Sort Command: ", cmd
    # index bam file
    cmd2 = '%s index %s/%s_sorted.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, sampleTitle,  pipelineLog)
    #os.system(cmd)
    #os.system(cmd2)
 
 
   # 
    
def runGATK(sampleResultsDir,sampleTitle):
    print
    print "\033[34m Running GATK Realigner.. \033[0m"
    #run Target Interbal Creater command
    cmd1 = '%s -Xmx128m -jar %s -T RealignerTargetCreator -R %s -I %s/%s_sorted.bam -o %s/%s.intervals 2>&1 | tee -a %s' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, pipelineLog)
    #run Indel Realigner command
    cmd2 = '%s -Xmx4G -jar %s -T IndelRealigner -R %s -I %s/%s_sorted.bam -targetIntervals %s/%s.intervals -o %s/%s_realigned.bam 2>&1 | tee -a %s' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, pipelineLog)
    # index bam file
    cmd3 = '%s index %s/%s_realigned.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, sampleTitle, pipelineLog)
    # Detect covariates
    cmd4 = '%s -Xmx4G -jar %s -T BaseRecalibrator -R %s -knownSites %s -I %s/%s_realigned.bam -o %s/%s_recal.table 2>&1 | tee -a %s' %(javaPath, gatkPath, genomeFasta, knownSites, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, pipelineLog)
    # Adjust quality scores
    cmd5 = '%s -Xmx4G -jar %s -T PrintReads -R %s -I %s/%s_realigned.bam -BQSR %s/%s_recal.table -o %s/%s_recal.bam 2>&1 | tee -a %s' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, sampleResultsDir, sampleTitle, pipelineLog)
    
    print "Command GATK Interval Creater: ", cmd1
    #os.system(cmd1)
    print
    print "Command GATK Realigner: ", cmd2
    #os.system(cmd2)
    print
    print "Command GATK index BAM: ", cmd3
    #os.system(cmd3)
    print
    print "Command GATK BaseRecalibrator: ", cmd4
    #os.system(cmd4)
    print
    print "Command GATK PrintReads: ", cmd5
    #os.system(cmd5)
    print





    
def runMarkDuplicates(sampleResultsDir,sampleTitle,sampleFolder):
    print
    print "\033[34m Running Mark Duplicates.. \033[0m"
    libraryName = sampleTitle.split("_L")[0]
    # collect list of recalibrated bam files
    recalBams = glob.glob('results/%s/%s/*_recal.bam' %(organism,sampleFolder))
    print sampleFolder, recalBams
    metricsFile = sampleResultsDir+"/"+libraryName+'.metrics'
    # creaate command line parameter for each file
    bamList = []
    for i in recalBams:
        inputAdd = 'INPUT=%s' %(i)
        bamList.append(inputAdd)
    bamListJoined = " ".join(bamList)
    
    # MArk Duplicates
    cmd = '%s -Xmx4G -jar %s MarkDuplicates VALIDATION_STRINGENCY=LENIENT METRICS_FILE=%s %s OUTPUT=%s/%s_marked.bam 2>&1 | tee -a %s' %(javaPath, piccardPath, metricsFile, bamListJoined, sampleResultsDir, libraryName, pipelineLog)
    print "Mark Duplicated Command: ", cmd
    print
    # index bam file
    cmd2 = '%s index %s/%s_marked.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, libraryName, pipelineLog)
    print "Index bamfile Command: ", cmd2
    print 
    #os.system(cmd)
    #os.system(cmd2)            
    return libraryName




    
def variantCalling(sampleResultsDir,sampleTitle,libraryName): # With samtools
    print
    print "\033[34m Running SAMtools Variant Calling.. \033[0m"
    # Produce BCF file with all locations in the genome
    cmd = '%s mpileup -go %s/%s.bcf -f %s %s/%s_marked.bam 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, libraryName, genomeFasta, sampleResultsDir, libraryName, pipelineLog)
    # Reduce list of sites
    cmd2 = '%s call -vmO z -o %s/%s.vcf.gz %s/%s.bcf 2>&1 | tee -a %s' %(bcftoolsPath, sampleResultsDir, libraryName, sampleResultsDir, libraryName, pipelineLog)
    # Prepare vcf file for querying
    cmd3 = '%s -p vcf %s/%s.vcf.gz 2>&1 | tee -a %s' %(tabixPath, sampleResultsDir, libraryName, pipelineLog)
    #Filtering
    percentageString = "%"
    cmd4 = "%s filter -O z -o %s/%s.filtered.vcf.gz -s LOWQUAL -i '%sQUAL>10' %s/%s.vcf.gz 2>&1 | tee -a %s" %(bcftoolsPath, sampleResultsDir, libraryName, percentageString, sampleResultsDir, libraryName, pipelineLog)
    print "Variant Calling mpileup: ", cmd
    #os.system(cmd)
    print
    print "Variant Calling Reduce: ", cmd2
   # os.system(cmd2)
    print
    print "Variant Calling tabix: ", cmd3
    #os.system(cmd3)
    print
    print "Variant Calling filtering: ", cmd4
    #os.system(cmd4)
    





def varscan(sampleResultsDir,sampleTitle,libraryName): # with varscan
    print
    print "\033[34m Running Varscan.. \033[0m"
    # varscan for snps
    cmd0 = '%s mpileup -B -f %s -o %s/%s.pileup %s/%s_marked.bam' %(samtoolsPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # varscan for snps
    cmd = '%s -Xmx128m -jar %s mpileup2snp %s/%s.pileup --output-vcf 1 --pvalue 0.05 > %s/%s_varscan_snp.vcf' %(javaPath, varscanPath, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    # varscan for indels
    cmd2 = '%s -Xmx128m -jar %s mpileup2indel %s/%s.pileup --output-vcf 1 --pvalue 0.05 > %s/%s_varscan_indel.vcf' %(javaPath, varscanPath, sampleResultsDir, libraryName, sampleResultsDir, libraryName)
    print "samtools Mpileup: ", cmd0
    #os.system(cmd0)
    print "Varscan for SNPs: ", cmd
    #os.system(cmd)
    print
    print "Varscan for INDELS: ", cmd2
    #os.system(cmd2)
    
    
    
def haplotypeCaller(sampleResultsDir,libraryName): #with GATK HaploTypeCaller
    print
    print "\033[34m Running GATK Haplotype Variant Caller.. \033[0m"
    # haplotype command
    cmd1 = '%s -Xmx128m -jar %s -T HaplotypeCaller -R %s -I %s/%s_sorted.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o %s/%s_gatk-variants-raw.vcf 2>&1 | tee -a %s' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName, pipelineLog)
    # Select snp variants
    cmd2 = '%s -Xmx128m -jar %s -T SelectVariants -R %s -V %s/%s_gatk-variants-raw.vcf -selectType SNP -o %s/%s_gatk-variants-raw-snps.vcf 2>&1 | tee -a %s' %(javaPath, gatkPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName, pipelineLog)
    # Apply filters to SNPs
    cmd3 = "%s -Xmx128m -jar %s -T VariantFiltration -R %s -V %s/%s_gatk-variants-raw-snps.vcf --filterExpression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0' --filterName 'my_snp_filter' -o %s/%s_gatk-variants-filtered-snps.vcf 2>&1 | tee -a %s" %(javaPath, gatkPath, genomeFasta, sampleResultsDir, libraryName, sampleResultsDir, libraryName, pipelineLog)
    
    print "GATK HaplotypeCaller Comnand: ", cmd1
    os.system(cmd1)
    print "Select SNP Variants: ", cmd2
    os.system(cmd2)
    print
    print "Applying filters for SNPs: ", cmd3
    os.system(cmd3)    
    
    
    

def runSnpEff(sampleResultsDir, libraryName):
    print
    print "\033[34m Running SnpEff.. \033[0m"

    cmd1='%s -Xmx2g -jar %s -v -o gatk -s %s/%s_snps.annotated.summary.html %s %s/%s_varscan_snp.vcf > %s/%s_varscan_snp.annotated.vcf' %(javaPath, snpEffPath, sampleResultsDir, libraryName, snpEffDatabase, sampleResultsDir, libraryName,  sampleResultsDir, libraryName)

    cmd2='%s -Xmx2g -jar %s -v -o gatk -s %s/%s_indels.annotated.summary.html %s %s/%s_varscan_indel.vcf > %s/%s_varscan_indel.annotated.vcf' %(javaPath, snpEffPath, sampleResultsDir, libraryName, snpEffDatabase, sampleResultsDir, libraryName,  sampleResultsDir, libraryName)
    print cmd1
    print
    #os.system(cmd1)
    print
    #os.system(cmd2)



def runSnpSift(sampleResultsDir, libraryName):
    print
    print "\033[34m Running SnpSift.. \033[0m"

    path2script=snpEffPath.split('/snpEff.jar')[0]

    cmd1='cat %s/%s_varscan_snp.annotated.vcf | perl %s/scripts/vcfEffOnePerLine.pl | %s -Xmx128m -jar %s/SnpSift.jar extractFields - CHROM POS REF ALT AF AC DP MQ "(FILTER = \'PASS\')" "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "ANN[*].ERRORS" > %s/%s_varscan_snps_final.txt'%(sampleResultsDir, libraryName, path2script, javaPath, path2script, sampleResultsDir, libraryName)

    cmd2='cat %s/%s_varscan_indel.annotated.vcf | perl %s/scripts/vcfEffOnePerLine.pl | %s -Xmx128m -jar %s/SnpSift.jar extractFields - CHROM POS REF ALT AF AC DP MQ "(FILTER = \'PASS\')" "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "ANN[*].ERRORS" "EFF[*].AA_LEN" > %s/%s_varscan_indels_final.txt'%(sampleResultsDir, libraryName, path2script, javaPath, path2script, sampleResultsDir, libraryName)
    print cmd1
    print
    #os.system(cmd1)
    print
    print cmd2
    #os.system(cmd2)
    print

    return None






def bamStats(sampleResultsDir, sampleTitle):
    print
    print "\033[34m Running BamStats.. \033[0m"
    libraryName = sampleTitle.split("_L")[0]
    # MArk Duplicates
    cmd = '%s -Xmx4g -jar %s -i %s/%s_marked.bam -f %s/%s_genomic.bed -o %s/%s_stats.html -v html 2>&1 | tee -a %s' %(javaPath,bamstatsPath,sampleResultsDir, libraryName, genomeDir, organism, sampleResultsDir, libraryName, pipelineLog)
    print "BamStats Command: ", cmd
    print
    #os.system(cmd)
  
  
    
    
def flagstats(sampleResultsDir, libraryName):
    print
    print "\033[34m Running Samtools Flagstat.. \033[0m"
    # Samtools flagstat
    cmd = '%s flagstat %s/%s_marked.bam > %s/%s.flagstats 2>&1 | tee -a %s' %(samtoolsPath, sampleResultsDir, libraryName, sampleResultsDir, libraryName, pipelineLog)
    print "Flagstat Command: ", cmd
    print
    #os.system(cmd)    
   

        

def genomeIndexes():
    print
    print "\033[34m Creating Genome indexes... \033[0m"
    # indexing genomes
    cmd1=bwaPath+' index '+genomeFasta
    cmd2=samtoolsPath+' faidx '+genomeFasta
    cmd3='%s -Xmx4G -jar %s CreateSequenceDictionary R=%s O=%s.dict 2>&1 | tee -a %s' %(javaPath,piccardPath,genomeFasta,genomeFasta, pipelineLog)
    print "BWA Genome index command: ", cmd1
    print
    print "Samtools Genome index Command: ", cmd2
    print
    print "GATK CreateDictionary Command: ", cmd3

    #os.system(cmd1)
    #os.system(cmd2)
    #os.system(cmd3)

    
################# Main #################    
# Input files    
organism = "dvh"
dataDir = "data" 
genomeDir = "reference"
genomeFasta = '%s/%s_genomic.fna' %(genomeDir, organism)
genomeGff = '%s/%s_genomic.gff' %(genomeDir, organism)
resultsDir = "results"
fastqcDir = '%s/fastqc' %(resultsDir)
pipelineLog = '%s/pipelineLog.txt' %(resultsDir)
knownSites = '%s-variants-compiled.vcf' %(organism)
# snpEff databases
if organism == "mmp":
    snpEffDatabase = "Methanococcus_maripaludis_S2_uid58035"
if organism == "dvh":
    snpEffDatabase = "Desulfovibrio_vulgaris_Hildenborough_uid57645"    

# Programs
bwaPath = "/users/sturkars/bwa-0.7.12/bwa" # path to STAR executable
samtoolsPath = "/users/sturkars/samtools-1.2/bin/samtools"   # path to Trimmomatic executable
bcftoolsPath = "/users/sturkars/bcftools-1.1/bcftools"
fastqc = "/users/sturkars/FastQC/fastqc" # path to fastqc executable
javaPath = "/usr/bin/java"
piccardPath = "/users/sturkars/picard-tools-1.139/picard.jar"
trimmomaticPath = "/users/sturkars/Trimmomatic-0.35/trimmomatic-0.35.jar"
gatkPath = "/users/sturkars/gatk/GenomeAnalysisTK.jar"
tabixPath = "/users/sturkars/htslib-1.1/tabix"
varscanPath = "/users/sturkars/VarScan.v2.3.9.jar"
snpEffPath = "/users/sturkars/snpEff/snpEff.jar"

# Run functions
#genomeIndexes()
firstPair, secondPair, libraryName, sampleResultsDir = runPipeline()
    