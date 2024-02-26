#imports
from Bio import SeqIO
from shutil import copyfile
import os,sys,subprocess,re,math,HTSeq

sample_dict = {'XXXIdentifier1XXX':'XXXSample1XXX',....
	      .....'XXXIdentifierNXXX':'XXXSampleNXXX'}


seq_file_directory = 'XXXFastqsPathXXX'
os.chdir(seq_file_directory)
file_list = os.listdir()

# just rename files

filtered_file_list = [k for k in file_list if '.gz' not in k]
filtered_file_list = [k for k in filtered_file_list if 'Undetermined' not in k]
filtered_file_list = [k for k in filtered_file_list if '.fastq' in k]

for item in filtered_file_list:
    file = item
    result = re.search('3045(.*)-1', item)
    current_sample_index = result.group(1)
    current_sample = sample_dict[current_sample_index]
    new_name = current_sample + '.fastq'
    rename_command = 'mv {old} {new}'.format(old = file, new = new_name)
    os.system(rename_command)

def trimmer(fileDirectory):
    path = fileDirectory
    for filename in os.listdir(path):     ###for each file name in a specific director
        if filename.endswith('fastq') and 'trimmed' not in filename:   ###make this more specific----use ends with
            inputFile = os.path.join(path, filename)  ###joins path string and filename string
            outputFile = inputFile[:-5] + 'trimmed.fastq' ##explicitly find the dot from the end
            adapter = 'CTGTCTCTTATA'   ####
            function = 'cutadapt'
            minimum = '-m 10'
            cores = '--cores=8'
            separator = " "   #spelling
            cmdLnOutPutFile = outputFile + '.trimReport.txt'
            commandList = (function,'-a',adapter,minimum,cores,inputFile,'>',outputFile,'2>', cmdLnOutPutFile) 
            os.system(separator.join(commandList))

trimmer(seq_file_directory)

os.chdir(seq_file_directory)

starCommand = """for i in *.trimmed.fastq; do\n\
data_repository/STAR-2.7.1a/bin/Linux_x86_64_static/STAR --runMode alignReads --outSAMtype BAM SortedByCoordinate Unsorted --genomeDir XXXGenomePathXXX --readFilesIn $i --runThreadN 16 --outFileNamePrefix ${i%trimmed.fastq}\n\
done"""

subprocess.call(starCommand,shell=True)

os.chdir(seq_file_directory)

opener = 'htseq-count'

filtered_gtf_path = 'XXXPathXXX/genes.filter.gtf'

GTF = filtered_gtf_path

fileFormat = '-f bam'

stranded = '-s no'

idAttribute = '--idattr=gene_name'

fileList = os.listdir(seq_file_directory)

bamList = []

separator = " "

for name in fileList:
    if name.endswith('.Aligned.out.bam') == True:   #ignores the sorted bams
        bamList.append(name)

for fileName in bamList:
    outputFile = fileName[:-15]+'count'
    commandList = (opener, fileFormat, stranded,idAttribute,fileName,GTF,'>',outputFile)
    os.system(separator.join(commandList))
