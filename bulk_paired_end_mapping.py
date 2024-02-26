import os, sys, subprocess

sample_dict = {'XXXIdentifier1XXX':'XXXSample1XXX',....
	      .....'XXXIdentifierNXXX':'XXXSampleNXXX'}


seq_file_directory = 'XXX'

os.chdir(seq_file_directory)
file_list = os.listdir()

#rename files

filtered_file_list = [k for k in file_list if '.gz' not in k]
filtered_file_list = [k for k in file_list if '.tsv' not in k]
filtered_file_list = [k for k in filtered_file_list if 'Undetermined' not in k]

sorted_pairs = [sorted(filtered_file_list)[i:i + 2] for i in range(0, len(sorted(filtered_file_list)), 2)]   

for item in sorted_pairs:
    first_file = item[0]
    second_file = item[1]
    result = re.search('3045(.*)-1', item[0])
    current_sample_index = result.group(1)
    current_sample = sample_dict[current_sample_index]
    read1_file = current_sample + '_R1' + '.r.fastq'
    read2_file = current_sample + '_R2' + '.r.fastq'
    separator = " "
    commandList = ('cp',first_file,read1_file)
    os.system(separator.join(commandList))
    commandList = ('cp',second_file,read2_file)
    os.system(separator.join(commandList))

#trim the reads

function = 'cutadapt'
adapter = 'CTGTCTCTTATACACATCT'
minimum = '-m 10'
cores = '--cores=8'
separator = " "

file_list = os.listdir()
filtered_file_list = [k for k in file_list if 'r.fastq' in k]
sorted_pairs = [sorted(filtered_file_list)[i:i + 2] for i in range(0, len(sorted(filtered_file_list)), 2)]   

for item in sorted_pairs:
    reads1 = item[0]
    reads2 = item[1]
    outFile1 = reads1[:-(len('r.fastq'))] + 'r.trimmed.fastq'
    outFile2 = reads2[:-(len('r.fastq'))] + 'r.trimmed.fastq'
    reportFile = item[0][0:4] + '.trimReport.txt'
    trimCommandList = (function,'-a',adapter,'-A',adapter,minimum,cores,'-o',outFile1,'-p',outFile2,reads1,reads2,'>', reportFile)
    trimCommand = separator.join(trimCommandList)
    subprocess.call(trimCommand,shell=True)

os.chdir(seq_file_directory)
file_list = os.listdir()

filtered_file_list = [k for k in file_list if 'r.trimmed.fastq' in k]
filtered_file_list = [k for k in filtered_file_list if '.txt' not in k]
sorted_pairs = [sorted(filtered_file_list)[i:i + 2] for i in range(0, len(sorted(filtered_file_list)), 2)]

#map the reads   

star_directory = 'datarepository/Linux_x86_64_static/STAR'
command_piece_1 = '--runMode alignReads --outSAMtype BAM SortedByCoordinate Unsorted --genomeDir XXXGenomePathXXX --readFilesIn'
command_piece_2 = '--runThreadN 16 --outFileNamePrefix'
separator = " "   #spelling

for item in sorted_pairs:
    reads1 = item[0]
    reads2 = item[1]
    outFileNamePrefix = reads1[:-(len('XXX.r.trimmed.fastq'))] + '.'
    starCommandList = (star_directory,command_piece_1,reads1,reads2,command_piece_2,outFileNamePrefix)
    starCommand = separator.join(starCommandList)
    subprocess.call(starCommand,shell=True)


# get the counts

import HTSeq
import os

os.chdir(seq_file_directory)

opener = 'htseq-count'

filtered_gtf_path = 'XXXGTFPathXXX'

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



