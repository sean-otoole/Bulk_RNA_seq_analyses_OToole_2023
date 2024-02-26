#imports
from Bio import SeqIO
from shutil import copyfile
import os,sys,subprocess,re,math,HTSeq


# vector files
agmat_vector = 'XXXpath/vector.gbXXX'
adamts2_vector = 'XXXpath/vector.gbXXX'
baz1a_vector = 'XXXpath/vector.gbXXX'

# vector info format: TSS, chromosome, strand, chromosome name

agmat_info = [141473983,'4','+','TPD36',agmat_vector]
adamts2_info = [50492911,'11','+','TPD130',adamts2_vector]
baz1a_info = [55033292,'12','-','FLG47',baz1a_vector]

vector_info_list = [agmat_info,adamts2_info,baz1a_info]

#original files
original_fasta = 'XXXPathXXX/Mus_musculus.GRCm39.dna.primary_assembly.fa'
original_gtf = 'XXXPathXXX/Mus_musculus.GRCm39.104.gtf'

#modified files
mod_fasta_path = 'XXXPathXXX/fasta.fa'
mod_gtf_path = 'XXXPathXXX/genes.gtf'
edited_gtf_path = 'XXXPathXXX/genes.edited.gtf'

filtered_gtf_path = 'XXXPathXXX/genes.filter.gtf'
bed_path = 'XXXPathXXX/genes.filter.bed'
genome_name = 'artificial_promoters_genome'

# ### amend the fasta file

i = 0

for vector in vector_info_list:
    # get the plasmid
    gb_file = vector[4]
    for gb_record in SeqIO.parse(open(gb_file,"r"), "genbank") :
        # now do something with the record
        print("Name %s, %i features" % (gb_record.name, len(gb_record.features)))
        print(repr(gb_record.seq))    
        vector_raw_seq = str(getattr(gb_record,'seq'))

    #use the original fasta to build new headers
    with open (original_fasta, "r") as myfile:
        genome_list=myfile.readlines()
    
    matched_items = [item for item in genome_list if ">" in item] #list comprehension to find the headers
    template_header = matched_items[20]
    virus_header = template_header.replace('X',vector[3])
    virus_header = virus_header.replace('171031299', str(len(vector_raw_seq)))

    ### break the plasmid string into 61 character increments for each line
    def chunkstring(string, length):
        return (string[0+i:length+i] for i in range(0, len(string), length))
    broken_dna = list(chunkstring(vector_raw_seq, 61))  ### this is where we've stored our list
    broken_dna = ['{0}\n'.format(element) for element in broken_dna]  #add new line characters

    #now place the plasmid sequences/lines sequence after the header
    fasta_input_list = []
    fasta_input_list.extend(broken_dna) 
    fasta_input_list.insert(0, virus_header)

    #copy the original file during the first iteration
    if i < 1:
        copyfile(original_fasta, mod_fasta_path)
    else: 
        pass
    
    i = i + 1
    
    #write to that file
    with open(mod_fasta_path, "a") as myfile:
        for item in fasta_input_list:
            myfile.write(item)

copyfile(original_fasta, mod_fasta_path)

#write Ns where the original sequence surrounding the TSS sites in the genome

modified_genome_file = open(mod_fasta_path, "r")
list_of_lines = modified_genome_file.readlines()

chrom_4_indices = [i for i, s in enumerate(list_of_lines) if '>4 dna' in s][0] + 1
chrom_11_indices = [i for i, s in enumerate(list_of_lines) if '>11 dna' in s][0] + 1
chrom_12_indices = [i for i, s in enumerate(list_of_lines) if '>12 dna' in s][0] + 1

chrom_indices_list = [chrom_4_indices,chrom_11_indices,chrom_12_indices]

for current_gene in [0,1,2]:
    starting_index = (vector_info_list[current_gene][0]-1000)/60
    starting_base = round((starting_index % 1) * 60) - 1
    starting_index = math.floor(starting_index) + chrom_indices_list[current_gene]

    ending_index = (vector_info_list[current_gene][0]+1000)/60
    ending_base = round((ending_index % 1) * 60) -1 
    ending_index = math.floor(ending_index) + chrom_indices_list[current_gene]

    # operation on first index

    string_to_replace = list_of_lines[starting_index][starting_base:]
    retaining_string = list_of_lines[starting_index][:starting_base]
    replacement_string  = re.sub('[A-Z]','N',string_to_replace)
    new_string = retaining_string + replacement_string
    list_of_lines[starting_index] = new_string

    # operation on middle indices
    
    middle_indices = list(range((starting_index+1), (ending_index)))
    
    for index in middle_indices:
        new_string  = re.sub('[A-Z]','N',list_of_lines[index])
        list_of_lines[index] = new_string

    # operation on last index
    string_to_replace = list_of_lines[ending_index][:ending_base]
    retaining_string = list_of_lines[ending_index][ending_base:]
    replacement_string  = re.sub('[A-Z]','N',string_to_replace)
    new_string = replacement_string + retaining_string
    list_of_lines[ending_index] = new_string

modified_genome_file = open(mod_fasta_path, "w")
modified_genome_file.writelines(list_of_lines)
modified_genome_file.close()


copyfile(original_gtf, mod_gtf_path)

copyfile(original_gtf, mod_gtf_path)

TPD36 = ['TPD36','TPD36','exon','1','6774','.','+','.',  ###first 8 elements, tab delimited
            'gene_id "TPD36"',   ### remaining elements, semicolon delimited
            ' gene_version "1"',
            ' transcript_id "TPD36"',
            ' transcript_version "1"',
            ' gene_name "TPD36"',
            ' gene_source "FMI"',
            ' gene_biotype "protein_coding"',
            ' transcript_name "TPD36"',
            ' transcript_source "custom"',
            ' transcript_biotype "protein_coding"',
            ' protein_id "TPD36"',
            ' tag "basic";']
TPD36_part1 = '\t'.join((TPD36[0:8]))
TPD36_part2 = ';'.join((TPD36[8:20]))
TPD36_sense = '\t'.join((TPD36_part1,TPD36_part2))  ## there is one tab separation
TPD36_anti_sense = TPD36_sense.replace("+", "-")

FLG0010 = ['FLG0010','FLG0010','exon','1','6773','.','+','.',  ###first 8 elements, tab delimited
            'gene_id "FLG0010"',   ### remaining elements, semicolon delimited
            ' gene_version "1"',
            ' transcript_id "FLG0010"',
            ' transcript_version "1"',
            ' gene_name "FLG0010"',
            ' gene_source "FMI"',
            ' gene_biotype "protein_coding"',
            ' transcript_name "FLG0010"',
            ' transcript_source "custom"',
            ' transcript_biotype "protein_coding"',
            ' protein_id "TPD130"',
            ' tag "basic";']
FLG0010_part1 = '\t'.join((FLG0010[0:8]))
FLG0010_part2 = ';'.join((FLG0010[8:20]))
FLG0010_sense = '\t'.join((FLG0010_part1,FLG0010_part2))  ## there is one tab separation
FLG0010_anti_sense = FLG0010_sense.replace("+", "-")

FLG0006 = ['FLG0006','FLG0006','exon','1','6774','.','+','.',  ###first 8 elements, tab delimited
            'gene_id "FLG0006"',   ### remaining elements, semicolon delimited
            ' gene_version "1"',
            ' transcript_id "FLG0006"',
            ' transcript_version "1"',
            ' gene_name "FLG0006"',
            ' gene_source "FMI"',
            ' gene_biotype "protein_coding"',
            ' transcript_name "FLG0006"',
            ' transcript_source "custom"',
            ' transcript_biotype "protein_coding"',
            ' protein_id "FLG0006"',
            ' tag "basic";']
FLG0006_part1 = '\t'.join((FLG0006[0:8]))
FLG0006_part2 = ';'.join((FLG0006[8:20]))
FLG0006_sense = '\t'.join((FLG0006_part1,FLG0006_part2))  ## there is one tab separation
FLG0006_anti_sense = FLG0006_sense.replace("+", "-")

final_genes = '\n'.join((TPD36_sense,TPD36_anti_sense,FLG0010_sense,FLG0010_anti_sense,FLG0006_sense,FLG0006_anti_sense))

with open(mod_gtf_path, "a") as myfile:
    myfile.write(final_genes)
    myfile.close()
# remove features that have no gene_name

i = 0
with open(mod_gtf_path, "r") as file_input:
    with open(edited_gtf_path, "w") as output:
        for line in file_input:
            if i < 5:
                i = i + 1
                output.write(line)
            elif (i >= 5):
                i = i + 1
                result = re.search('gene_name(.*); gene_source', line)
                if result == None:
                    continue
                else:
                    output.write(line)

# filter the gtf, removes unneccesary elements that are not needed to reduce overlap

filter_command = "cellranger mkgtf {} {} --attribute=gene_biotype:protein_coding \
                   --attribute=gene_biotype:lincRNA \
                   --attribute=gene_biotype:antisense \
                   --attribute=gene_biotype:IG_LV_gene \
                   --attribute=gene_biotype:IG_V_gene \
                   --attribute=gene_biotype:IG_V_pseudogene \
                   --attribute=gene_biotype:IG_D_gene \
                   --attribute=gene_biotype:IG_J_gene \
                   --attribute=gene_biotype:IG_J_pseudogene \
                   --attribute=gene_biotype:IG_C_gene \
                   --attribute=gene_biotype:IG_C_pseudogene \
                   --attribute=gene_biotype:TR_V_gene \
                   --attribute=gene_biotype:TR_V_pseudogene \
                   --attribute=gene_biotype:TR_D_gene \
                   --attribute=gene_biotype:TR_J_gene \
                   --attribute=gene_biotype:TR_J_pseudogene \
                   --attribute=gene_biotype:TR_C_gene".format(edited_gtf_path,filtered_gtf_path)
os.system(filter_command)


# construct the bed file

bed_command = "gtf2bed < {} > {}".format(filtered_gtf_path,bed_path)
os.system(bed_command)

# construct the mapping genome
os.chdir('/home_fmi/01/otoosean/home_faimsrv01/Sequel/opt_items/artificial_promoters_genome_simple')
make_genome_command = "cellranger mkref --genome={} --fasta={} --genes={} --nthreads=30".format(genome_name,mod_fasta_path,filtered_gtf_path)
os.system(make_genome_command)
