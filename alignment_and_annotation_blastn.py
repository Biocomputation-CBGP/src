#Program to  align and annotate the reads that come from sanger sequencing to a reference genome with available annotations
#This script has been tested with bacterial genomes

#Import the neccesary python packages
import sys
import pandas as pd
import os
import subprocess



#Control of arguments (Reminder of arguments and order for input, all are directories)
if len(sys.argv) != 5:
    print("""---------------------------------------------------------------------------------------------------
Need 4 arguments:
    Directory of sequencing results (reads)
    Reads file type (lower case) {txt|seq}
    Whole genome file (fasta format ex faa or fna)
    Genome Annotation file (csv format)
---------------------------------------------------------------------------------------------------""")
    quit()
    
else:
    directory_files = sys.argv[1]
    type_files = sys.argv[2] #txt and seq for now
    database = sys.argv[3]
    file_annotation = sys.argv[4]
    pass

#Checking if the output directory is there so we do not overwrite a directory with the same name
#Output folder is name results_script_blast
if os.path.isdir('./results_script_blast'):
    txt = input("results_script_blast is going to be overwrited, are you sure? [y/n]:  ")
    if txt == "y":
        os.system("mkdir -p results_script_blast")
    elif txt == "n":
        quit() #We exit from the program in general
    else:
        print("""------------------------------
Assumed NO, Exiting program
------------------------------""")
        quit() #We exit from the program in general
else:
    os.system("mkdir results_script_blast") #Directory in which we are going to store the results


#-------------------------------------------
#Creation of compressed sequencing reads file to align to reference genome with blastn
#Format other than .seq or .txt can need other processing, to be included in this section

if type_files == 'seq':
    os.system("find "+directory_files+"/. -regex '.*[^_raw]\.seq' -exec cat {} \; -exec echo \; > all_reads_merged.fna")
    #This expression takes all the .seq files of the given directory, it prints them, makes a new line and put all this into all_reads_merged.fna
elif type_files == 'txt':
    os.system("find "+directory_files+"/. -regex '.*\.txt' -exec cat {} \; -exec echo \; > all_reads_merged_tabs.txt")
    os.system("sed -e 's/\t/ /g' all_reads_merged_tabs.txt > all_reads_merged.fna")
    os.system("rm all_reads_merged_tabs.txt")
    #These expresions takes the .txt files and later it removes the \t (tab) from the first lines (description lines in fasta)
else:
    print("Only .txt and .seq files supported")
    quit() #quits the program



try: #We check if the database already exists, it works not only with the name of the db but also with the path to it
    subprocess.check_call(["blastdbcmd", "-info", "-db", database], stdout = subprocess.DEVNULL, stderr = subprocess.DEVNULL)
    #check_call makes the expression in the list and checks if there is a exit code 0 (success) or other (error)
    #stdout and stderr are th eoutput and error message to a file that is temporal and not stored
    #If it is a success it means that the indexed genome works
    print("""------------------------
------------------------
Database already exists!
Info of the DB:""")
    os.system("blastdbcmd -info -db "+database)
    print("""------------------------
------------------------""")
    #gives the info for the database
except:
    #Create the index so blastn can do the alignment
    os.system('makeblastdb -in '+database+' -dbtype nucl')
    #in case that check_all gives and error it creates the dabatase from the file given

#-------------------------------------------
#-------------------------------------------
#Selected blastn output columns and genome annotations columns and significant score hit to elaborate final annotated table
## IMPORTANT!  VARIABLE PART OF THE SCRIPT (MODULAR)

columns_output_blastn = ["std", "sstrand"]#This will be filled in case that we do not want the standard output
                          #If empty will be the std output (check blastn manual for columns in standard output)
                          #It is needed that the query acc. and the bit score is there. The others are optional

#Determine the columns that we want
#Columns in the sequenece alignment file that we want (needs to be typed exactly as in the output of blastn)
columns_seq_alig = ["query acc.", "s. start", "% identity", "alignment length", "mismatches", "gap opens", "evalue", "bit score", "subject strand"] #According to the variable columns_output_blastn, same order as blastn output
#It will give the output in the order of naming but without repetitions (if we put std score it will only print std)
#IMPORTANT: query acc. and bit score is always needed
#Columns of the gene annotation file that we want in the final annotated table
columns_ann = ["Locus Tag","Feature Type","Start","End","Strand","Gene Name","Product Name","Subcellular Localization [Confidence Class]"]
#IMPORTANT:Locus Tag column are always needed

range_value = 0.01 #This is a variable we can set to identify significant alignments, score has been assesed according to our experimental results
#This range value is the proportion that the score can be lower than the best hit in the allignment that we allow as a "valid" allignment
#This variable affects which reads have multiple allignments or not and the locus that we take as a hit


#-------------------------------------------
#------------------------------------------- REST OF SCRIPT IS NOT MODULAR (MAKE CHANGES AT YOUR OWN RISK)


#We set the headers for the blastn output file (with -outfmt 6 there are no headers) with columns_output_blastn variable
#Extra columns added are formatted adequately to proceed with the script

if len(columns_output_blastn) > 0:
    command_of_6_17 = " ".join(columns_output_blastn)+"'"
    print(command_of_6_17)
    header_output_6 = "query acc.\tsubject acc.\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score"+"\tsubject strand"
    #header_output_6 = "\t".join(columns_output_blastn)
else:
    command_of_6_17 = "std'"
    header_output_6 =  "query acc.\tsubject acc.\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score"

#Perform Blastn
command_blastn_tabular_output = "blastn -query all_reads_merged.fna -db "+database+" -out ./results_script_blast/all_seq_aligned.txt -outfmt '6 "+command_of_6_17
print(command_blastn_tabular_output)
command_blastn_SAM_output = "blastn -query all_reads_merged.fna -db "+database+" -out ./results_script_blast/all_seq_aligned.sam -outfmt '17 "+command_of_6_17
#It is possible to adjust with other arguments these expressions (check blastn manual)
#The way of changing the output is -outfmt "6 std" for the standard output, other example  -outfmt "6 std staxid" this will give us the standard output and the tax id of the subject
#It will give the output in the order of naming but without repetitions (if we put std score it will only print std)

os.system(command_blastn_SAM_output)
os.system(command_blastn_tabular_output)
os.system("rm all_reads_merged.fna") #We remove the file all_reads_merged but we can keep it deleting this command (or commenting)

# ----------------------------------

#Formatting of blastn output file 
command_add_header_txt = "echo '"+header_output_6+"' | cat - ./results_script_blast/all_seq_aligned.txt > temp && mv temp ./results_script_blast/all_seq_aligned.txt"
os.system(command_add_header_txt)
#add the header expression to the output of blastn

#state the files of blast, genome annotation file and the final table
file_magicblast = "./results_script_blast/all_seq_aligned.txt"
final_table_name = "./results_script_blast/table_reads_genes_description.csv"

#Load tables as pandas dataframes
table_ann = pd.read_csv(file_annotation)

table_seq = pd.read_table(file_magicblast) #Tochange to blastn

#Create a cross join table from sequences alignment and annotation to select rows with matches

#We create a column called key to be able to do a cross table (all possible combinations) with both tables
table_seq['key'] = 1
table_ann['key'] = 1
table_cross = pd.merge(table_seq, table_ann, on ='key').drop("key", 1)
#after doing this cross table we delete the column key


#Filter rows by only taking the rows that the subject (reference genome in this case) starting position is between the annotated genes
#this filtering is needed because there are lines in the table that the subject positions does not correspond to the annotated genes positions
table_matches = table_cross[(table_cross["s. start"] >= table_cross["Start"]) & (table_cross["End"] >= table_cross["s. start"])]

#filter columns saving only the ones that have been selected previously in columns_seq_alig and columns_ann 
table_matches = table_matches[columns_seq_alig + columns_ann]


#Create a table with allignments that have not been matched with any annoted genes
table_not_matches = table_seq[~table_seq["query acc."].isin(table_matches["query acc."])][columns_seq_alig]

#Create the final table that is the alignments matched with genes with the ones not matched with any annotated gene
final_table = pd.concat([table_matches, table_not_matches], sort = False)

#------------------------------------------------------------

#We remove the allignments that are not within the threshold that we put so we discard these alignments according to the variable range_value
#We are adding a column with the score of the best hit of that read and then we delete the ones that are lower than score*threshold
final_table["Highest bit score"] = final_table.groupby('query acc.', sort=False)['bit score'].transform('max')
final_table.drop(final_table[final_table["bit score"] < final_table["Highest bit score"]*(1-range_value)].index, inplace=True)
#now we can erase this column because it has already served its purpose
del final_table["Highest bit score"]

#Fill the Na (non available) in pandas with the character - for stetic purposes
final_table.fillna("-",inplace=True)

#We add the "warning" column in which we say if there are duplicates or not
final_table["Multiple Allignments"] = final_table.duplicated(subset = ["query acc."], keep = False)

#Now we create the multiple locus column in case there are multiple allignments
locus_associated = [] 
for group_q_acc in final_table.groupby('query acc.', sort = False):#we group by query acc.
    locus_associated_query = list(final_table[final_table["query acc."] == group_q_acc[0]]["Locus Tag"].values) #list of the locuses attached to the query read
    if len(locus_associated_query) > 1:
        locus_associated_query = locus_associated_query[1:]
    else:
        locus_associated_query = "-"
    locus_associated.append(locus_associated_query)

#Now we drop the duplicates only keeping the best alignment
#Warning: duplicates will also be dropped for alignments with the same score or withing the threshold
final_table.drop_duplicates(subset ="query acc.", inplace=True, keep = "first")
#This is the moment in which we only keep the best hit for the allignment

#We add to the table the locus column
final_table["Rest of Locus Tag Associated"] = locus_associated

#------------------

#We export the table that we have created as final_table_name.csv
final_table.to_csv(final_table_name, index = False)

#------------------
#Final encouraing message
print("Done :)")
#------------------