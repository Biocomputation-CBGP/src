###ESTE SCRIPT NO ESTA SUBIDO A GITHUB Y ES EN EL QUE USAMOS BLAST EN VEZ DE MAGICBLAST PERO SEGUIMOS QUITANDO DUPLICADOS
###ESTO PODRIA SER UN PROBLEMA AHORA, ESTAMOS ESPERANDO LA DECISION DE LOREA PARA SABER QUE DARLE AL FINAL


#Program for  align and annotate the reads that come from a sequentation

#Import the neccesary python packages
import sys
import pandas as pd
import os



# Control of arguments (So we do not give less or more arguments than needed)
if len(sys.argv) != 5:
    print("""---------------------------------------------------------------------------------------------------
Need 4 arguments:
    Directory of reads
    Reads file type (lower case)
    Database to do MagicBLAST (fasta format)
    Annotation file (csv format)
---------------------------------------------------------------------------------------------------""")
    quit()
    
else:
    directory_files = sys.argv[1]
    type_files = sys.argv[2] #txt and seq for now
    database = sys.argv[3]
    file_annotation = sys.argv[4]
    pass

#Checking if the directory is there so we do not overwritted
if os.path.isdir('./results_script_magicblast'):
    txt = input("results_script_magicblast is going to be overwrited, are you sure? [y/n]:  ")
    if txt == "y":
        os.system("mkdir -p results_script_magicblast")
    elif txt == "n":
        quit()
    else:
        print("""------------------------------
Assumed NO, Exiting program
------------------------------""")
        quit()
else:
    os.system("mkdir results_script_magicblast") #Directory in which we are going to store the results


#-------------------------------------------
#Creation of compressed reads file to align with magicblast

if type_files == 'seq':
    os.system("find "+directory_files+"/. -regex '.*[^_raw]\.seq' -exec cat {} \; -exec echo \; > all_reads_merged.txt")
    
    
elif type_files == 'txt':
    os.system("find "+directory_files+"/. -regex '.*\.txt' -exec cat {} \; -exec echo \; > all_reads_merged_tabs.txt")
    os.system("sed -e 's/\t/ /g' all_reads_merged_tabs.txt > all_reads_merged.txt")
    os.system("rm all_reads_merged_tabs.txt")

    
else:
    print("Only .txt and .seq files supported")
    quit()



#-------------------------------------------
#-------------------------------------------

#Create the index so magicblast can do the alignment
os.system('makeblastdb -in '+database+' -dbtype nucl')
#Magicblast
command_magicblast_2 = "blastn -query all_reads_merged.txt -db "+database+" -out ./results_script_magicblast/all_seq_aligned.txt -outfmt 6"
command_magicblast = "blastn -query all_reads_merged.txt -db Putida -out ./results_script_magicblast/all_seq_aligned.sam -outfmt 17"

#Si quisieramos poner mas columnas seria como lo siguiente '-outfmt "6 std qlen"'
#https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

os.system(command_magicblast)
os.system(command_magicblast_2)
# = "sed -i '1s/^/query\sacc.\tsubject\sacc.\t%\sidentity\talignment\slength\tmismatches\tgap\sopens\tq.\sstart\tq.\send\ts.\sstart\ts.\send\tevalue\tbit\sscore\n/' -f results_script_magicblast\/all_seq_aligned.txt"
command_add_header_txt = "echo 'query acc.\tsubject acc.\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score' | cat - ./results_script_magicblast/all_seq_aligned.txt > temp && mv temp ./results_script_magicblast/all_seq_aligned.txt"

os.system(command_add_header_txt)
#sed: -e expression #1, char 126: unterminated `s' command
#os.system("rm all_reads_merged.txt")

#-------------------------------------------
#-------------------------------------------




#------------------#ONLY CHANGABLE PART OF SCRIPT!!!!!
#Determine the columns that we want
#Columns in the sequenece alignment file that we want
columns_seq_alig = ["query acc.", "s. start", "% identity", "alignment length", "mismatches", "gap opens", "evalue", "bit score"]
#Estos no estan en el de blast "query strand","reference strand"
#IMPORTANT: query acc. is always needed
#Columns of the gene annotation file that we want
columns_ann = ["Locus Tag","Feature Type","Start","End","Strand","Gene Name","Product Name","Subcellular Localization [Confidence Class]"]
#------------------

#------------------
file_magicblast = "./results_script_magicblast/all_seq_aligned.txt"
final_table_name = "./results_script_magicblast/table_reads_genes_description.csv"

#Load tables and fix them to be pandas tables
table_ann = pd.read_csv(file_annotation)
# print(table_ann)

table_seq = pd.read_table(file_magicblast)

# headers_seq = file_seq[0].split("\t") #Text processing needed to create the pands table
# #file_seq = file_seq[3:]
# for i, read in enumerate(file_seq):
#     file_seq[i] = read.strip().split("\t")
# table_seq = pd.DataFrame(file_seq, columns=headers_seq)
#Remove duplicates in case of existing
table_seq.drop_duplicates(subset ="query acc.", inplace=True)
#Change type of some columns from str to numeric
table_seq[["s. start"]] = table_seq[["s. start"]].apply(pd.to_numeric) #To perform numeric comparation


#------------------

#------------------
#Create a cross join table to select rows with matches

table_seq['key'] = 1
table_ann['key'] = 1
table_cross = pd.merge(table_seq, table_ann, on ='key').drop("key", 1)


#Filter rows
table_matches = table_cross[(table_cross["End"] >= table_cross["s. start"]) & (table_cross["s. start"] >= table_cross["Start"])]
table_matches = table_matches[columns_seq_alig + columns_ann]

table_matches['Start'] = table_matches['Start'].astype(str).apply(lambda x: x.replace('.0',''))
table_matches['End'] = table_matches['End'].astype(str).apply(lambda x: x.replace('.0',''))

#Create a table with reads that have not been matched with annoted genes
table_not_matches = table_seq[~table_seq["query acc."].isin(table_matches["query acc."])][columns_seq_alig]

#Create the final csv file
final_table = pd.concat([table_matches, table_not_matches],sort = False)
final_table.to_csv(final_table_name,index=False)
#------------------

#------------------
#Final message
print("Done :)")
#------------------