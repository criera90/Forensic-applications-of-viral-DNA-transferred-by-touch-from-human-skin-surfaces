#!/util/opt/anaconda/deployed-conda-envs/packages/python/envs/python-3.9/bin/python

# functions
def usage():
    usage="""
        Usage:
        ./get_taxonomy_from_taxid.py -d <file.csv> -k <number_of_clusters> -e <epsilon> -i <number_of_iterations> -o <output_dir>
        Options:
           -d (directory)    path to a file with data. A tab delimited file
                             containing contigID, subject, p_value, pident, taxid
           -t (integer)      path to the directory containing nodes.dmp and 
                             rankedlineage.dmp filed from NCBI databases
           -o (directory)    output_dir with minimum an output file name
           -h                display this help and exit
    """
    print(usage)
def command_line_interaction():
    import sys, getopt
    # get arguments
    # stores in a tuple
    # https://www.cyberciti.biz/faq/python-command-line-arguments-argv-example/
    opts, args = getopt.getopt(sys.argv[1:], 'd:t:o:h')
    # initialize values
    dataset = ""
    n_clusters = -1
    epsilon = -1
    iterations = -1
    output_dir = ""
    for o,a in opts:
        if o == '-d':
            dataset = a     # a csv_file
        if o == '-t':
            tax_db_dir = a  # a path to a directory contating nodes and rankedlineage
        if o == '-o':
            output_dir = str(a)
        if o == '-h':
            usage()
            sys.exit(2)     # exits after -h is used
    return([dataset,tax_db_dir,output_dir])

def get_data(dataset):
    import pandas as pd
    import math
    data=[]
    # to add column names to the blast output
    # base it on the -outfmt "6 .." option on blastn
    col_names=["qaccver", "saccver", "sallseqid", "sscinames", "evalue", "pident", "qcovs", "staxids", "sskingdoms"]
    df= pd.read_csv(dataset,sep="\t", header=None, names=col_names)
    return df


# import libraries
import os, sys, getopt, time
import pandas as pd
import taxidTools

# getting input and parameters
inputs=command_line_interaction()
dataset=inputs[0]
tax_db_dir=inputs[1]
output_dir=inputs[2]
#print("Looking for %d clusters in dataset %s.\nUsing %d iterations with epsilon=%.2f" % (n_clusters, dataset, iterations, epsilon))

df=get_data(dataset)
contigs=df.qaccver.unique()
print(contigs)

# taxonomy ranks to consider form the taxonomy tree
rankFilter=['species', 'genus', 'subfamily', 'family', 'order', 'class', 'phylum', 'kingdom', 'clade', 'superkingdom']
# get taxonomy db into a taxidTools object
nodes=os.path.join(tax_db_dir,"nodes.dmp")
rankedlineage=os.path.join(tax_db_dir,"rankedlineage.dmp")
tax=taxidTools.Taxonomy.from_taxdump(nodes,rankedlineage)

entries=[]
for contig in contigs:
    df_temp=df.query('qaccver == @contig & evalue < 1e-17 & qcovs >= 30 & sskingdoms.notna()').sort_values(by=["qcovs"],ascending = False)
    if df_temp.empty == False:
        print(df_temp)
        contigID=df_temp["qaccver"].iloc[0]
        subjectACC=df_temp["saccver"].iloc[0]
        taxID=df_temp["staxids"].iloc[0]
        print(contigID)
        print(taxID)

        dic={}

        dic["contigID"]=contigID
        dic["subjectACC"]=subjectACC
        dic["taxID"]=taxID
        print(dic)

        # get taxonomy information
        l=tax.getAncestry(taxID)
        rankList=[name.rank for name in l]
        namesList=[name.name for name in l]
        print(rankList)
        print(namesList)

        n=0
        for rankL in rankList:
            if rankL in rankFilter:
                dic[rankL]=namesList[n]
            n+=1
        print(dic)
        entries.append(dic)

print(entries)

for ith_dict in entries:
    intersection_rankFilter_ith_entry=set.intersection(set(rankFilter),set(ith_dict))
    if set(rankFilter) == set(intersection_rankFilter_ith_entry):
        break
col_names=ith_dict.keys()
print(ith_dict.keys())

df=pd.DataFrame()

for col_name in col_names:
    temp_vector=[]
    for entry in entries:
        if col_name in entry:
            element=entry[col_name]
        else:
            element=str('unknown_'+col_name)
        temp_vector.append(element)
    temp_df=pd.DataFrame(temp_vector,columns=[col_name])
    df=pd.concat([df,temp_df],axis=1)
print(df)
df.to_csv(output_dir, index=False)


