#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 23:04:23 2022

@author: bengi
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 23:02:31 2022

@author: bengi
"""
#PATH: /Users/bengi/Desktop/CommsBioRevision_Final_9Dec/REVISION#2/DATA_CODES_new
import pandas as pd
import os

class StatisticalAnalysis_GenieTcga:
    def __init__(self, file_genei, file_maf):
        
#read GENIE and TCGA parsed mutation files
mut_file_genie = "GENIE_vol6_Mutation_Cleaned_Parsed_File__new.txt"

mut_file_tcga = "TCGA_MC3callSet_Mutation_Cleaned_Parsed_File__new.txt" 


SpecificColumns=['Hugo_Symbol','Tumor_Sample_Barcode',"WT_Residue",'HGVSp_Short', "VAF","Variant_Classification","WildTypeResidue","ResidueNumber",'MutantResidue']

df_genie= pd.read_csv(mut_file_genie,sep="\t",comment='#',usecols=SpecificColumns)

df_tcga= pd.read_csv(mut_file_tcga ,sep="\t",comment='#',usecols=SpecificColumns)

# drop rows which have same 'Hugo_Symbol','Tumor_Sample_Barcode', 'HGVSp_Short', "VAF",
#"Variant_Classification" and keep latest entry
new_tcga = df_tcga.drop_duplicates(
  subset = ['Hugo_Symbol','Tumor_Sample_Barcode', 'HGVSp_Short', "VAF","Variant_Classification"],
  keep = 'last').reset_index(drop = True)

new_genie = df_genie.drop_duplicates(
  subset = ['Hugo_Symbol','Tumor_Sample_Barcode', 'HGVSp_Short', "VAF","Variant_Classification"],
  keep = 'last').reset_index(drop = True)  

    
    

    
df_genie1=new_genie[new_genie["VAF"]>0.125]
df_tcga1=new_tcga[new_tcga["VAF"]>0.125] 

df_genie1.reset_index(inplace=True, drop=True)
df_tcga1.reset_index(inplace=True, drop=True)

    
#############   

# #check nonsense mutation
# df_tcga1.columns
# nonsense=df_tcga1[df_tcga1["Variant_Classification"]=="Nonsense_Mutation"]
# nonsense1=nonsense[nonsense["Hugo_Symbol"]=="ARID1A"]

#mutation occurs ~25% of the cells of sequenced tumor 

AllTumorCount=len(set(df_genie1["Tumor_Sample_Barcode"].to_list()).union(set(df_tcga1["Tumor_Sample_Barcode"].to_list()))) 
 #in both data 62,567 tumors with point mutations
# and  VAF>0.125
AllTumorCount
#
df_merged=pd.concat([df_genie1,df_tcga1], ignore_index=True)
#df_merged.to_csv("Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt",index=False,sep="\t",header=True)
########
df_merged["WildTypeResidue"]
Genes=set(df_genie1['Hugo_Symbol'].to_list()).union(set(df_tcga1['Hugo_Symbol'].to_list())  )                    
len(Genes) # # 19,415 genes with point mutations in both data sets

# Returns   #Gene:Mut:Sample nested dictionary for genie and tcga

def GENIE_TCGA_Gene_Mutation_Sample_Dictionary(df_genie1 = df_genie1 , df_tcga1 = df_tcga1):
    
    #Return {Gene:{Mut:Sample} }nested dictionary for genie and tcga
    
    Gene_Mutation_Sample_dict={gene:dict() for gene in Genes}
    
    for ind in df_genie1.index:
        gene=df_genie1.iloc[ind]['Hugo_Symbol']
        mutation=df_genie1.iloc[ind]['WT_Residue']
        tumor=df_genie1.iloc[ind]['Tumor_Sample_Barcode']
        if mutation not in Gene_Mutation_Sample_dict[gene].keys():
            Gene_Mutation_Sample_dict[gene][mutation]=set()
            Gene_Mutation_Sample_dict[gene][mutation].add(tumor)
        elif mutation in Gene_Mutation_Sample_dict[gene].keys():
            Gene_Mutation_Sample_dict[gene][mutation].add(tumor)
        else:
            print((gene,mutation))
            
    
    for ind in df_tcga1.index:
        gene=df_tcga1.iloc[ind]['Hugo_Symbol']
        mutation=df_tcga1.iloc[ind]['WT_Residue']
        tumor=df_tcga1.iloc[ind]['Tumor_Sample_Barcode']
        if mutation not in Gene_Mutation_Sample_dict[gene].keys():
            Gene_Mutation_Sample_dict[gene][mutation]=set()
            Gene_Mutation_Sample_dict[gene][mutation].add(tumor)
        elif mutation in Gene_Mutation_Sample_dict[gene].keys():
            Gene_Mutation_Sample_dict[gene][mutation].add(tumor)
        else:
            print((gene,mutation)) 
    return Gene_Mutation_Sample_dict
            
Gene_Mutation_Sample_dict = GENIE_TCGA_Gene_Mutation_Sample_Dictionary(df_genie1 = df_genie1 , df_tcga1 = df_tcga1)

len(Gene_Mutation_Sample_dict["PIK3CA"]["E726"])
len(Gene_Mutation_Sample_dict["PIK3CA"]["M1004"])
########
######
        #for each gene write mutations observed on >=3 tumors to file with "," separated tumor list and its size
Genes_Sorted=sorted(list(Genes)) #19415 genes
###write gene, mutation, tumors, and tumor_count to file
for gene in Genes_Sorted:
    for mutation in Gene_Mutation_Sample_dict[gene].keys():
        if len(Gene_Mutation_Sample_dict[gene][mutation])>=3: #use as paramete, initially we used 5 as a threshold
            AllTumor=",".join(list(Gene_Mutation_Sample_dict[gene][mutation]))
            TumorCount=str(len(Gene_Mutation_Sample_dict[gene][mutation]))
            try:
                with open("Gene_Mutation_TumorList__geq3_new.txt","a") as outfile:
                    outfile.write(str(gene)+"\t"+str(mutation)+"\t"+str(TumorCount)+"\t"+AllTumor+"\n")
            except TypeError:
                print((gene, mutation))

#####
###write gene, mutation, tumors, and tumor_count to file for the mutations observed on 2 tumors
for gene in Genes_Sorted:
    for mutation in Gene_Mutation_Sample_dict[gene].keys():
        if len(Gene_Mutation_Sample_dict[gene][mutation])<=2: #mutations observed on less than 2 tumors
            AllTumor=",".join(list(Gene_Mutation_Sample_dict[gene][mutation]))
            TumorCount=str(len(Gene_Mutation_Sample_dict[gene][mutation]))
            try:
                with open("Gene_Mutation_TumorList__LEQ_2_new.txt","a") as outfile:
                    outfile.write(str(gene)+"\t"+str(mutation)+"\t"+str(TumorCount)+"\t"+AllTumor+"\n")
            except TypeError:
                print((gene, mutation))

 #   
#Form the potential doublets from the mutations observed on >=3 tumors
#FolderName: /Users/bengi/Desktop/CommsBioRevision_Final_9Dec/REVISION#2/DATA_CODES_new
Gene_Mutation_List={} #each gene keep track of mutations with #tumors>=3 --12,724 key genes with 
#at least one mutation observed on more than 2 tumors
GeneMut_TumorDict={} #{Gene_Mutation:Tumor_set} dictionary, 65,872 mutations in total
with open("Gene_Mutation_TumorList__geq3_new.txt","r") as infile: 
    for line in infile:
        splitted=line.rstrip("\n").split("\t")
        gene=splitted[0]
        mutation=splitted[1]
        genemut=gene+"_"+mutation
        patient=set(splitted[-1].split(","))
        GeneMut_TumorDict[genemut]=patient
        if gene not in Gene_Mutation_List.keys():
            Gene_Mutation_List[gene]=[]
            Gene_Mutation_List[gene].append(mutation)
        elif gene in Gene_Mutation_List.keys():
            Gene_Mutation_List[gene].append(mutation)
len(GeneMut_TumorDict) 
len(Gene_Mutation_List)          
##Genes with 2 or more mutations that or observed on at least 3 tumors
#we need at least 2 mutations on a gene to form potential doublets
Genes_MutationFor_Combinations={}   #8189 genes with >=2 mutations
for gene in Gene_Mutation_List.keys():
    if len(Gene_Mutation_List[gene]) >=2:
        Genes_MutationFor_Combinations[gene]=Gene_Mutation_List[gene]
len(Genes_MutationFor_Combinations["PIK3CA"])   
len(Genes_MutationFor_Combinations)   
Genes_MutationFor_Combinations
#
####
#Form the potential doublets, be careful about (a,b)-(b,a) type duplicates
Doublets=[]
#2230203, doublets are formed from mutations observed on >=tumors
for gene in Genes_MutationFor_Combinations:
    for i in range(0,len(Genes_MutationFor_Combinations[gene])-1):
        for j in range(i+1,len(Genes_MutationFor_Combinations[gene])):
            Doublets.append((gene+"_"+Genes_MutationFor_Combinations[gene][i],gene+"_"+Genes_MutationFor_Combinations[gene][j]))
len(Doublets)
Doublets
#19,029 of the doublets observed on at least one tumor (intersection >0)
Doublets_Nonzero=[] 
for doublet in Doublets:
    if len(GeneMut_TumorDict[doublet[0]].intersection(GeneMut_TumorDict[doublet[1]]))>0:
        Doublets_Nonzero.append(doublet)
len(Doublets_Nonzero)
Doublets_Nonzero
##########
##########
#the doublets matching with our initial 228 doublets        
df_original_doublets=pd.read_excel("/Users/bengi/Desktop/CommsBioRevision_Final_9Dec/Supplementray_Table_1.xlsx",sheet_name="Statistics_")
OriginalDoublets=[]
df_original_doublets.columns
for inds in df_original_doublets.index:
    mut1=df_original_doublets.iloc[inds]['Gene']+"_"+df_original_doublets.iloc[inds]['Residue_1']
    mut2=df_original_doublets['Gene'].iloc[inds]+"_"+df_original_doublets.iloc[inds]['Residue_2']
    OriginalDoublets.append((mut1,mut2))
common=[] #146 of the original doublets exist now 
for double in Doublets_Nonzero:
    m1=double[0]
    m2=double[1]
    if ((m1,m2) in OriginalDoublets) or ((m2,m1) in OriginalDoublets):
        common.append(double)

#among 228 doublets which not exist now, eaxclude while uploading
###############
not_common=[]
for double in OriginalDoublets:
    m1=double[0]
    m2=double[1]
    if ((m1,m2) not in common) and ((m2,m1) not in common):
        not_common.append(double)
        
        #the doublets in the following file contains 82 doublets out of 228 initial significant doublets
        #now either the components' intersection is empty tumor set
        #or at least one component observed on less than 3 tumors, hence not taken into account
for double in not_common:
    m1=double[0]
    m2=double[1]
    try:
        count=len(GeneMut_TumorDict[m1].intersection(GeneMut_TumorDict[m2]))
        with open("In228butNotExistNow_geq3_new.txt","a") as outfile:
            outfile.write(m1+"\t"+m2+"\t"+str(count)+"\n")
    except KeyError:
        with open("In228butNotExistNow_geq3_new.txt","a") as outfile:
            outfile.write(m1+"\t"+m2+"\t"+"AtLeastOneComponentNotExist"+"\n")
        print(double)
 
#GeneMut_TumorDict["BRAF_I710"]
#GeneMut_TumorDict['TSC1_E878'].intersection(GeneMut_TumorDict['TSC1_R892'])
#which original doublets not exist
#this part could be excluded
######################
        
###
#Obtain the number of gene mutant tumors (mutations observed on >=3 tumors)
Gene_MutantTumors_Dict={}
with open("Gene_Mutation_TumorList__geq3_new.txt","r") as infile: 
    for line in infile:
        splitted=line.rstrip("\n").split("\t")
        gene=splitted[0]
        patient=set(splitted[-1].split(","))
        if gene not in Gene_MutantTumors_Dict.keys():
            Gene_MutantTumors_Dict[gene]=set()
            Gene_MutantTumors_Dict[gene]=Gene_MutantTumors_Dict[gene].union(patient)
        elif gene in Gene_MutantTumors_Dict.keys():
            Gene_MutantTumors_Dict[gene]=Gene_MutantTumors_Dict[gene].union(patient)
len(Gene_MutantTumors_Dict["PIK3CA"] )           
### 
###
#write the statistics for the potential doublets observed on >0 tumors without any further filtering
            #after this step read the file and compute p values for 4 different contingency tables
with open("SameGeneDoubletsTumorCountsNoFilter__geq3_new.txt","a") as outfile:
    outfile.write("Gene"+"\t"+"Mut1"+"\t"+"Mut2"+"\t"+"Among228Doubles"+"\t"+"DoubleMut#"+"\t"+"OnlyMut1"+"\t"+"OnlyMut2"+"\t"+"GeneMutant#"+"\n")

for doublet in Doublets_Nonzero:
    if doublet in common:
        tag="Yes"
    else:
        tag="No"
    gene=doublet[0].split("_")[0]
    mut1=doublet[0]
    mut2=doublet[1]
    doubleMutant=GeneMut_TumorDict[mut1].intersection(GeneMut_TumorDict[mut2])
    only1=GeneMut_TumorDict[mut1].difference(GeneMut_TumorDict[mut2])
    only2=GeneMut_TumorDict[mut2].difference(GeneMut_TumorDict[mut1])
    geneMutant=Gene_MutantTumors_Dict[gene]
    with open("SameGeneDoubletsTumorCountsNoFilter__geq3_new.txt","a") as outfile:
        outfile.write(gene+"\t"+mut1.split("_")[1]+"\t"+mut2.split("_")[1]+"\t"+tag+"\t"+str(len(doubleMutant))+"\t"+\
                     str(len(only1))+"\t"+str(len(only2))+"\t"+str(len(geneMutant))+"\n")
#SameGeneDoubletsTumorCountsNoFilter__.txt does not contains nan values on the columns mut1, mut2     
GeneMut_TumorDict.keys()        
import pandas as pd
from scipy.stats import fisher_exact
import numpy as np

df=pd.read_csv("SameGeneDoubletsTumorCountsNoFilter__geq3_new.txt",sep="\t")


#first write header
with open("SameGeneStatisticsAllDoublets_VariousContingencyTables__geq3_new.txt","a") as outfile:
    outfile.write("Gene"+"\t"+"Mut1"+"\t"+"Mut2"+"\t"+"Among228Doubles"+"\t"+"DoubleMut#"+"\t"+"OnlyMut1"+"\t"+"OnlyMut2"+"\t"+"GeneMutant#"+"\t"+\
                  "AllTumor#VAF>"+"\t"+"p1"+"\t"+"p2"+"\t"+"p3"+"\t"+"OddsR3"+"\t"+"p4"+"\t"+"OddsR4"+"\n")
#calculate p-values with fisher exact test for different contingency tables
for i in df.index:
    tag=df.iloc[i]['Among228Doubles']
    a=df.iloc[i]['DoubleMut#']
    b=df.iloc[i]["OnlyMut1"]
    c=df.iloc[i]["OnlyMut2"]
    d1=b+c
    d2=0
    d3=df.iloc[i]['GeneMutant#']-(a+b+c)
    d4=AllTumorCount-(a+b+c)
    table1=np.array([[a,b],[c,d1]]) #values are dependentt
    table2=np.array([[a,b],[c,d2]]) #values indep., evaluated only among A or B mutants
    table3=np.array([[a,b],[c,d3]]) #values indep., evaluated  among gene mutants
    table4=np.array([[a,b],[c,d4]]) #values indep., evaluated  among all tumors with vaf>0.125
    oddsr1, p1 = fisher_exact(table1, alternative='two-sided')
    oddsr2, p2 = fisher_exact(table2, alternative='two-sided')
    oddsr3, p3 = fisher_exact(table3, alternative='two-sided')
    oddsr4, p4 = fisher_exact(table4, alternative='two-sided')
    try:
        with open("SameGeneStatisticsAllDoublets_VariousContingencyTables__geq3_new.txt","a") as outfile:
            outfile.write(df.iloc[i]["Gene"]+"\t"+df.iloc[i]['Mut1']+"\t"+df.iloc[i]['Mut2']+"\t"+tag+"\t"+str(a)+"\t"+str(b)+"\t"+str(c)+"\t"+str(df.iloc[i]['GeneMutant#'])\
                          +"\t"+str(AllTumorCount)+"\t"+str(p1)+"\t"+str(p2) +"\t"+str(p3)+"\t"+ str(oddsr3)+"\t"+str(p4)+"\t"+ str(oddsr4) +"\n")
    except IndexError:
        print(i)
###double mutant tumor list        
for doublet in Doublets_Nonzero:
    if doublet in common:
        tag="Yes"
    else:
        tag="No"
    gene=doublet[0].split("_")[0]
    mut1=doublet[0]
    mut2=doublet[1]
    doubleMutant=GeneMut_TumorDict[mut1].intersection(GeneMut_TumorDict[mut2])
    Tumors=",".join(doubleMutant)     
    with open("AllDoubletsObservedAtLeastOneTumorDoubleMutantTumorsInfo_new.txt","a") as outfile:
        outfile.write(gene+"\t"+mut1+"\t"+mut2+"\t"+tag+"\t"+str(len(doubleMutant))+"\t"+Tumors+"\n")
        
