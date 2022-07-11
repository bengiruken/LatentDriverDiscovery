#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 11:18:49 2022

@author: bengi
"""
###1..TCGA to PICKLE
from Preprocessing_TCGA import TCGA_MAF

from Preprocessing_TCGA import Convert_Maf_to_Txt

#Convert_Maf_to_Txt()

mut_file_tcga = "MC3_TCGA_Mutations.txt" 
path = ".Data"       
maf_tcga = TCGA_MAF(path,mut_file_tcga )  

 


DF_tcga=maf_tcga.PreprocessingTcgaMutationFile()

import pickle
#pickle.dump(DF_tcga, open("./Data/Preprocessed_TCGA_df.p", "wb"))    
#############3

#2.GENIE to PICKLE

from Preprocessing_GENIE import GENIE_MAF

path_genie="./Data"

"""data_mutations_extended.txt (AACR GENIE vol. 6.1) can be downloaded from https://doi.org/10.7303/syn20333031 into the Data folder"""
mut_file_genie = "data_mutations_extended.txt"
maf_genie = GENIE_MAF(path_genie,mut_file_genie )  

df_genie=maf_genie.PreprocessingGenieMutationFile()

df_genie_tumor=maf_genie.return_GENIE_patient_tumor_df()

# import pickle
# pickle.dump(df_genie, open("./Data/Preprocessed_GENIE_df.p", "wb")) 

# pickle.dump(df_genie_tumor, open("./Data/Preprocessed_GENIE_Patient_Sample_Info_df.p", "wb"))


##########
#3_Get necessary files dictionaries from preprocessed pickle files of tcga
#write dictionaries as pickle
#
path_merged = "./Data"
tcga_file_name  = "Preprocessed_TCGA_df.p"
genie_file_name = "Preprocessed_GENIE_df.p"
import pickle

from MERGED_GENIE_TCGA_FILES_preprocessed import TCGA_GENIE_MERGED_DATA

merged_obj = TCGA_GENIE_MERGED_DATA(path_merged, tcga_file_name, genie_file_name, Vaf_Value=0.125)
merged_df = merged_obj.Merge_TCGA_GENIE_VAF_filtered()

all_tumors = merged_obj.AllTumors()
Mutation_List_GEQ3 = merged_obj.MutationsList_GEQ3()

###########3
#tumor same gene double mutations
"""
def Tumor_Same_Gene_Doublets():
    
    Tumor_Double_Dict = {}
    import pickle
    Tumor_GeneMutation_Dictionary=pickle.load(open("Merged_TCGA_GENIE_Tumor_GeneMutation_GEQ3_Dictionary.p", "rb"))    #merged_obj.Merged_TCGA_GENIE_Tumor_GeneMutation_GEQ3_Dictionary()
    for tumor in Tumor_GeneMutation_Dictionary.keys():
        if len(Tumor_GeneMutation_Dictionary[tumor])>1:
            if tumor not in Tumor_Double_Dict.keys():
                Tumor_Double_Dict[tumor]=set()
                
                for mut1 in Tumor_GeneMutation_Dictionary[tumor]:
                    for mut2 in Tumor_GeneMutation_Dictionary[tumor]:
                        if (mut1!=mut2) and (mut1.split("_")[0]==mut2.split("_")[0]):
                            if ((mut1,mut2) not in Tumor_Double_Dict[tumor]) and ((mut2,mut1) not in Tumor_Double_Dict[tumor]):
                                Tumor_Double_Dict[tumor].add((mut1,mut2))
            elif tumor in Tumor_Double_Dict.keys():
                print(tumor)
                for mut1 in Tumor_GeneMutation_Dictionary[tumor]:
                    for mut2 in Tumor_GeneMutation_Dictionary[tumor]:
                        if (mut1!=mut2) and (mut1.split("_")[0]==mut2.split("_")[0]):
                            if ((mut1,mut2) not in Tumor_Double_Dict[tumor]) and ((mut2,mut1) not in Tumor_Double_Dict[tumor]):
                                Tumor_Double_Dict[tumor].add((mut1,mut2))
                    
    return Tumor_Double_Dict

Tumor_Same_Gene_Double_Dictionary = Tumor_Same_Gene_Doublets()

pickle.dump(Tumor_Same_Gene_Double_Dictionary, open("Tumor_Potential_Same_Gene_Doubles_Dictionary.p", "wb")) 
"""
Tumor_Same_Gene_Double_Dictionary=pickle.load(open("Tumor_Potential_Same_Gene_Doubles_Dictionary.p", "rb"))    #merged_obj.Merged_TCGA_GENIE_Tumor_GeneMutation_GEQ3_Dictionary()
############3
no_double_tumor = 0
All_Doubles = set()
for tumor in Tumor_Same_Gene_Double_Dictionary.keys():
    if len(Tumor_Same_Gene_Double_Dictionary[tumor])>0:
        for double in Tumor_Same_Gene_Double_Dictionary[tumor]:
            if (double not in All_Doubles) and  ((double[1],double[0]) not in All_Doubles):
                All_Doubles.add(double)
    else:
        no_double_tumor+=1
##output
print("There are {} tumors with at least two mutations observed on at least three tumors".format(len(Tumor_Same_Gene_Double_Dictionary)) )        
print("Among these, {} tumors do not have any same-gene double mutation".format(no_double_tumor ))  

len(Tumor_Same_Gene_Double_Dictionary)-no_double_tumor   


#pickle.dump(All_Doubles, open("All_Potential_Same_Gene_Doubles_Set.p", "wb")) 
##################
#All gene mutants (mutations >=3 tumors are included only)
Gene_Mutation_Tumor_Dictionary=pickle.load(open("./Data/Merged_GENIE_TCGA_Gene_Mutation_Tumor_Dict.p", "rb"))    #merged_obj.Merged_TCGA_GENIE_Tumor_GeneMutation_GEQ3_Dictionary()
Mutation_List_GEQ3 

All_Genes = sorted(list(Gene_Mutation_Tumor_Dictionary.keys()))
"""
Gene_All_Mutants_Dict = {gene:set() for gene in All_Genes} 
#get mutations in Mutation_List_GEQ3  observed on more than two tumors
for gene in Gene_All_Mutants_Dict.keys():
    for mut in Gene_Mutation_Tumor_Dictionary[gene].keys():
        gene_mut = gene+"_"+mut
        if gene_mut in Mutation_List_GEQ3:
            Gene_All_Mutants_Dict[gene] = Gene_All_Mutants_Dict[gene].union(Gene_Mutation_Tumor_Dictionary[gene][mut])
            
len(Gene_All_Mutants_Dict["PIK3CA"])            
pickle.dump(Gene_All_Mutants_Dict, open("./Data/Gene_All_Mutants_Dictionary.p", "wb")) 
"""  
Gene_All_Mutants_Dictionary =  pickle.load(open("./Data/Gene_All_Mutants_Dictionary.p", "rb")) 
len(Gene_All_Mutants_Dictionary["PIK3CA"])         
#################3
#FISHER EXACT TEST AMONG ALL TUMORS with at least one mutation of VAF>0.125
from scipy.stats import fisher_exact
import numpy as np
#first write header
with open("./Data/Potential_SameGeneDoubles_Fisher_Exact_Test.txt","a") as outfile:
    outfile.write("Gene"+"\t"+"Mut1"+"\t"+"Mut2"+"\t"+"DoubleMut#"+"\t"+"OnlyMut1"+"\t"+"OnlyMut2"+"\t"+"GeneMutant#"+"\t"+\
                  "AllTumor#VAF>"+"\t"+"p4"+"\t"+"OddsR4"+"\n")
#calculate p-values with fisher exact test for different contingency tables
for double in All_Doubles:
    double1 = double[0]
    double2 = double[1]
    gene = double1.split("_")[0]
    mut1 = double1.split("_")[1]
    mut2 = double2.split("_")[1]
    gene_all_mutants = len(Gene_All_Mutants_Dictionary[gene])
    a = len(Gene_Mutation_Tumor_Dictionary[gene][mut1].intersection(Gene_Mutation_Tumor_Dictionary[gene][mut2]))
    b = len(Gene_Mutation_Tumor_Dictionary[gene][mut1].difference(Gene_Mutation_Tumor_Dictionary[gene][mut2]))
    c = len(Gene_Mutation_Tumor_Dictionary[gene][mut2].difference(Gene_Mutation_Tumor_Dictionary[gene][mut1]))
    
    d4=len(all_tumors)-(a+b+c)
 
    table4=np.array([[a,b],[c,d4]]) #values indep., evaluated  among all tumors with vaf>0.125

    oddsr4, p4 = fisher_exact(table4, alternative='two-sided')

    with open("./Data/Potential_SameGeneDoubles_Fisher_Exact_Test.txt","a") as outfile:
        outfile.write(gene+"\t"+mut1+"\t"+mut2+"\t"+str(a)+"\t"+str(b)+"\t"+str(c)+"\t"+str(gene_all_mutants)+"\t"+str(len(all_tumors))+"\t"+str(p4)+"\t"+ str(oddsr4) +"\n")

#################3
#MULTIPLE TESTING-FDR CORRECTION
####Label OG-TSG, KnownDriver[D], driver[d] from the files


##################
########
#FDR corrected p val with benjamini hochberg method
#!pip install statsmodels
def FDR_Correction():
    import statsmodels.stats.multitest
    import pandas as pd
        
    from pandas import ExcelWriter
    import numpy as np
    
    OG_TSG_Dict =  pickle.load(open("./Data/OG_TSG_Dictionary.p", "rb")) 
    OncogenicMutations =  pickle.load(open("OncogenicMutations_Dictionary.p", "rb")) 
    
    df_stats=pd.read_csv("Potential_SameGeneDoubles_Fisher_Exact_Test.txt",sep="\t")
    df_stats.columns
    """
    ###write nonsignificant doublets to file
    df_stats_nonsig=df_stats[df_stats['p4']>=0.05]
    
    from pandas import ExcelWriter
    from pandas import ExcelFile
    import numpy as np
    
    
    writer = ExcelWriter('P4_SameGeneDoublet_NonSignificant__geq3_new.xlsx', engine='xlsxwriter')
    df_stats_nonsig.to_excel(writer,'P4NonSignificant',index=False,header=True)
    # Close the Pandas Excel writer and output the Excel file.
    writer.save()
    """
    
    ##significant doublets and multiple testing
    #P2<0.05 
    
    df_stats["OG_TSG"]=["None"]*len(df_stats)
    
    for ind in df_stats.index:
        gene=df_stats.iloc[ind]["Gene"]
        if gene in OG_TSG_Dict.keys():
            df_stats.at[ind,"OG_TSG"]=OG_TSG_Dict[gene]
        else:
            df_stats.at[ind,"OG_TSG"]="NotSpecified"
    ###check oncogenic drivers
    df_stats["OncogenicDriver_Mut1"]=["None"]*len(df_stats)
    df_stats["OncogenicDriver_Mut2"]=["None"]*len(df_stats)
    
    for ind in df_stats.index:
        gene=df_stats.iloc[ind]["Gene"]
        mut1=df_stats.iloc[ind]["Mut1"]
        genemut=gene+"_"+mut1
        if genemut in OncogenicMutations:
            df_stats.at[ind,"OncogenicDriver_Mut1"]="KnownDriver[D]"
        else:
            df_stats.at[ind,"OncogenicDriver_Mut1"]="driver[d]"
    
    for ind in df_stats.index:
        gene=df_stats.iloc[ind]["Gene"]
        mut2=df_stats.iloc[ind]["Mut2"]
        genemut2=gene+"_"+mut2
        if genemut2 in OncogenicMutations:
            df_stats.at[ind,"OncogenicDriver_Mut2"]="KnownDriver[D]"
        else:
            df_stats.at[ind,"OncogenicDriver_Mut2"]="driver[d]"
    
    ####
    #Multiple Testing Benjamini hochberg
    P4=df_stats["p4"].to_list()
    q4=statsmodels.stats.multitest.multipletests(P4, alpha=0.1, method='fdr_bh', is_sorted=False, returnsorted=False)
    reject4=q4[0]    
    qval4=q4[1]
    
    df_stats["Reject4"]=reject4
    df_stats["Qval4"]=qval4
    df_stats["Log2_OR4"]=np.log2(df_stats["OddsR4"])
    
    df_stats["Tendency"]=["None"]*len(df_stats)
    
    for ind in df_stats.index:
        if df_stats.iloc[ind]["Log2_OR4"]>0:
            df_stats.at[ind,"Tendency"]="Co-occurence"
        else:
            df_stats.at[ind,"Tendency"]="Mutual Exclusivity"
        
    
    df_stats.to_csv("./Data/Potential_SameGeneDoubles_Fisher_Exact_Test_FDR_Multiple_Testing.txt",sep="\t",index=False,header=True)
    

    ####
    writer = ExcelWriter('./Data/Potential_SameGeneDoubles_Fisher_Exact_Test_FDR_Multiple_Testing.xlsx', engine='xlsxwriter')
    df_stats.to_excel(writer,'Statistics_AmongAllTumors',index=False,header=True)
    # Close the Pandas Excel writer and output the Excel file.
    writer.save()

FDR_Correction()
#############3


####3 LATENT DRIVER LABELING



import pandas as pd

df_significant = pd.read_csv("./Data/Potential_SameGeneDoubles_Fisher_Exact_Test_FDR_Multiple_Testing.txt", sep="\t")
df_significant.columns

df_significant = df_significant[(df_significant['p4']<0.05)&(df_significant['Qval4']<0.1)]
df_significant.reset_index(drop=True,inplace=True)

df_significant["Mut1_All"] = df_significant['OnlyMut1'] + df_significant["DoubleMut#"]
df_significant["Mut2_All"] = df_significant['OnlyMut2'] + df_significant["DoubleMut#"]

df_significant["Fraction1"] = (df_significant["Mut1_All"] /  df_significant['GeneMutant#'])*100
df_significant["Fraction2"] = (df_significant["Mut2_All"] /  df_significant['GeneMutant#'])*100

df_significant.head()
df_significant["Class1"] = [None]*len(df_significant)
df_significant["Class2"] = [None]*len(df_significant)
####
#Class 1&2: strong driver if known driver with fraction among gene mutants >10%
#weak driver if known driver with fraction among gene mutants <=10%
#strong latent driver if driver with fraction among gene mutants >1%
#weak latent driver if driver with fraction among gene mutants <=1%

for ind in df_significant.index:
    if df_significant.iloc[ind]["OncogenicDriver_Mut1"]  == "KnownDriver[D]" : #"driver[d]"
        if df_significant.iloc[ind]["Fraction1"] > 10:
            df_significant.at[ind,"Class1"] = "StrongDriver"
        elif df_significant.iloc[ind]["Fraction1"] <= 10:
            df_significant.at[ind,"Class1"] = "WeakDriver"
    elif df_significant.iloc[ind]["OncogenicDriver_Mut1"]  == "driver[d]" : #"driver[d]"
        if df_significant.iloc[ind]["Fraction1"] > 1:
            df_significant.at[ind,"Class1"] = "StrongLatentDriver"
        elif df_significant.iloc[ind]["Fraction1"] <= 1:
            df_significant.at[ind,"Class1"] = "WeakLatentDriver"
for ind in df_significant.index:
    if df_significant.iloc[ind]["OncogenicDriver_Mut2"]  == "KnownDriver[D]" : #"driver[d]"
        if df_significant.iloc[ind]["Fraction2"] > 10:
            df_significant.at[ind,"Class2"] = "StrongDriver"
        elif df_significant.iloc[ind]["Fraction2"] <= 10:
            df_significant.at[ind,"Class2"] = "WeakDriver"
    elif df_significant.iloc[ind]["OncogenicDriver_Mut2"]  == "driver[d]" : #"driver[d]"
        if df_significant.iloc[ind]["Fraction2"] > 1:
            df_significant.at[ind,"Class2"] = "StrongLatentDriver"
        elif df_significant.iloc[ind]["Fraction2"] <= 1:
            df_significant.at[ind,"Class2"] = "WeakLatentDriver"



############### 
#Specific Ordering wrt Activation
df_new1 = df_significant
df_new1['Activation'] = 0 # add a class column with 0 as default value

# driver+wd veya d+d olanlar k
df_new1.loc[(df_new1['Class1']=="StrongDriver") & #  
       (df_new1['Class2']=="StrongDriver"), 'Activation'] = 120 # then set class to 120

df_new1.loc[((df_new1['Class1']=="StrongDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="WeakDriver"))|((df_new1['Class1']=="WeakDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="StrongDriver")), 'Activation'] = 100 # then set class to 1

# find all rows that fulfills your conditions and set class to 1
df_new1.loc[((df_new1['Class1']=="StrongDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="StrongLatentDriver"))|((df_new1['Class1']=="StrongLatentDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="StrongDriver")), 'Activation'] = 100 # then set class to 1
# find all rows that fulfills your conditions and set class to 1
df_new1.loc[((df_new1['Class1']=="StrongDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="WeakLatentDriver"))|((df_new1['Class1']=="WeakLatentDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="StrongDriver")), 'Activation'] = 90 # then set class to 1

# driver+wd veya d+d olanlar k
df_new1.loc[((df_new1['Class1']=="WeakDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="WeakDriver")), 'Activation'] = 85 # then set class to 1
            
df_new1.loc[((df_new1['Class1']=="WeakDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="StrongLatentDriver"))|((df_new1['Class1']=="StrongLatentDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="WD")), 'Activation'] = 75 # then set class to 1
            
df_new1.loc[((df_new1['Class1']=="WeakDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="WeakLatentDriver"))|((df_new1['Class1']=="WeakLatentDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="WeakDriver")), 'Activation'] = 50 # then set class to 1            
            
    
# driver+wd veya d+d olanlar k
df_new1.loc[((df_new1['Class1']=="WeakLatentDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="WeakLatentDriver")), 'Activation'] = 30 # then set class to 1
    
df_new1.loc[((df_new1['Class1']=="StrongLatentDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="StrongLatentDriver")), 'Activation'] = 30 # then set class to 1
    

df_new1.loc[((df_new1['Class1']=="WeakLatentDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="StrongLatentDriver"))|((df_new1['Class1']=="StrongLatentDriver") & # if discount is more than .2 of total 
       (df_new1['Class2']=="WeakLatentDriver")), 'Activation'] = 30 # then set class to 1

########
######   Filtering 1, number of double mutant 
df_new2 = df_new1[df_new1['DoubleMut#']>=3]  

df_new2.to_csv("./Data/Significant_Double_Mutations.txt",sep="\t",index=False,header=True)                                                 
###########3

"""OG_TSG_Dict={}
with open("./Data/OncoKB_CancerGeneList_OG_TSG_Label.txt","r") as infile:
    for line in infile:
        if "Hugo" not in line:
            splitted=line.rstrip("\n").split("\t")
            OG_TSG_Dict[splitted[0]]=splitted[-1]
OG_TSG_Dict["PIK3CA"]
pickle.dump(OG_TSG_Dict , open("./Data/OG_TSG_Dictionary.p", "wb")) 
    
OG_TSG_Dict        
#get the list of oncogenic drivers        
OncogenicMutations=[]
with open("./Data/Oncogenic_Driver_Mutations_from_CGI.txt","r") as infile:
    for line in infile:
        if  "Residue" not in line:
            splitted=line.rstrip("\n").split("\t")    
            OncogenicMutations.append(splitted[0]+"_"+splitted[-1])
pickle.dump(OncogenicMutations , open("./Data/OncogenicMutations_Dictionary.p", "wb")) """



        