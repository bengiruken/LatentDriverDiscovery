#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: bengi
"""
##########
#Upload necessary files and dictionaries from preprocessed pickle files of TCGA and GENIE cohorts

PATH = "./Data"

tcga_file_name  = "Preprocessed_TCGA_df.p"
genie_file_name = "Preprocessed_GENIE_df.p"
import pickle

from MERGED_GENIE_TCGA_FILES_preprocessed import TCGA_GENIE_MERGED_DATA

merged_obj = TCGA_GENIE_MERGED_DATA(PATH, tcga_file_name, genie_file_name, Vaf_Value=0.125)
merged_df = merged_obj.Merge_TCGA_GENIE_VAF_filtered()
merged_df.head()
all_tumors = merged_obj.AllTumors()
Mutation_List_GEQ3 = merged_obj.MutationsList_GEQ3()


gene_mut_tumor_dict =pickle.load(open(PATH+"Gene_Mutation_Tumor_GEQ3_Dictionary.p","rb"))


Gene_Mutation_Tumor_Dictionary=pickle.load(open(PATH+"Merged_GENIE_TCGA_Gene_Mutation_Tumor_Dict.p", "rb"))    #merged_obj.Merged_TCGA_GENIE_Tumor_GeneMutation_GEQ3_Dictionary()

Gene_All_Mutants_Dictionary =  pickle.load(open(PATH+"Gene_All_Mutants_Dictionary.p", "rb")) 
############################


GeneList = sorted(list(set([mut.split("_")[0] for mut in gene_mut_tumor_dict.keys()])))
#find the list of all potential double mutations

def Gene_Mutation_Dictionary(geneList = GeneList,gene_mut_tumor_dict=gene_mut_tumor_dict):
    Gene_Mutation_Dict = {gene:set() for gene in geneList}
    
    for mut in gene_mut_tumor_dict.keys():
        #print(mut)
        gene = mut.split("_")[0]
        Gene_Mutation_Dict[gene].add(mut)
    return Gene_Mutation_Dict

Gene_Mutation_Dictionary()
###########
################Get the list of potential double mutations on each gene
def All_Potential_Doubles(geneList = GeneList):
    GeneMutationDict = Gene_Mutation_Dictionary()
    all_Doubles = set()
    Count=0
    for gene in geneList:
        for mut1 in GeneMutationDict[gene]:
            for mut2 in GeneMutationDict[gene]:
                if ((mut1,mut2) not in all_Doubles) and (((mut2,mut1) not in all_Doubles)):
                    if mut1!=mut2:
                        all_Doubles.add((mut1,mut2))
                        Count+=1
    print(Count) 
    return  sorted(list(all_Doubles)) 
#All_Doubles  = All_Potential_Doubles()  
########### 

def Fisher_Exact_Test(All_tumors=all_tumors):
    #FISHER EXACT TEST AMONG ALL TUMORS with at least one mutation of VAF>0.125
    from scipy.stats import fisher_exact
    import numpy as np
    All_Doubles  = All_Potential_Doubles()  
    #LOAD the pickle files
    import pickle
   
    #####################
    Gene_Mutation_Tumor_Dictionary=pickle.load(open(PATH+"Merged_GENIE_TCGA_Gene_Mutation_Tumor_Dict.p", "rb"))    #merged_obj.Merged_TCGA_GENIE_Tumor_GeneMutation_GEQ3_Dictionary()
    Gene_All_Mutants_Dictionary =  pickle.load(open(PATH+"Gene_All_Mutants_Dictionary.p", "rb")) 
    ############################
    #first write header
    with open(PATH+"Potential_SameGeneDoubles_Fisher_Exact_Test.txt","a") as outfile:
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
        #all_tumors = merged_obj.AllTumors()
        
        d4=len(All_tumors)-(a+b+c) #len(all_tumors)
        
        table4=np.array([[a,b],[c,d4]]) #values indep., evaluated  among all tumors with vaf>0.125

        oddsr4, p4 = fisher_exact(table4, alternative='two-sided')

        with open(PATH+"Potential_SameGeneDoubles_Fisher_Exact_Test.txt","a") as outfile:
            outfile.write(gene+"\t"+mut1+"\t"+mut2+"\t"+str(a)+"\t"+str(b)+"\t"+str(c)+"\t"+str(gene_all_mutants)+"\t"+str(len(All_tumors))+"\t"+str(p4)+"\t"+ str(oddsr4) +"\n")
    #################3
Fisher_Exact_Test()


#FDR corrected p val with benjamini hochberg method
#!pip install statsmodels
def FDR_Correction():
    import statsmodels.stats.multitest
    import pandas as pd
        
    import numpy as np
    
    df_stats=pd.read_csv(PATH+"Potential_SameGeneDoubles_Fisher_Exact_Test.txt",sep="\t")
    df_stats.columns
    print(len(df_stats))
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
        
    
    df_stats.to_csv(PATH+"Potential_SameGeneDoubles_Fisher_Exact_Test_FDR_Multiple_Testing.txt",sep="\t",index=False,header=True)
    

    # ####
    # writer = ExcelWriter(PATH+'Potential_SameGeneDoubles_Fisher_Exact_Test_FDR_Multiple_Testing.xlsx', engine='xlsxwriter')
    # df_stats.to_excel(writer,'Statistics_AmongAllTumors',index=False,header=True)
    # # Close the Pandas Excel writer and output the Excel file.
    # writer.save()

"""return p and q values of each potential doublet with fisher exact test and benjamini-hochberg, respectively   """
FDR_Correction()



############
#
"""Filtering with q<0.1, and number of double mutatn tumors >=3"""
import pandas as pd
df_all = pd.read_csv(PATH+"Potential_SameGeneDoubles_Fisher_Exact_Test_FDR_Multiple_Testing.txt",sep="\t")
df_all.columns




df_significant = df_all[df_all['Qval4']<0.1]
df_significant.to_csv(PATH+"SignificantDoubles_newStatistics_FDR_q<01.txt",sep="\t",index=False)

df_significant1 = df_significant[df_significant['DoubleMut#']>=3]


"""Filtering with VAF values and variant classess of each double mutation constituent"""

##############
"""Double mutations-VAF Filtering """
#1.prepare a file with mut1 mut2 tumor vaf values and mut class
merged_df.head()
merged_df.columns
#######change the column names
Mutation_Tumor_VAF_Dict = {mut:{} for mut in gene_mut_tumor_dict.keys()}
Mutation_Tumor_VariantClass_Dict = {mut:{} for mut in gene_mut_tumor_dict.keys()}

merged_df.reset_index(inplace=True,drop=True)
for j in merged_df.index:
    mut = merged_df.loc[j,"Hugo_Symbol"] + "_" + merged_df.loc[j,"WT_Residue"]
    if mut in gene_mut_tumor_dict.keys():
        tumor = merged_df.loc[j,"Tumor_Sample_Barcode"]
        vaf = merged_df.loc[j,"VAF"]
        variant = merged_df.loc[j,"Variant_Classification"]
        
        Mutation_Tumor_VAF_Dict[mut][tumor] = vaf
        Mutation_Tumor_VariantClass_Dict[mut][tumor] = variant
   


#####Get the significant doubles, q<0.1 and add vaf and variant class columns
df_significant = df_all[df_all['Qval4']<0.1]
df_significant.reset_index(drop=True,inplace=True)


with open(PATH+"Double_Mutants_VAF_total_VariantClass_For_Filtering.txt","a") as outfile:
    outfile.write("Mut1"+"\t"+"Mut2"+"\t"+"Tumor"+"\t"+"VAF1"+"\t"+"VAF2"\
        +"\t"+"VarClass1"+"\t"+"VarClass2"+"\n")


for ind in df_significant.index:
    mut1 = df_significant.loc[ind,"Gene"]+"_"+df_significant.loc[ind,"Mut1"]
    mut2 = df_significant.loc[ind,"Gene"]+"_"+df_significant.loc[ind,"Mut2"]
    
    DoubleMutant = gene_mut_tumor_dict[mut1].intersection(gene_mut_tumor_dict[mut2])
    if len(DoubleMutant)>0:
        for tumor in DoubleMutant:
            with open(PATH+"Double_Mutants_VAF_total_VariantClass_For_Filtering.txt","a") as outfile:
                outfile.write(mut1+"\t"+mut2+"\t"+tumor+"\t"+str(Mutation_Tumor_VAF_Dict[mut1][tumor])+"\t"+str(Mutation_Tumor_VAF_Dict[mut2][tumor])\
                    +"\t"+Mutation_Tumor_VariantClass_Dict[mut1][tumor]+"\t"+Mutation_Tumor_VariantClass_Dict[mut2][tumor]+"\n")
            
    

###########
#Doublets where at least one component is nonsense and their fraction

"""FILTERING NONSNESE MUTATIONS"""
## Doublets with nonsense mutation constituent (at least one)
import pandas as pd
df = pd.read_csv(PATH+"Double_Mutants_VAF_total_VariantClass_For_Filtering.txt",sep="\t")
df["Nonsense_Status"]=[None]*len(df)
df["Doublet"] = df["Mut1"]+"+"+df["Mut2"]
for ind in df.index:
    doublet=df.loc[ind,"Doublet"]
    var1=df.loc[ind,"VarClass1"]
    var2=df.loc[ind,"VarClass2"]
    if (var2=="Nonsense_Mutation") and (var1=="Nonsense_Mutation"):
        df.loc[ind,"Nonsense_Status"]="both_nonsense"
    elif (var1!=var2) and ((var2=="Nonsense_Mutation") or (var1=="Nonsense_Mutation")):
        df.loc[ind,"Nonsense_Status"]="one_nonsense"
    elif (var2!="Nonsense_Mutation") and (var1!="Nonsense_Mutation"):
        df.loc[ind,"Nonsense_Status"]="no_nonsense"
    
df.columns
############3
### both nonsense/ one nonsense fraction

Doublets=set(df["Doublet"].to_list()) 

df["Percentage_Both_Nonsense"] =[None]*len(df)  
##############3
for double in Doublets:
    same=0
    List=[]
    total=len(df[df["Doublet"]==double])
    for ind in df.index:
        if (df.loc[ind,"Doublet"]==double):
            List.append(ind)
            if (df.loc[ind,"Nonsense_Status"]=="both_nonsense"):
                same+=1
    for i in List:
        df.loc[i,"Percentage_Both_Nonsense"]=(same/total)*100
        
###
##
"""Filter doublets if nonsense mutation is in downstream from another mutation """
 

df["Nonsense_Position"]=["None"]*len(df)
for ind in df.index:
    doublet=df.loc[ind,"Doublet"]
    var1=df.loc[ind,"VarClass1"]
    var2=df.loc[ind,"VarClass2"]
    mut1=df.loc[ind,"Mut1"].split("_")[1]
    mut2=df.loc[ind,"Mut2"].split("_")[1]
    res1=int(mut1[1:])
    res2=int(mut2[1:])
    if (var1=="Nonsense_Mutation") and (var2!="Nonsense_Mutation"):
        if res1<res2:
            print((res1,res2))
            df.loc[ind,"Nonsense_Position"]="Smaller"
        else:
            df.loc[ind,"Nonsense_Position"]="Larger"
    elif (var1!="Nonsense_Mutation") and (var2=="Nonsense_Mutation"):
        if res2<res1:
            df.loc[ind,"Nonsense_Position"]="Smaller"
        else:
            df.loc[ind,"Nonsense_Position"]="Larger"
            
#if position of nonsense mutation smaller           
df["Percentage_One_Nonsense"] =[None]*len(df)  
##############3
for double in Doublets:
    same=0
    List=[]
    total=len(df[df["Doublet"]==double])
    for ind in df.index:
        if (df.loc[ind,"Doublet"]==double):
            List.append(ind)
            if (df.loc[ind,"Nonsense_Position"]=="Smaller"):
                same+=1
    for i in List:
        df.loc[i,"Percentage_One_Nonsense"]=(same/total)*100


###########
        
##########
df["Percentage_No_Nonsense"] =[None]*len(df)  
##############3
for double in Doublets:
    same=0
    List=[]
    total=len(df[df["Doublet"]==double])
    for ind in df.index:
        if (df.loc[ind,"Doublet"]==double):
            List.append(ind)
            if (df.loc[ind,"Nonsense_Status"]=="no_nonsense"):
                same+=1
    for i in List:
        df.loc[i,"Percentage_No_Nonsense"]=(same/total)*100
########

#Filter with respect to VAF
df["VAF_total"] = df["VAF1"]+df["VAF2"]

df["Tag"]=[None]*len(df)

for ind in df.index:
    if df.loc[ind,"VAF_total"]>0.5:
        df.loc[ind,"Tag"]="SameCell"
    else:
        df.loc[ind,"Tag"]="DifferentCell"
        

df["Percentage_Same"] =[None]*len(df)  
##############3
for double in Doublets:
    same=0
    List=[]
    total=len(df[df["Doublet"]==double])
    for ind in df.index:
        if (df.loc[ind,"Doublet"]==double):
            List.append(ind)
            if (df.loc[ind,"Tag"]=="SameCell"):
                same+=1
    for i in List:
        df.loc[i,"Percentage_Same"]=(same/total)*100

#df.to_csv(PATH+"Rev#5_Double_Mutants_NonsenseComponents_VariantClass_Total_VAF_For_Filtering.txt",sep="\t",index=False)



#Filtering part
########
import pandas as pd

df1 = pd.read_csv(PATH+"RDouble_Mutants_NonsenseComponents_VariantClass_Total_VAF_For_Filtering.txt",sep="\t")   
df1.columns
#
df1=df1[['Doublet','Percentage_Same','Percentage_Both_Nonsense', 'Percentage_One_Nonsense','Percentage_No_Nonsense']]
###

df2=df1.drop_duplicates(keep="first",ignore_index=True)
#########3

Both_Nonsense=set()
for i in df2.index:
    if df2.loc[i,"Percentage_Both_Nonsense"]>=50:
        Both_Nonsense.add(df2.loc[i,"Doublet"])
        
g = set([gene.split("_")[0] for gene in Both_Nonsense ])
###########
df3=df2[df2["Percentage_Both_Nonsense"]<=50]        
################
One_Nonsense=set() #nonsense position is samller in more than half of tumors
for i in df3.index:
    if df3.loc[i,"Percentage_One_Nonsense"]>=50:
        One_Nonsense.add(df3.loc[i,"Doublet"])
###########
df4=df3[df3["Percentage_One_Nonsense"]<50]    
df5 = df4[df4["Percentage_Same"]>=20]




######################3
##Remaining Doublets after filtering
Filtered_Doublets=set(df5["Doublet"].to_list()) 

df_significants=pd.read_csv(PATH+"Rev#5_SignificantDoubles_newStatistics_FDR_q<01.txt",sep="\t")

df_significants["M1"]=df_significants["Gene"]+"_"+df_significants["Mut1"]
df_significants["M2"]=df_significants["Gene"]+"_"+df_significants["Mut2"]

df_significants["Mut1+Mut2"] = df_significants["M1"]+"+"+df_significants["M2"]
df_significants["Mut2+Mut1"] = df_significants["M2"]+"+"+df_significants["M1"]

df_sign_filtered=df_significants[(df_significants["Mut1+Mut2"].isin(Filtered_Doublets))|(df_significants["Mut2+Mut1"].isin(Filtered_Doublets))]

df_sign_filtered.to_csv(PATH+"Significant_doubles_VAF_20percent_AND_Nonsense_filtered.txt",index=False,sep="\t")

 
