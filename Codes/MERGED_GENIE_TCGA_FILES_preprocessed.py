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
#PATH: ./Data
import pandas as pd
import os

####### Merging point mutations and tumor information of GENIE and TCGA cohortd

import pandas as pd
import os
import pickle

class TCGA_GENIE_MERGED_DATA:
    import pickle
    
    def __init__(self, path, tcga_file_name, genie_file_name, Vaf_Value=0.125): #obtin by writing geni-tcga preprocessed files
        self.path = path
        self.Vaf_Value = Vaf_Value
        
        self.tcga_file_name = tcga_file_name
        self.file_path_tcga = os.path.join(self.path, self.tcga_file_name)
        
        self.genie_file_name = genie_file_name
        self.file_path_genie = os.path.join(self.path, self.genie_file_name)

        self.name = "TCGA-Mutation Annotation File Preprocessing"
        #read preprocessed TCGA and GENIE files from --TCGA,GENIE preprocessing-- with the following columns listed in --SpecificColumns--
        SpecificColumns=['Hugo_Symbol','Tumor_Sample_Barcode',"WT_Residue",'HGVSp_Short', "VAF","Variant_Classification","WildTypeResidue","ResidueNumber",'MutantResidue']
        #select specific columns to read
        self.SpecificColumns = SpecificColumns 
        self.tcga_df = pickle.load(open(self.file_path_tcga , "rb"))
        self.genie_df = pickle.load(open(self.file_path_genie , "rb"))
        self.tcga_df = self.tcga_df[SpecificColumns]
        self.genie_df = self.genie_df[SpecificColumns]
        # self.tcga_df = pd.read_csv( self.file_path_tcga, sep = "\t", usecols = SpecificColumns, comment="#")
        # self.genie_df = pd.read_csv( self.file_path_genie, sep = "\t", usecols = SpecificColumns, comment="#")
        
        
        
    #######CLASS END
    def TCGA_Vaf_Filtering(self):
        
        #df_tcga= pd.read_csv(self.file_path_tcga ,sep="\t",comment='#',usecols=self.SpecificColumns)
        df_tcga = pickle.load(open(self.file_path_tcga , "rb"))
        # drop rows which have same 'Hugo_Symbol','Tumor_Sample_Barcode', 'HGVSp_Short', "VAF","Variant_Classification" and keep last
        new_tcga = df_tcga.drop_duplicates(
          subset = ['Hugo_Symbol','Tumor_Sample_Barcode', 'HGVSp_Short', "VAF","Variant_Classification"],
          keep = 'last').reset_index(drop = True)

        df_tcga1=new_tcga[new_tcga["VAF"] > self.Vaf_Value] 
        df_tcga1.reset_index(inplace=True, drop=True)
        
        self._tcga_df_vaf_filtered = df_tcga1
        return self._tcga_df_vaf_filtered

    def GENIE_Vaf_Filtering(self):
        
        #df_genie= pd.read_csv(self.file_path_genie ,sep="\t",comment='#',usecols=self.SpecificColumns)
        df_genie = pickle.load(open(self.file_path_genie , "rb"))
        # drop rows which have same 'Hugo_Symbol','Tumor_Sample_Barcode', 'HGVSp_Short', "VAF","Variant_Classification" and keep last
        new_genie = df_genie.drop_duplicates(
          subset = ['Hugo_Symbol','Tumor_Sample_Barcode', 'HGVSp_Short', "VAF","Variant_Classification"],
          keep = 'last').reset_index(drop = True)
   
        df_genie1=new_genie[new_genie["VAF"] > self.Vaf_Value] 
        df_genie1.reset_index(inplace=True, drop=True)
        
        self._genie_df_vaf_filtered = df_genie1
        return self._genie_df_vaf_filtered

        
        
     #######CLASS END
        
    def Merge_TCGA_GENIE_VAF_filtered(self):
        tcga_df = self.TCGA_Vaf_Filtering()
        
        genie_df = self.GENIE_Vaf_Filtering()
            
        df_merged=pd.concat([tcga_df, genie_df], ignore_index=True)
        df_merged.reset_index(inplace=True, drop=True)
        self._merged_df = df_merged
        return self._merged_df 
        #############   
        
        ##aşağıdaki kısmı class dışındaki fonksiyonlardan çekmek gerekiyor
        ##############
        ################3
    # def Write_Merged_File_to_csv(self):
    #     df_merged = self.Merge_TCGA_GENIE_VAF_filtered()
    #     df_merged.to_csv("OOP_Generated_Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt",sep="\t",index=False,header=True)
    
    #mutation occurs ~25% of the cells of sequenced tumor 
    def AllTumors(self):
        df_merged = self.Merge_TCGA_GENIE_VAF_filtered()
        total_tumors = set(df_merged["Tumor_Sample_Barcode"].to_list())
        self._total_tumor = total_tumors
        return self._total_tumor
     #in both data 62,567 tumors with point mutations
    # and  VAF>0.125
    def AllGenes(self):
        
        df_merged = self.Merge_TCGA_GENIE_VAF_filtered()
        
        self._all_mutated_genes = set(df_merged['Hugo_Symbol'].to_list()) 
        return self._all_mutated_genes
    
    #mutation list present in at least three tumors
    #"""
    def Merged_TCGA_GENIE_Gene_Mutation_Tumor_Dictionary(self):
        
        df_merged = self.Merge_TCGA_GENIE_VAF_filtered()
        AllGenes = self.AllGenes()
        #['Hugo_Symbol','Tumor_Sample_Barcode', 'HGVSp_Short', "VAF","Variant_Classification"]
        Gene_Mutation_Tumor_dict = { gene:dict() for gene in AllGenes }
        
        for ind in df_merged.index:
            tumor = df_merged.loc[ind,"Tumor_Sample_Barcode"]
            gene = df_merged.loc[ind,"Hugo_Symbol"]
            mutation = df_merged.loc[ind,"WT_Residue"]
            
            if mutation not in Gene_Mutation_Tumor_dict[gene].keys():
                Gene_Mutation_Tumor_dict[gene][mutation]=set()
                Gene_Mutation_Tumor_dict[gene][mutation].add(tumor)
            elif mutation in Gene_Mutation_Tumor_dict[gene].keys():
                Gene_Mutation_Tumor_dict[gene][mutation].add(tumor)
        self._Gene_Mutation_Tumor_dict = Gene_Mutation_Tumor_dict
        pickle.dump(Gene_Mutation_Tumor_dict, open("Merged_GENIE_TCGA_Gene_Mutation_Tumor_Dict.p", "wb"))
        return self._Gene_Mutation_Tumor_dict
    #"""
    def MutationsList_GEQ3(self):
        """gene_mutation observed on at least three or less than three tumors"""
        MutationList_GEQ3=set()
        GeneMutTumor = pickle.load(open("Merged_GENIE_TCGA_Gene_Mutation_Tumor_Dict.p", "rb"))
        for gene in GeneMutTumor.keys():
            for mut in GeneMutTumor[gene].keys():
                if len(GeneMutTumor[gene][mut]) >= 3:
                    MutationList_GEQ3.add(gene+"_"+mut)
               
        self._MutationList_GEQ3 = MutationList_GEQ3
        
        return self._MutationList_GEQ3
                    
    def MutationsList_LEQ2(self):
        """gene_mutation observed on at least three or less than three tumors"""
        MutationList_LEQ2=set()
        GeneMutTumor = pickle.load(open("Merged_GENIE_TCGA_Gene_Mutation_Tumor_Dict.p", "rb"))
        for gene in GeneMutTumor.keys():
            for mut in GeneMutTumor[gene].keys():
                if len(GeneMutTumor[gene][mut]) <= 2:
                    MutationList_LEQ2.add(gene+"_"+mut)
        self._MutationList_LEQ2 = MutationList_LEQ2
        
        return self._MutationList_LEQ2                
                       
    #Tumor Gene_Mutation dictionary, mutations observed on at elat three tumors
    """
    def Merged_TCGA_GENIE_Tumor_GeneMutation_GEQ3_Dictionary(self):
        df_merged = self.Merge_TCGA_GENIE_VAF_filtered()
        Tumors = set(df_merged["Tumor_Sample_Barcode"].to_list())
        
        MUTlist = self.MutationsList_GEQ3()
        
        Tumor_GeneMutation_dict = { tumor:set() for tumor in Tumors }
        
        for ind in df_merged.index:
            tumor = df_merged.loc[ind,"Tumor_Sample_Barcode"]
            gene = df_merged.loc[ind,"Hugo_Symbol"]
            mutation = df_merged.loc[ind,"WT_Residue"]
            
            genemutation = gene+"_"+mutation 
            
            if genemutation in MUTlist:
                Tumor_GeneMutation_dict[tumor].add(genemutation)
        pickle.dump(Tumor_GeneMutation_dict, open("Merged_TCGA_GENIE_Tumor_GeneMutation_GEQ3_Dictionary.p", "wb"))
        return Tumor_GeneMutation_dict
        """        
        
        # #['Hugo_Symbol','Tumor_Sample_Barcode', 'HGVSp_Short', "VAF","Variant_Classification"]
        # Tumor_Gene_Mutation_dict = { tumor:dict() for tumor in Tumors }
        # MUTlist = self.MutationsList_GEQ3()
        # for ind in df_merged.index:
        #     tumor = df_merged.loc[ind,"Tumor_Sample_Barcode"]
        #     gene = df_merged.loc[ind,"Hugo_Symbol"]
        #     mutation = df_merged.loc[ind,"WT_Residue"]
        #     if gene+"_"+mutation in MUTlist: ###!!!!!! bunu fonksiyonu çağırmadan yapabiliyor muyuz bilmiyorum
        #         if gene not in Tumor_Gene_Mutation_dict[tumor].keys():
        #             Tumor_Gene_Mutation_dict[tumor][gene] = set()
        #             Tumor_Gene_Mutation_dict[tumor][gene].add(mutation)
        #         elif gene in Tumor_Gene_Mutation_dict[tumor].keys():
        #             Tumor_Gene_Mutation_dict[tumor][gene].add(mutation)
        # self._Tumor_Gene_Mutation_dict = Tumor_Gene_Mutation_dict
        # return self._Tumor_Gene_Mutation_dict
    
    # def TumorSpecific_PotentialDoubleMutationdictionary(self):
    #     Tumor_Gene_Mutation_Dict = self.Merged_TCGA_GENIE_Tumor_Gene_Mutation_GEQ3_Dictionary()
        
        
    #     Tumor_Gene_Doublet_Dict = {tumor:dict() for tumor in Tumor_Gene_Mutation_Dict.keys()}
    #     for tumor in Tumor_Gene_Doublet_Dict.keys():
    #         if Tumor_Gene_Doublet_Dict[tumor]!={}:
    #             for gene in Tumor_Gene_Mutation_Dict[tumor].keys():
    #                 if len(Tumor_Gene_Doublet_Dict[tumor][gene])>=2:
    #                     #Tumor_Gene_Doublet_Dict[tumor][gene] = set()
    #                     MutList = list(Tumor_Gene_Mutation_Dict[tumor][gene])
    #                     Tumor_Gene_Doublet_Dict[tumor][gene] = MutList
    #                     # for i in range(0,len(MutList)-1):
    #                     #     for j in range(i+1,len(MutList)):
    #                     #         Tumor_Gene_Doublet_Dict[tumor][gene].add((MutList[i],MutList[j]))
                       
    #     self._Tumor_Gene_Doublet_Dict = Tumor_Gene_Doublet_Dict
    #     return self._Tumor_Gene_Doublet_Dict
                            
# mut_file_tcga = "OOP_Generated_TCGA_MC3callSet_Mutation_Cleaned_Parsed_File__new.txt" 
# mut_file_genie = "OOP_Generated_GENIE_vol6_Mutation_Cleaned_Parsed_File__new.txt"            
# path = "/Users/bengi/Desktop/CommsBioRevision_Final_9Dec/REVISION#2/DATA_CODES_new/OOP_OUTPUT"       
# merged = TCGA_GENIE_MERGED_DATA(path,mut_file_tcga, mut_file_genie )  

# #merged.MutationsList_GEQ3()

# x=merged.Merged_TCGA_GENIE_Tumor_GeneMutation_GEQ3_Dictionary()
# #merged.TotalTumorCount()    
# #need to call this to write preprocessed and parsed file to file

# len(merged.MutationsList_GEQ3())
# #mut_dict=merged.Merged_TCGA_GENIE_Tumor_Gene_Mutation_GEQ3_Dictionary()
# muts=set()
# with open("/Users/bengi/Desktop/CommsBioRevision_Final_9Dec/REVISION#2/DATA_CODES_new/Gene_Mutation_TumorList__geq3_new.txt","r") as infile:
#     for line in infile:
#         splitted=line.rstrip("\n").split("\t")
#         muts.add(splitted[0]+"_"+splitted[1])
# len(muts)        
# #len(merged.AllGenes())

#         #Form the potential doublets, be careful about (a,b)-(b,a) type duplicates
#         Doublets=[]
#         #2230203, doublets are formed from mutations observed on >=tumors
#         for gene in Genes_MutationFor_Combinations:
#             for i in range(0,len(Genes_MutationFor_Combinations[gene])-1):
#                 for j in range(i+1,len(Genes_MutationFor_Combinations[gene])):
#                     Doublets.append((gene+"_"+Genes_MutationFor_Combinations[gene][i],gene+"_"+Genes_MutationFor_Combinations[gene][j]))
#         len(Doublets)
        
            
# #Genes=merged.AllGenes()           

# Possible=set()
# for tumor in x.keys():
#     for mut1 in x[tumor]:
#         for mut2 in x[tumor]:
#             if (mut1!=mut2) and (mut1.split("_")[0]==mut2.split("_")[0]):
#                 if ((mut1,mut2) not in Possible) and ((mut2,mut1) not in Possible):
#                     Possible.add((mut1,mut2))
                
        
            
            
        
# DF=merged.tcga_df

# import pickle
# pickle.dump(DF, open("TCGA_df.p", "wb"))    
 
# network_dict = pickle.load(open("TCGA_df.p", "rb"))
  
# def read_patient_networks():
#     network_dict = pickle.load(open("TCGA_df.p", "rb"))
#     return network_dict

# # pickle.dump(merged,open("TCGA_GENIE_merged_object.p", "wb"))


# merged_pickled = pickle.load(open("TCGA_GENIE_merged_object.p", "rb"))




#         #for each gene write mutations observed on >=3 tumors to file with "," separated tumor list and its size
# # Genes_Sorted=sorted(list(Genes)) #19415 genes
# # ###write gene, mutation, tumors, and tumor_count to file
# # for gene in Genes_Sorted:
# #     for mutation in Gene_Mutation_Sample_dict[gene].keys():
# #         if len(Gene_Mutation_Sample_dict[gene][mutation])>=3: #use as paramete, initially we used 5 as a threshold
# #             AllTumor=",".join(list(Gene_Mutation_Sample_dict[gene][mutation]))
# #             TumorCount=str(len(Gene_Mutation_Sample_dict[gene][mutation]))
# #             try:
# #                 with open("Gene_Mutation_TumorList__geq3_new.txt","a") as outfile:
# #                     outfile.write(str(gene)+"\t"+str(mutation)+"\t"+str(TumorCount)+"\t"+AllTumor+"\n")
# #             except TypeError:
# #                 print((gene, mutation))

# # #####
# # ###write gene, mutation, tumors, and tumor_count to file for the mutations observed on 2 tumors
# # for gene in Genes_Sorted:
# #     for mutation in Gene_Mutation_Sample_dict[gene].keys():
# #         if len(Gene_Mutation_Sample_dict[gene][mutation])<=2: #mutations observed on less than 2 tumors
# #             AllTumor=",".join(list(Gene_Mutation_Sample_dict[gene][mutation]))
# #             TumorCount=str(len(Gene_Mutation_Sample_dict[gene][mutation]))
# #             try:
# #                 with open("Gene_Mutation_TumorList__LEQ_2_new.txt","a") as outfile:
# #                     outfile.write(str(gene)+"\t"+str(mutation)+"\t"+str(TumorCount)+"\t"+AllTumor+"\n")
# #             except TypeError:
# #                 print((gene, mutation))

#  #   
# #Form the potential doublets from the mutations observed on >=3 tumors
# #FolderName: /Users/bengi/Desktop/CommsBioRevision_Final_9Dec/REVISION#2/DATA_CODES_new
# Gene_Mutation_List={} #each gene keep track of mutations with #tumors>=3 --12,724 key genes with 
# #at least one mutation observed on more than 2 tumors
# GeneMut_TumorDict={} #{Gene_Mutation:Tumor_set} dictionary, 65,872 mutations in total
# with open("Gene_Mutation_TumorList__geq3_new.txt","r") as infile: 
#     for line in infile:
#         splitted=line.rstrip("\n").split("\t")
#         gene=splitted[0]
#         mutation=splitted[1]
#         genemut=gene+"_"+mutation
#         patient=set(splitted[-1].split(","))
#         GeneMut_TumorDict[genemut]=patient
#         if gene not in Gene_Mutation_List.keys():
#             Gene_Mutation_List[gene]=[]
#             Gene_Mutation_List[gene].append(mutation)
#         elif gene in Gene_Mutation_List.keys():
#             Gene_Mutation_List[gene].append(mutation)
# len(GeneMut_TumorDict) 
# len(Gene_Mutation_List)          
# ##Genes with 2 or more mutations that or observed on at least 3 tumors
# #we need at least 2 mutations on a gene to form potential doublets
# Genes_MutationFor_Combinations={}   #8189 genes with >=2 mutations
# for gene in Gene_Mutation_List.keys():
#     if len(Gene_Mutation_List[gene]) >=2:
#         Genes_MutationFor_Combinations[gene]=Gene_Mutation_List[gene]
# len(Genes_MutationFor_Combinations["PIK3CA"])   
# len(Genes_MutationFor_Combinations)   
# Genes_MutationFor_Combinations
# #
# ####
# #Form the potential doublets, be careful about (a,b)-(b,a) type duplicates
# Doublets=[]
# #2230203, doublets are formed from mutations observed on >=tumors
# for gene in Genes_MutationFor_Combinations:
#     for i in range(0,len(Genes_MutationFor_Combinations[gene])-1):
#         for j in range(i+1,len(Genes_MutationFor_Combinations[gene])):
#             Doublets.append((gene+"_"+Genes_MutationFor_Combinations[gene][i],gene+"_"+Genes_MutationFor_Combinations[gene][j]))
# len(Doublets)
# Doublets
# #19,029 of the doublets observed on at least one tumor (intersection >0)
# Doublets_Nonzero=[] 
# for doublet in Doublets:
#     if len(GeneMut_TumorDict[doublet[0]].intersection(GeneMut_TumorDict[doublet[1]]))>0:
#         Doublets_Nonzero.append(doublet)
# len(Doublets_Nonzero)
# Doublets_Nonzero
# ##########
# ##########
# #the doublets matching with our initial 228 doublets        
# df_original_doublets=pd.read_excel("/Users/bengi/Desktop/CommsBioRevision_Final_9Dec/Supplementray_Table_1.xlsx",sheet_name="Statistics_")
# OriginalDoublets=[]
# df_original_doublets.columns
# for inds in df_original_doublets.index:
#     mut1=df_original_doublets.iloc[inds]['Gene']+"_"+df_original_doublets.iloc[inds]['Residue_1']
#     mut2=df_original_doublets['Gene'].iloc[inds]+"_"+df_original_doublets.iloc[inds]['Residue_2']
#     OriginalDoublets.append((mut1,mut2))
# common=[] #146 of the original doublets exist now 
# for double in Doublets_Nonzero:
#     m1=double[0]
#     m2=double[1]
#     if ((m1,m2) in OriginalDoublets) or ((m2,m1) in OriginalDoublets):
#         common.append(double)

# #among 228 doublets which not exist now, eaxclude while uploading
# ###############
# not_common=[]
# for double in OriginalDoublets:
#     m1=double[0]
#     m2=double[1]
#     if ((m1,m2) not in common) and ((m2,m1) not in common):
#         not_common.append(double)
        
#         #the doublets in the following file contains 82 doublets out of 228 initial significant doublets
#         #now either the components' intersection is empty tumor set
#         #or at least one component observed on less than 3 tumors, hence not taken into account
# for double in not_common:
#     m1=double[0]
#     m2=double[1]
#     try:
#         count=len(GeneMut_TumorDict[m1].intersection(GeneMut_TumorDict[m2]))
#         with open("In228butNotExistNow_geq3_new.txt","a") as outfile:
#             outfile.write(m1+"\t"+m2+"\t"+str(count)+"\n")
#     except KeyError:
#         with open("In228butNotExistNow_geq3_new.txt","a") as outfile:
#             outfile.write(m1+"\t"+m2+"\t"+"AtLeastOneComponentNotExist"+"\n")
#         print(double)
 
# #GeneMut_TumorDict["BRAF_I710"]
# #GeneMut_TumorDict['TSC1_E878'].intersection(GeneMut_TumorDict['TSC1_R892'])
# #which original doublets not exist
# #this part could be excluded
# ######################
        
# ###
# #Obtain the number of gene mutant tumors (mutations observed on >=3 tumors)
# Gene_MutantTumors_Dict={}
# with open("Gene_Mutation_TumorList__geq3_new.txt","r") as infile: 
#     for line in infile:
#         splitted=line.rstrip("\n").split("\t")
#         gene=splitted[0]
#         patient=set(splitted[-1].split(","))
#         if gene not in Gene_MutantTumors_Dict.keys():
#             Gene_MutantTumors_Dict[gene]=set()
#             Gene_MutantTumors_Dict[gene]=Gene_MutantTumors_Dict[gene].union(patient)
#         elif gene in Gene_MutantTumors_Dict.keys():
#             Gene_MutantTumors_Dict[gene]=Gene_MutantTumors_Dict[gene].union(patient)
# len(Gene_MutantTumors_Dict["PIK3CA"] )           
# ### 
# ###
# #write the statistics for the potential doublets observed on >0 tumors without any further filtering
#             #after this step read the file and compute p values for 4 different contingency tables
# with open("SameGeneDoubletsTumorCountsNoFilter__geq3_new.txt","a") as outfile:
#     outfile.write("Gene"+"\t"+"Mut1"+"\t"+"Mut2"+"\t"+"Among228Doubles"+"\t"+"DoubleMut#"+"\t"+"OnlyMut1"+"\t"+"OnlyMut2"+"\t"+"GeneMutant#"+"\n")

# for doublet in Doublets_Nonzero:
#     if doublet in common:
#         tag="Yes"
#     else:
#         tag="No"
#     gene=doublet[0].split("_")[0]
#     mut1=doublet[0]
#     mut2=doublet[1]
#     doubleMutant=GeneMut_TumorDict[mut1].intersection(GeneMut_TumorDict[mut2])
#     only1=GeneMut_TumorDict[mut1].difference(GeneMut_TumorDict[mut2])
#     only2=GeneMut_TumorDict[mut2].difference(GeneMut_TumorDict[mut1])
#     geneMutant=Gene_MutantTumors_Dict[gene]
#     with open("SameGeneDoubletsTumorCountsNoFilter__geq3_new.txt","a") as outfile:
#         outfile.write(gene+"\t"+mut1.split("_")[1]+"\t"+mut2.split("_")[1]+"\t"+tag+"\t"+str(len(doubleMutant))+"\t"+\
#                      str(len(only1))+"\t"+str(len(only2))+"\t"+str(len(geneMutant))+"\n")
# #SameGeneDoubletsTumorCountsNoFilter__.txt does not contains nan values on the columns mut1, mut2     
# GeneMut_TumorDict.keys()        
# import pandas as pd
# from scipy.stats import fisher_exact
# import numpy as np

# df=pd.read_csv("SameGeneDoubletsTumorCountsNoFilter__geq3_new.txt",sep="\t")


# #first write header
# with open("SameGeneStatisticsAllDoublets_VariousContingencyTables__geq3_new.txt","a") as outfile:
#     outfile.write("Gene"+"\t"+"Mut1"+"\t"+"Mut2"+"\t"+"Among228Doubles"+"\t"+"DoubleMut#"+"\t"+"OnlyMut1"+"\t"+"OnlyMut2"+"\t"+"GeneMutant#"+"\t"+\
#                   "AllTumor#VAF>"+"\t"+"p1"+"\t"+"p2"+"\t"+"p3"+"\t"+"OddsR3"+"\t"+"p4"+"\t"+"OddsR4"+"\n")
# #calculate p-values with fisher exact test for different contingency tables
# for i in df.index:
#     tag=df.iloc[i]['Among228Doubles']
#     a=df.iloc[i]['DoubleMut#']
#     b=df.iloc[i]["OnlyMut1"]
#     c=df.iloc[i]["OnlyMut2"]
#     d1=b+c
#     d2=0
#     d3=df.iloc[i]['GeneMutant#']-(a+b+c)
#     d4=AllTumorCount-(a+b+c)
#     table1=np.array([[a,b],[c,d1]]) #values are dependentt
#     table2=np.array([[a,b],[c,d2]]) #values indep., evaluated only among A or B mutants
#     table3=np.array([[a,b],[c,d3]]) #values indep., evaluated  among gene mutants
#     table4=np.array([[a,b],[c,d4]]) #values indep., evaluated  among all tumors with vaf>0.125
#     oddsr1, p1 = fisher_exact(table1, alternative='two-sided')
#     oddsr2, p2 = fisher_exact(table2, alternative='two-sided')
#     oddsr3, p3 = fisher_exact(table3, alternative='two-sided')
#     oddsr4, p4 = fisher_exact(table4, alternative='two-sided')
#     try:
#         with open("SameGeneStatisticsAllDoublets_VariousContingencyTables__geq3_new.txt","a") as outfile:
#             outfile.write(df.iloc[i]["Gene"]+"\t"+df.iloc[i]['Mut1']+"\t"+df.iloc[i]['Mut2']+"\t"+tag+"\t"+str(a)+"\t"+str(b)+"\t"+str(c)+"\t"+str(df.iloc[i]['GeneMutant#'])\
#                           +"\t"+str(AllTumorCount)+"\t"+str(p1)+"\t"+str(p2) +"\t"+str(p3)+"\t"+ str(oddsr3)+"\t"+str(p4)+"\t"+ str(oddsr4) +"\n")
#     except IndexError:
#         print(i)
# ###double mutant tumor list        
# for doublet in Doublets_Nonzero:
#     if doublet in common:
#         tag="Yes"
#     else:
#         tag="No"
#     gene=doublet[0].split("_")[0]
#     mut1=doublet[0]
#     mut2=doublet[1]
#     doubleMutant=GeneMut_TumorDict[mut1].intersection(GeneMut_TumorDict[mut2])
#     Tumors=",".join(doubleMutant)     
#     with open("AllDoubletsObservedAtLeastOneTumorDoubleMutantTumorsInfo_new.txt","a") as outfile:
#         outfile.write(gene+"\t"+mut1+"\t"+mut2+"\t"+tag+"\t"+str(len(doubleMutant))+"\t"+Tumors+"\n")
        
