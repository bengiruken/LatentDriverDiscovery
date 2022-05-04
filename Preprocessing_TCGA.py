#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 21:45:05 2022

@author: bengi
"""

import pandas as pd
import os

class TCGA_MAF:
    
    def __init__(self, path, file_name):
        self.path = path
        self.file_name = file_name
        
        self.file_path = os.path.join(self.path, self.file_name)

        self.name = "TCGA-Mutation Annotation File Preprocessing"
        #read a MAF file with the following columns listed in --SpecificColumns--
        SpecificColumns = ['Hugo_Symbol','Center','Variant_Classification', 'Variant_Type', "Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2", \
                         'Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode', 'HGVSp_Short',\
                         "t_ref_count","t_alt_count","Codons","SWISSPROT","t_depth"] #select specific columns to read
        
        self.SpecificColumns = SpecificColumns   
        self.mutation_df = pd.read_csv( self.file_path, sep = "\t", usecols = SpecificColumns, comment="#")
    
    # def Convert_Maf_to_Txt(self):
    #     """convert maf file to txt file"""
    #     maf_file = os.path.join(self.path, "mc3.v0.2.8.PUBLIC.maf")
    #     with open(maf_file, 'rb') as fd: #3600964 line in the file, wc -l mc3.v0.2.8.PUBLIC.maf
    #         d = pd.read_csv(fd, sep="\t",usecols= self.SpecificColumns)
    #     d.to_csv("./OOP_OUTPUT/MC3_TCGA_Mutations.txt",sep="\t",header=True,index=False) 
    
    def GetTumorSamples(self):
        """Select the rows where tumor barcode ends with -01"""
        df_mutation = self.mutation_df
        df_mutation = df_mutation[df_mutation["HGVSp_Short"].notnull()==True] # tumors with nonnull hgvsp-short
        
        df_mutation.reset_index(inplace=True,drop=True)
        
        #keep only tumors 
        df_mutation['Tumor_Sample_Barcode']=df_mutation['Tumor_Sample_Barcode'].astype(str).str[:15]
        
        len(set(df_mutation["Tumor_Sample_Barcode"].to_list())) #10295 sample barcodes ends with 01,03...
        
        df_mutation=df_mutation[df_mutation['Tumor_Sample_Barcode'].str.endswith("-01")]
        len(set(df_mutation["Tumor_Sample_Barcode"].to_list())) #9718 tumors ends with 01
        self._tumors_samples = list(set(df_mutation["Tumor_Sample_Barcode"].to_list())) #9718 tumors ends with 01
        df_mutation.reset_index(inplace=True,drop=True)
        return self._tumors_samples
    
    def LengthTumors(self):
        TumorsList = self.GetTumorSamples()
        self._LengthTumorList = len(TumorsList)
        return self._LengthTumorList
    def PreprocessingTcgaMutationFile(self):
        df_mutation = self.mutation_df
        
        len(set(df_mutation["Tumor_Sample_Barcode"].to_list()))  #10295 tumors
                                
        #keep non_null values for "HGVSp_Short" column
        
        df_mutation = df_mutation[df_mutation["HGVSp_Short"].notnull()==True] # tumors with nonnull hgvsp-short
        
        df_mutation.reset_index(inplace=True,drop=True)
        len(set(df_mutation["Tumor_Sample_Barcode"].to_list())) 
        
        
        #keep only tumors 
        df_mutation['Tumor_Sample_Barcode']=df_mutation['Tumor_Sample_Barcode'].astype(str).str[:15]
        
        len(set(df_mutation["Tumor_Sample_Barcode"].to_list())) #10295 sample barcodes ends with 01,03...
        
        df_mutation=df_mutation[df_mutation['Tumor_Sample_Barcode'].str.endswith("-01")]
        len(set(df_mutation["Tumor_Sample_Barcode"].to_list())) #9718 tumors ends with 01
        
        df_mutation=df_mutation[df_mutation['HGVSp_Short']!="."]
        df_mutation=df_mutation[df_mutation['HGVSp_Short']!="?"]
        len(set(df_mutation["Tumor_Sample_Barcode"].to_list())) #9716 tumors 
        set(df_mutation["Tumor_Sample_Barcode"].to_list())
        #keep specific variant classifications
        df_mutation_final = df_mutation[df_mutation["Variant_Classification"].isin(["Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"])]
        
        #parsing the 'HGVSp_Short' column 
        pd.options.mode.chained_assignment = None  # default='warn'
        #
        df_mutation_final['HGVSp_Short_'] = df_mutation_final['HGVSp_Short'] .str.replace(r'p.', '',regex=True)
        df_mutation_final = df_mutation_final[~df_mutation_final['HGVSp_Short_'].str.contains('_', na=False)] #remove mutations affecting more than one position
        len(set(df_mutation["Tumor_Sample_Barcode"].to_list())) #9716 tumors 
        
        df_mutation_final['HGVSp_Short_'] = df_mutation_final['HGVSp_Short_'].str.split('fs*').str[0]
        df_mutation_final['HGVSp_Short_'] = df_mutation_final['HGVSp_Short_'].str.split('ext*').str[0]
        
        #derive new columns for wild type, mutant residues and residue number
        df_mutation_final= df_mutation_final.assign(ResidueNumber = lambda x: x['HGVSp_Short_'].str.extract('(\d+)'))
        df_mutation_final['WildTypeResidue'] = df_mutation_final['HGVSp_Short_'].astype(str).str[0]
        df_mutation_final['MutantResidue'] = df_mutation_final['HGVSp_Short_'].astype(str).str[-1]
        df_mutation_final["WT_Residue"]=df_mutation_final['WildTypeResidue']+df_mutation_final["ResidueNumber"]
        
        #remove rows where wt or mutant residues are not known
        #remove rows with . and ?
        df_mutation_final_=df_mutation_final[df_mutation_final['WildTypeResidue']!="?"]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['WildTypeResidue']!="."]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['MutantResidue']!="?"]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['MutantResidue']!="."]
        df_mutation_final_=df_mutation_final_[df_mutation_final_["WT_Residue"]!="nan"]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['WildTypeResidue']!="nan"]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['MutantResidue']!="nan"]
        
        
        print(len(set(df_mutation_final_["Tumor_Sample_Barcode"].to_list()))) #9703 tumors with point mutations out of 9718 tumors
        
        df_mutation_final_["VAF"]=df_mutation_final_["t_alt_count"]/df_mutation_final_["t_depth"]
        ########write Genie cleaned mutation file as txt file
        #write to file
        df_mutation_final_.to_csv("./OOP_OUTPUT/OOP_Generated_TCGA_MC3callSet_Mutation_Cleaned_Parsed_File__new.txt",sep="\t",index=False,header=True)
        self._preprocessed_mutation_df = df_mutation_final_
        return self._preprocessed_mutation_df



path1=PATH_NAME
file_name1 = "MC3_TCGA_Mutations.txt"

df = TCGA_MAF(path=path1,file_name = file_name1) 
df.PreprocessingTcgaMutationFile()
df_mut = df.mutation_df       
# Tumors=len(df.GetTumorSamples()) #9718 tumors ending with -01
# df.LengthTumors()
