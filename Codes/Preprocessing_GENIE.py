#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 21:06:49 2022

@author: bengi
"""
import pandas as pd
import os


class GENIE_MAF:
    
    def __init__(self, path, file_name):
        self.path = path
        self.file_name = file_name
        
        self.file_path = os.path.join(self.path, self.file_name)

        self.name = "GENIE-Mutation Annotation File Preprocessing"
        #read a MAF file with the following columns listed in --SpecificColumns--
        SpecificColumns = ['Hugo_Symbol','Center','Variant_Classification', 'Variant_Type', "Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2", \
                         'Tumor_Sample_Barcode','Matched_Norm_Sample_Barcode', 'HGVSp_Short',\
                         "t_ref_count","t_alt_count","Codons","SWISSPROT","t_depth"] #select specific columns to read
        
            
        self.mutation_df = pd.read_csv( self.file_path, sep = "\t", usecols = SpecificColumns, comment="#")
    
    def name(self):
        return self.name
    
    def tumor_df(self, tumor_file = "data_clinical_sample.txt"):
        self._tumor_df = pd.read_csv(os.path.join(self.path, tumor_file), sep = "\t",comment = "#")
        return self._tumor_df
    
    def patient_df(self, patient_file = "data_clinical_patient.txt"):
        self._patient_df = pd.read_csv(os.path.join(self.path, patient_file), sep = "\t",comment = "#")
        return self._patient_df
        
    def all_patients(self):
        
        df_patient = self.patient_df()
        #returns {patient:{tumor:sample_type}} dictionary to check patients with mutltiple samples
        self._all_patient_list =  df_patient["PATIENT_ID"].to_list()
        return self._all_patient_list
    
    def all_tumors(self):
        df_tumor = self.tumor_df()
        #returns {patient:{tumor:sample_type}} dictionary to check patients with mutltiple samples
        self._all_tumor_list =  df_tumor["SAMPLE_ID"].to_list()
        return  self._all_tumor_list
    
    def PatientSampleDictionary(self):
 
        Patient_Tumor_SampType_Dict={ patient:{} for patient in self.all_patients() }
        df_tumor = self.tumor_df()
        for ind in df_tumor.index:
            patient = df_tumor.iloc[ind]['PATIENT_ID']
            tumor = df_tumor.iloc[ind]['SAMPLE_ID']
            sampType = df_tumor.iloc[ind]['SAMPLE_TYPE']
            if patient in Patient_Tumor_SampType_Dict.keys():
                Patient_Tumor_SampType_Dict[patient][tumor]=sampType
            else:
               print((ind,patient))
        
        self._Patient_Tumor_SampType_Dict = Patient_Tumor_SampType_Dict
        return self._Patient_Tumor_SampType_Dict
    
    def PatientsWithMultipleSamples(self):
        Patient_Tumor_SampType_Dict = self.PatientSampleDictionary()
        multiple_tumors = []
        counter=0
        for pat in Patient_Tumor_SampType_Dict.keys():
            if len(Patient_Tumor_SampType_Dict[pat])>=2:
                multiple_tumors.append((pat,len(Patient_Tumor_SampType_Dict[pat])))
                counter+=1
        
        self._multiple_tumors = multiple_tumors
        return self._multiple_tumors
    #match each patient with unique sample
    def PatientTumorSampleOneToOneMatch(self):
        
        c_multiple=0 #keep the number of patients with multiple tumor samples
        c_primary=0
        c_metastasis=0
        c_notspecified=0
        c_else=0
        c_one=0

        Patient_Tumor_SampType_Dict = self.PatientSampleDictionary()

        cleaned_Patient_Tumor_SampType_Dict={}
        for patient in Patient_Tumor_SampType_Dict.keys():
            if len(Patient_Tumor_SampType_Dict[patient])==1:
                c_one+=1
                cleaned_Patient_Tumor_SampType_Dict[patient]=Patient_Tumor_SampType_Dict[patient]
            elif len(Patient_Tumor_SampType_Dict[patient])>1:
                c_multiple+=1
                cleaned_Patient_Tumor_SampType_Dict[patient]={}
                for tumor, sampType in Patient_Tumor_SampType_Dict[patient].items():  
                    if "Primary" in Patient_Tumor_SampType_Dict[patient].values():
                        if sampType == "Primary":
                            c_primary+=1
                            cleaned_Patient_Tumor_SampType_Dict[patient][tumor]=Patient_Tumor_SampType_Dict[patient][tumor]
                            break
                    else:
                        if sampType == "Metastasis":
                            c_metastasis+=1
                            cleaned_Patient_Tumor_SampType_Dict[patient][tumor]=Patient_Tumor_SampType_Dict[patient][tumor]
                            break
                        elif (sampType == "Unspecified") or (sampType == "Not Collected"):
                            c_notspecified+=1
                            cleaned_Patient_Tumor_SampType_Dict[patient][tumor]=Patient_Tumor_SampType_Dict[patient][tumor]
                            break
                        else: #category: Not Applicable or Heme
                            cleaned_Patient_Tumor_SampType_Dict[patient][tumor]=Patient_Tumor_SampType_Dict[patient][tumor]
                            #print((patient,sampType))
                            c_else+=1
                            break
        self._cleaned_Patient_Tumor_SampType_Dict = cleaned_Patient_Tumor_SampType_Dict
        #write the number opatients with multiple mutations, how many of them 
        print("Distribution of sample types among multiple tumors and also the number of patients with one tumor sample")
        print("Multiple:{}, Primary:{},Metastasis:{},Not_spec:{},Else:{},One_tumor:{}".format(c_multiple,c_primary,c_metastasis,c_notspecified,c_else,c_one))
            #write patient, matching tumor and sample type to file
        ####### optional
        # PATH = "/Users/bengi/Desktop/CommsBioRevision_Final_9Dec/REVISION#2/DATA_CODES_new"
        # for patient in cleaned_Patient_Tumor_SampType_Dict.keys():
        #     for tumor in cleaned_Patient_Tumor_SampType_Dict[patient]:
        #         with open(PATH+"/OOP_OUTPUT/Genie_PAtient_Unique_Tumor_SampleType_Dict.txt","a") as outfile:
        #             outfile.write(patient+"\t"+tumor+"\t"+cleaned_Patient_Tumor_SampType_Dict[patient][tumor]+"\n")
        return self._cleaned_Patient_Tumor_SampType_Dict        
              
    #62,471 patient with exactly one matching tumor
    #2930 patient matches with more than one tumor, among these 2019 have primary tumor
                        #757 metastatic, 48 not specified, 106 not specified or heme

    #######
    
    #write to file
    #df_tumor_unique_.to_csv("GENIE_vol6_Unique_Tumor_Oncotree_SampleType__new.txt",sep="\t",index=False,header=True)

    ###65,401 unique tumors 
    #####Result of the above code: 2019 primary, 757 metastatic, 48 not specified, 106 Not Applicable or Heme _> total 2930

    #Preprocessing of the Mutation File "data_mutations_extended.txt"
    #mutation data file path
    # PATH = "/Users/bengi/Desktop/CommsBioRevision_Final_9Dec/REVISION#2/GENIE_DATA_VOL_6/data_mutations_extended.txt"

    #Read GENIE mutation data
 
       
    def PreprocessingGenieMutationFile(self) :   
        #from df_tumor, select the rows only with the tumor sample_id's from the cleaned_Patient_Tumor_SampType_Dict's inner dict keys
        cleaned_Patient_Tumor_SampType_Dict = self.PatientTumorSampleOneToOneMatch()
        unique_tumor_ids=[tumor for k, v in cleaned_Patient_Tumor_SampType_Dict.items()  for tumor in v.keys() ]
        df_tumor = self.tumor_df()
        df_tumor_unique = df_tumor[df_tumor["SAMPLE_ID"].isin(unique_tumor_ids)]
        df_tumor_unique_=df_tumor_unique[["PATIENT_ID", 'SAMPLE_ID','ONCOTREE_CODE','SAMPLE_TYPE']]

        df_mutation = self.mutation_df
        #after filtering t_depth there reamins 62081 tumors
        df_mutation1 = df_mutation[df_mutation["t_depth"]!=0] #will be used to calculate VAF
        df_mutation1.reset_index(inplace=True,drop=True)
        print(len(set(df_mutation1["Tumor_Sample_Barcode"].to_list())))
        #df_mutation.iloc[68]["HGVSp_Short"]   
        #keep non_null values for "HGVSp_Short" column

        df_mutation2 = df_mutation1[df_mutation1["HGVSp_Short"].notnull()==True] #61788 tumors with nonnull hgvsp-short

        df_mutation2.reset_index(inplace=True,drop=True)
        print(len(set(df_mutation2["Tumor_Sample_Barcode"].to_list()))) 

        ####
        #take the rows with unique tumor samples
        df_mutation_=df_mutation2[df_mutation2['Tumor_Sample_Barcode'].isin(unique_tumor_ids)]
        print(len(set(df_mutation_["Tumor_Sample_Barcode"].to_list())))#58769 unique tumor-patient remain

        #take specific variant classifications
        df_mutation_final = df_mutation_[df_mutation_["Variant_Classification"].isin(["Frame_Shift_Del","Frame_Shift_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation"])]
        print(len(set(df_mutation_final["Tumor_Sample_Barcode"].to_list()))) #57926 unique tumor-patient remain with the above variant classifications
        
        ############
        #parsing the 'HGVSp_Short' column 
        pd.options.mode.chained_assignment = None  # default='warn'
        #
        df_mutation_final['HGVSp_Short_'] = df_mutation_final['HGVSp_Short'].str.replace(r'p.', '', regex=True)

        df_mutation_final = df_mutation_final[~df_mutation_final['HGVSp_Short_'].str.contains('_', na=False)] #remove mutations affecting more than one position


        df_mutation_final['HGVSp_Short_'] = df_mutation_final['HGVSp_Short_'].str.split('fs*').str[0]
        df_mutation_final['HGVSp_Short_'] = df_mutation_final['HGVSp_Short_'].str.split('ext*').str[0]

        #derive new columns for wild type, mutant residues and residue number
        df_mutation_final= df_mutation_final.assign(ResidueNumber = lambda x: x['HGVSp_Short_'].str.extract('(\d+)'))

        df_mutation_final['WildTypeResidue'] = df_mutation_final['HGVSp_Short_'].astype(str).str[0]
        df_mutation_final['MutantResidue'] = df_mutation_final['HGVSp_Short_'].astype(str).str[-1]
        df_mutation_final["WT_Residue"]=df_mutation_final['WildTypeResidue']+df_mutation_final["ResidueNumber"]

        print(len(set(df_mutation_final["Tumor_Sample_Barcode"].to_list()))) #57,921 tumors with point mutations (out of 65401 unique tumors)

        #x=df_mutation_final[df_mutation_final["Hugo_Symbol"]=="B2M"]
        #x.reset_index(inplace=True,drop=True)
        #x.iloc[49]["HGVSp_Short"]
        #remove rows where wt or mutant residues are not known
        #remove rows with . and ?
        df_mutation_final_=df_mutation_final[df_mutation_final['WildTypeResidue']!="?"]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['WildTypeResidue']!="."]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['MutantResidue']!="?"]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['MutantResidue']!="."]

        df_mutation_final_=df_mutation_final_[df_mutation_final_["WT_Residue"]!="nan"]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['WildTypeResidue']!="nan"]
        df_mutation_final_=df_mutation_final_[df_mutation_final_['MutantResidue']!="nan"]

        print(len(set(df_mutation_final_["Tumor_Sample_Barcode"].to_list()))) #57,921 tumors with point mutations (out of 65401 unique tumors)

        df_mutation_final_["VAF"]=df_mutation_final_["t_alt_count"]/df_mutation_final_["t_depth"]
        ########write Genie cleaned mutation file as txt file
  
        self._preprocessed_mutation_df = df_mutation_final_
        
        return self._preprocessed_mutation_df
    
    def return_GENIE_patient_tumor_df(self):
        #from df_tumor, select the rows only with the tumor sample_id's from the cleaned_Patient_Tumor_SampType_Dict's inner dict keys
        cleaned_Patient_Tumor_SampType_Dict = self.PatientTumorSampleOneToOneMatch()
        unique_tumor_ids=[tumor for k, v in cleaned_Patient_Tumor_SampType_Dict.items()  for tumor in v.keys() ]
        df_tumor = self.tumor_df()
        df_tumor_unique = df_tumor[df_tumor["SAMPLE_ID"].isin(unique_tumor_ids)]
        df_tumor_unique_=df_tumor_unique[["PATIENT_ID", 'SAMPLE_ID','ONCOTREE_CODE','SAMPLE_TYPE']]
        
        self._GENIE_patient_tumor_df = df_tumor_unique_
        return self._GENIE_patient_tumor_df
        #df_tumor_unique_.to_csv("OOP_Generated_GENIE_vol6_CleanedPAtient_Tumor_Info.txt",sep="\t",index=False,header=True)   
        #pickle.dump(df_tumor_unique_, open("Preprocessed_GENIE_Patient_Sample_Info_df.p", "wb")) 


        
        # def return_Genie_DF(self):
        #     df_genie_preprocessed = self.PreprocessingGenieMutationFile()
        #     self.Genie_DF = df_genie_preprocessed
        #     return self.Genie_DF
        #     #returns {patient:{tumor:sample_type}} dictionary to check patients with mutltiple samples
        #     #self._preprocessed_genie = df_tumor
        #     #write to file
        #     #PATH = "/Users/bengi/Desktop/CommsBioRevision_Final_9Dec/REVISION#2/DATA_CODES_new"
        #     #OPTIONAL 
        #     """change file names   """
        #     #pickle.dump(df_genie_preprocessed, open("Preprocessed_GENIE_df.p", "wb")) 
        #     #df_genie_preprocessed.to_csv("OOP_Generated_GENIE_vol6_Mutation_Cleaned_Parsed_File__new.txt",sep="\t",index=False,header=True)
        
        
        
 
    

 