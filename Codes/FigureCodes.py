#!/usr/bin/env python3
# -*- coding: utf-8 -*-


"""Codes for reproducing main figures, the input files for each panel 
can be found in 10.6084/m9.figshare.21788192.      """


"""Figure 1A: Windrose Plot"""


# import pickle
# pickle.dump(Tissue_Tumor_Dict ,open("GENIE_TCGA_ClinicalVAF>12_5_Tissue_Tumor_Dictionary.p", "wb"))
# ##########

PATH = "./Data"
def Figure1A_Windrose_Plot():
    #Load mutation_tumor, gene_tumor dictionaries
    import pickle, pandas as pd
    Gene_Mutation_Tumor_GEQ3_Dict = pickle.load(open(PATH + "Gene_Mutation_Tumor_GEQ3_Dictionary.p", "rb"))
    Gene_Mutation_Tumor_LEQ2_Dict = pickle.load(open(PATH +"Gene_Mutation_Tumor_LEQ2_Dictionary.p", "rb"))
    Gene_All_Mutants = pickle.load(open(PATH +"Gene_All_Mutants_Dictionary.p", "rb"))
    Tissue_Tumor_Dict = pickle.load(open(PATH +"GENIE_TCGA_ClinicalVAF>12_5_Tissue_Tumor_Dictionary.p", "rb"))
    #load filtered double mutation tab delimited file
    
    df_doubles = pd.read_excel(PATH +"Supplementary_Table_1", sep = "\t")
    #194 doubles, 54 genes
    Genes = set(df_doubles["Gene"].to_list()) #54 genescarries at least one doublet observed on 3 tumors
    sorted(list(Genes))
    print("There are {} genes bearing at least one double mutation".format(len(Genes)))
    
    #Gene_
    Doublet_List = []
    
    for ind in df_doubles.index:
        Doublet_List.append((df_doubles.iloc[ind]["Gene"]+"_"+df_doubles.iloc[ind]["Mut1"],df_doubles.iloc[ind]["Gene"]+"_"+df_doubles.iloc[ind]["Mut2"]))

    #Tumors caarrying these 194 doublets
    Doublet_Tumor_Dict = { } #(gene_mut1,gene_mut2) are keys

    for double in Doublet_List:
        Doublet_Tumor_Dict[double] = Gene_Mutation_Tumor_GEQ3_Dict[double[0]].intersection(Gene_Mutation_Tumor_GEQ3_Dict[double[1]])
#############3
    #################Number of all double mutant tumors : 
    All_Double_Mutant_Tumors = set()

    for d in Doublet_Tumor_Dict.keys():
        All_Double_Mutant_Tumors = All_Double_Mutant_Tumors.union(Doublet_Tumor_Dict[d])
    print("Number of all double mutant tumors is {}".format(len(All_Double_Mutant_Tumors)))   
    #############3
    ###########For each tissue find the number of double mutant tumors
    Tissue_Double_Mutant_Tumor_Count = {}  

    for tissue in Tissue_Tumor_Dict.keys() :
        Tissue_Double_Mutant_Tumor_Count[tissue] = len(Tissue_Tumor_Dict[tissue].intersection(All_Double_Mutant_Tumors))
            
        ####
        ####
    #### Windrose Plot input file
    """
    with open("./FIGURES/Figure1A_Windrose_Plot_Input.txt" , "a") as outfile:
        outfile.write("Tissue"+"\t"+"CountType"+"\t"+"Count"+"\n")
        
    for tissue in Tissue_Double_Mutant_Tumor_Count.keys():
        Length = str(len(Tissue_Tumor_Dict[tissue]))
        Tissue = tissue + " (N = {})".format(Length)
        DoubleMut = str(Tissue_Double_Mutant_Tumor_Count[tissue])
        NotDoubleMut = str(len(Tissue_Tumor_Dict[tissue])-Tissue_Double_Mutant_Tumor_Count[tissue])
        with open("./FIGURES/Figure1A_Windrose_Plot_Input.txt" , "a") as outfile:
            outfile.write(Tissue+"\t"+"DoubleMut"+"\t"+str(DoubleMut)+"\n"+ \
                          Tissue+"\t"+"NotDoubleMut"+"\t"+str(NotDoubleMut)+"\n")
    """
    ###############
    #Plotting Windrose Plot

    !pip install windrose
    !pip install plotly
    !pip install kaleido
    from windrose import plot_windrose 
    import numpy as np
    np.log10(0)
    np.log(1256)
    np.log(5073)
    from matplotlib import cm
    import plotly.express as px
    import matplotlib.pyplot as plt
    from plotly import graph_objects
    import math
    !pip install -U orca
    !pip install orca
    !pip install plotly
    from kaleido.scopes.plotly import PlotlyScope
    import plotly.graph_objects as go
    import plotly.express as px
    #df = px.data.wind()
    df2 = pd.read_csv("./FIGURES/Figure1A_Windrose_Plot_Input.txt",sep="\t")
    #df2.columns=["direction","#Patients","frequency"]
    df2 = df2.sort_values(["Tissue"], ascending=True) \
        .groupby(["Tissue"], sort=False) \
        .apply(lambda x: x.sort_values(['CountType'], ascending=True)) \
        .reset_index(drop=True)


    scope = PlotlyScope()
    fig = px.bar_polar(df2, r="Count", theta="Tissue",
                       color="CountType",log_r=True)

    fig.update_layout(
        title='Sample Size Distributions with and without Double Mutations (LogScaleAxis)\n \
            Among only 62566 tumors with tissue info out of 62567 that have point mut and VAF>12.5',
        font_size=6,
        legend_font_size=12,legend_font_family ="Arial")

    #fig.update_layout(
    #    template=None,
    #    polar = dict(
    #        radialaxis = dict( showticklabels=True)))
    fig.show()
    fig.update_layout(width=600,height=600)
    fig.write_image("./FIGURES/Figure1A_TissueWindRose_LogScaleAxis.pdf")
    fig.write_image("./FIGURES/Figure1A_TissueWindRose_LogScaleAxis.svg")
    fig.show()
    #


#figure1A = Figure1A_Windrose_Plot()       



###Figure 1B, Cancer specifict Scatter Plot

def Figure1B_CancerSpecificity_Scatter_Plot():
    #Load mutation_tumor, gene_tumor dictionaries
    import pickle, pandas as pd
    Gene_Mutation_Tumor_GEQ3_Dict = pickle.load(open(PATH +"Gene_Mutation_Tumor_GEQ3_Dictionary.p", "rb"))
    Gene_Mutation_Tumor_LEQ2_Dict = pickle.load(open(PATH +"Gene_Mutation_Tumor_LEQ2_Dictionary.p", "rb"))
    Gene_AllMutants = pickle.load(open(PATH +"Gene_All_Mutants_Dictionary.p", "rb"))
    Tissue_Tumor_Dict = pickle.load(open(PATH +"GENIE_TCGA_ClinicalVAF>12_5_Tissue_Tumor_Dictionary.p", "rb"))
    #load filter
    
    #load filtered double mutation tab delimited file
    
    df_doubles = pd.read_excel(PATH +"Supplementary_Table_1", sep = "\t")
    #194 doubles, 54 genes
    Genes = set(df_doubles["Gene"].to_list()) #59 genescarries at least one doublet observed on 3 tumors
    sorted(list(Genes))
    print("There are {} genes bearing at least one double mutation".format(len(Genes)))
    
    #Gene_
    Doublet_List = []
    
    for ind in df_doubles.index:
        Doublet_List.append((df_doubles.iloc[ind]["Gene"]+"_"+df_doubles.iloc[ind]["Mut1"],df_doubles.iloc[ind]["Gene"]+"_"+df_doubles.iloc[ind]["Mut2"]))

    #Tumors caarrying these 194 doublets
    Doublet_Tumor_Dict = { } #(gene_mut1,gene_mut2) are keys

    for double in Doublet_List:
        Doublet_Tumor_Dict[double] = Gene_Mutation_Tumor_GEQ3_Dict[double[0]].intersection(Gene_Mutation_Tumor_GEQ3_Dict[double[1]])
#############3
    #################Number of all double mutant tumors : 
    All_Double_Mutant_Tumors = set()

    for d in Doublet_Tumor_Dict.keys():
        All_Double_Mutant_Tumors = All_Double_Mutant_Tumors.union(Doublet_Tumor_Dict[d])
    print("Number of all double mutant tumors is {}".format(len(All_Double_Mutant_Tumors)))   
    #############3
    ###########For each tissue find the number of double mutant tumors
    Tissue_Double_Mutant_Tumor_Count = {}  

    for tissue in Tissue_Tumor_Dict.keys() :
        Tissue_Double_Mutant_Tumor_Count[tissue] = len(Tissue_Tumor_Dict[tissue].intersection(All_Double_Mutant_Tumors))
            
        ####****
    #total number of double components 
    DoubleComponents=set()  #??? double mut components
    for double in Doublet_Tumor_Dict.keys():
        DoubleComponents.add(double[0])
        DoubleComponents.add(double[1])
    #############3
    print("There are {} double mutation components of {} doublets".format(len(DoubleComponents),len(Doublet_Tumor_Dict)))
    #######################
    #for the  genes, get the list of double mutant tumors
    Gene_DoubleMutants_Dict = { g:set() for g in Genes}
    
    #double  mutant all tumors, will be used for filtering tissues
    for d in Doublet_Tumor_Dict.keys():
        d1 = d[0]
        d2 = d[1]
        gen = d1.split("_")[0]
        Gene_DoubleMutants_Dict[gen] = Gene_DoubleMutants_Dict[gen].union(Doublet_Tumor_Dict[d])

    ##############take the percentage among the single mutations that are 
    #components of double mutations
   
    #divide the number of double mutantx to all gene mutants to get prevalance value
    Prevalance_Dict = {}

    for gene in Genes:
        Prevalance_Dict[gene] = len(Gene_DoubleMutants_Dict[gene]) / len(Gene_AllMutants[gene])
        

    ################Select tissues that carries at least 3 tumors with double muatations
    TissueFiltering =[]
    for tissue in Tissue_Tumor_Dict.keys():
        if len(All_Double_Mutant_Tumors.intersection(Tissue_Tumor_Dict[tissue])) >= 0:
            TissueFiltering.append(tissue)
            
    #For each gene find the intersection with each tissue for double mutant and all mutant tumors
    Gene_DoubleMutants_Tissue_Dict = {g:{} for g in Genes}

    for gene in Genes:
        for tissue in TissueFiltering:
            Gene_DoubleMutants_Tissue_Dict[gene][tissue] = len(Gene_DoubleMutants_Dict[gene].intersection(Tissue_Tumor_Dict[tissue]))
    ####
    #count the number of tissues that have >0 double mutants
    TDouble_dict = {gene : 0 for gene in Genes}
    for gene in Gene_DoubleMutants_Tissue_Dict.keys():
        for tissue in Gene_DoubleMutants_Tissue_Dict[gene].keys():
            if Gene_DoubleMutants_Tissue_Dict[gene][tissue] >0 :
                TDouble_dict[gene]+=1
                
    #########3
    #Select Tissues where at least 3 tumors carries a double mutation

    ##############
    #For each gene find the intersection with each tissue for double mutant and all mutant tumors
    Gene_AllMutants_Tissue_Dict = {g:{} for g in Genes}

    for gene in Genes:
        for tissue in TissueFiltering:
            Gene_AllMutants_Tissue_Dict[gene][tissue] = len(Gene_AllMutants[gene].intersection(Tissue_Tumor_Dict[tissue]))

    ################
    #count the number of tissues that have >0 single mutants
    TSingle_dict = {gene : 0 for gene in Genes}
    for gene in Genes:
        for tissue in Gene_AllMutants_Tissue_Dict[gene].keys():
            if Gene_AllMutants_Tissue_Dict[gene][tissue] >0 :
                TSingle_dict[gene]+=1
                
    #divide number of double mutant tissues to number of single mutant tissues            
    TDouble_TSingle_Ratio = {}            

    for gene in Genes:
        TDouble_TSingle_Ratio[gene] = TDouble_dict[gene] / TSingle_dict[gene]
        
        
    ############
    #Be careful about the thresholds
    with open("./FIGURES/Figure1B_CancerSpecificty_Input_OG_TSG.txt","a") as outfile:
        outfile.write("Gene"+"\t"+"TDouble_TSingle"+"\t"+"Prevalance"+"\t"+"TSingle"+"\t"+"Tdouble"+"\t"+"SingleMutant#"+"\t"+"DoubleMutant#"+"\n")              

    for gene in TDouble_TSingle_Ratio.keys():
        with open("./FIGURES/Figure1B_CancerSpecificty_Input_OG_TSG.txt","a") as outfile:
            outfile.write(gene+"\t"+str(TDouble_TSingle_Ratio[gene])+"\t"+str(Prevalance_Dict[gene])+"\t"+str(TSingle_dict[gene])+"\t"+str(TDouble_dict[gene])+"\t"+str(len(Gene_AllMutants[gene]))+"\t"+str(len(Gene_DoubleMutants_Dict[gene]))+"\n")              



    #######3
    #scatter plot from input file

    import pandas as pd
    df = pd.read_csv("./FIGURES/Figure1B_CancerSpecificty_Input_OG_TSG.txt", sep="\t")
    df.columns
    df=df[["Gene","TDouble_TSingle","Prevalance","TSingle"]]    
    df.sort_values(by=["TDouble_TSingle"],inplace=True) 
    df["Color"]  = [None]  * len(df)
    for ind in df.index:
        if df.loc[ind,"TDouble_TSingle"]<=0.20:
            df.at[ind,"Color"] = "Red"
        elif (df.loc[ind,"TDouble_TSingle"]<=0.40) and df.loc[ind,"TDouble_TSingle"]>0.20:
            df.at[ind,"Color"] = "Pink"  
        elif (df.loc[ind,"TDouble_TSingle"]>0.40):
            df.at[ind,"Color"] = "Skyblue"  


    # basic plot
    import seaborn as sns
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors as mcolors
    import numpy as np
    plt.figure(figsize=(9,6))
    #ax = fig.add_subplot(111)  
    ax = plt.gca()
    ax.set_facecolor('xkcd:white')
    palette={"Red":"red","Pink":"pink", "Skyblue":"skyblue" }
    p1=sns.scatterplot(data=df, x="TDouble_TSingle", y="Prevalance", marker="o", vmin=0,vmax=0.6, size="TSingle",sizes=(10, 300),hue="Color",palette=palette)

    # add annotations one by one with a loop
    for line in range(0,df.shape[0]):
         p1.text(df.TDouble_TSingle[line]-0.01, df.Prevalance[line]+0.002, df.Gene[line], horizontalalignment='left', size='smaller', color='gray', weight='bold')
    plt.xticks(np.arange(0,0.7, 0.1) )
    plt.yticks(np.arange(0, 0.6, 0.05)) 
    plt.xlabel("TDouble/TSingle")
    plt.ylabel("Prevalence of double mutations")
    plt.ylim(ymax = 0.20, ymin = 0)
    plt.xticks(fontname = "Arial",fontsize=13)            
    plt.yticks(fontname = "Arial",fontsize=13)  # This argument will change the font.

    plt.savefig("./FIGURES/Figure1B_CancerSpecificty_OG_TSG.svg",bbox_inches='tight',dpi=12000)

    plt.savefig("./FIGURES/Figure1B_CancerSpecificty_OG_TSG.pdf",bbox_inches='tight',dpi=12000)


    plt.show()


    ####write gene names to file

    df_red = df[df["Color"]=="Red"]
    df_red = df_red[["Gene"]]

    df_red.sort_values(by=["Gene"],inplace=True)
    df_red.to_csv("./FIGURES/Fig1B_Legend_RedGenes_OG_TSG.txt",sep="\t",index=False)


#figure1B = Figure1B_CancerSpecificity_Scatter_Plot()

# import pandas as pd
# df_doubles = pd.read_csv("Rev#3_Significant_doubles_VAF_Nonsense_filtered.txt", sep = "\t")
# #194 doubles, 54 genes
# Genes = set(df_doubles["Gene"].to_list()) #59 genescarries at least one doublet observed on 3 tumors
# sorted(list(Genes))

######Figure 1C, stacked bar plot

def Figure1C_StackedBarPlot_OG_TSG():
    
    #Stacked bar plot showing the number of DD,Dd,dd type doublets for OG and TSG's

    import pandas as pd
    import numpy as np

    df_new1 = pd.read_excel(PATH +"Figure1C_input", sep = "\t")

    df_new2 = df_new1[ df_new1['OG_TSG'].isin(["OG","TSG"])]

    df_new2 = df_new2[['Gene', 'Mut1', 'Mut2','OG_TSG','OncogenicDriver_Mut1', 'OncogenicDriver_Mut2']]

    df_new2["Constituents"] = df_new2['OncogenicDriver_Mut1'] + "+" + df_new2['OncogenicDriver_Mut2']

    df_og = df_new2[ df_new2["OG_TSG"] == "OG" ] #80 doublets
    og = len(set(df_og["Gene"].to_list())) #14 oncogenes
    set(df_og["Gene"].to_list())

    df_og_dd = df_og[ df_og["Constituents"] == "driver[d]+driver[d]" ]

    df_og_DD = df_og[ df_og["Constituents"] == "KnownDriver[D]+KnownDriver[D]" ]
    og_Dd = len(df_og)-(len(df_og_dd)+len(df_og_DD))

    ##############

    df_tsg = df_new2[ df_new2["OG_TSG"] == "TSG" ] #170 doublets
    tsg = len(set(df_tsg["Gene"].to_list())) # 30 tumor suppressor genes

    df_tsg_dd = df_tsg[ df_tsg["Constituents"] == "driver[d]+driver[d]" ]
    df_tsg_DD = df_tsg[ df_tsg["Constituents"] == "KnownDriver[D]+KnownDriver[D]" ]
    tsg_Dd = len(df_tsg)-(len(df_tsg_dd)+len(df_tsg_DD))

    DataDict = {"OG":{"dd":len(df_og_dd) , "DD":len(df_og_DD) , "Dd":og_Dd } ,\
                "TSG":{"dd":len(df_tsg_dd) , "DD": len(df_tsg_DD), "Dd": tsg_Dd}}
        
    ####
    #Draw the stacked bar chart
    import pandas as pd
    import matplotlib.pyplot as plt

    data={"OG":[(DataDict["OG"]["dd"]/len(df_og))*100, (DataDict["OG"]["Dd"]/len(df_og))*100, (DataDict["OG"]["DD"]/len(df_og))*100],\
        "TSG":[(DataDict["TSG"]["dd"]/len(df_tsg))*100 , (DataDict["TSG"]["Dd"]/len(df_tsg))*100 , (DataDict["TSG"]["DD"]/len(df_tsg))*100],\
          "Type":["dd","Dd","DD"]}
    d=pd.DataFrame.from_dict(data)
    d=d.set_index('Type',drop=True)
    d1=d.T
    import scipy.stats as stats
    # from scipy.stats import mannwhitneyu
    # U1, pvalue = mannwhitneyu(data["TSG"], data["OG"], method="auto")
    #statistic,pvalue=stats.ttest_ind()
    #ttest and mann whitney u test not returning significant p values for percentage or counts
    #I compared the number of dd and DD between og and tsg's
    vect1 = [DataDict["OG"]["dd"],DataDict["OG"]["DD"]] 
    vect2 = [DataDict["TSG"]["dd"],DataDict["TSG"]["DD"]]
    oddsratio, pvalue = stats.fisher_exact([vect1, vect2])


    ax = d1[["dd","Dd","DD"]].plot(kind="bar", stacked=True,width=0.5,figsize=(2,4))
    ax.legend(bbox_to_anchor=(2,0.3), loc="lower right")
    ax.set_xlabel("Percentages of Double Mutations on Same TSG/OG,\n FisherExact test, p_val= {:.12f})\n \
                  compare dd and DD numbers".format(pvalue))
    plt.xticks(fontname = "Arial",fontsize=14)            
    plt.yticks(fontname = "Arial",fontsize=14)  # This argument will change the font.

    fig = ax.get_figure()
    fig.savefig("./FIGURES/Figure1C_TSG_OG_Percentage_Stacked.pdf",bbox_inches='tight', dpi=12000)
    fig.savefig("./FIGURES/Figure1C_TSG_OG_Percentage_Stacked.svg",bbox_inches='tight', dpi=12000)
    pvalue


    # #write these to the f,gure
    # # OG:{'dd': 6, 'DD': 54, 'Dd': 20}
    # # TSG:{'dd': 47, 'DD': 68, 'Dd': 55}

    # {'TSG': [27.647058823529413, 32.35294117647059, 40.0], 'OG': [7.5, 25.0, 67.5], 'Type': ['dd', 'Dd', 'DD']}
    # ['dd', 'Dd', 'DD']

#figure1C = Figure1C_StackedBarPlot_OG_TSG()

#Figure 1D, Passenger Load

def Figure1D_BoxPlot_OG_TSG_PassengerLoad():
  
    ####Figure 1D Box Plot for og-tsg passenger loads
    #label known drivers and drivers - the rest will be passenger mutation
    #compare the passenger loads of tumors with at least one og or tsg double mutation

#Stacked bar plot showing the number of DD,Dd,dd type doublets for OG and TSG's

    import pandas as pd, numpy as np, pickle 
    
    from matplotlib import pyplot as plt, seaborn as sns
    
    #Read in data & create total column

    
    df_new1 = pd.read_csv(PATH+"Figure1D_PassengerLoad.txt", sep = "\t")
    
    df_new2 = df_new1[ df_new1['OG_TSG'].isin(["OG","TSG"])]
    df_new2.reset_index(inplace=True)
    df_new2 = df_new2[['Gene', 'Mut1', 'Mut2','OG_TSG','OncogenicDriver_Mut1', 'OncogenicDriver_Mut2']]

    ##### D:known driver d:driver
    Dd_Dict = {} 
    OG_Doublets = [] # 80 doublets-14 og
    TSG_Doublets = [] #170 doublets-30 tsg
    og = set()
    tsg=set()
    for ind in df_new2.index:
        mut1 = df_new2.iloc[ind]["Gene"] + "_" + df_new2.iloc[ind]["Mut1"]
        mut2 = df_new2.iloc[ind]["Gene"] + "_" + df_new2.iloc[ind]["Mut2"]
        type1 = df_new2.iloc[ind]["OncogenicDriver_Mut1"]
        type2 = df_new2.iloc[ind]["OncogenicDriver_Mut2"]
        Dd_Dict[mut1] = type1
        Dd_Dict[mut2] = type2
        if df_new2.iloc[ind]["OG_TSG"] == "OG":
            OG_Doublets.append((mut1,mut2))
            og.add(df_new2.iloc[ind]["Gene"])
        elif df_new2.iloc[ind]["OG_TSG"] == "TSG":
            TSG_Doublets.append((mut1,mut2))
            tsg.add(df_new2.iloc[ind]["Gene"])
           
            
           #251 doublet constituents
    d_count=0 #134 driver[d]
    D_count=0 #117 KnownDriver[D]
    for mut in Dd_Dict.keys():
        if Dd_Dict[mut]=="driver[d]" :
            d_count+=1
        else:
            D_count+=1
            
            
    ###get the list of all drivers
    
    
    All_Drivers = pickle.load(open("OncogenicMutations_Dictionary.p","rb"))        
    
    len(set(All_Drivers))  #4196 known driver + driver counts   
    ########
    #get the list of og|tsg-double mutant tumors

    OG_DoubleMutant = set() # 586 tumors
    
    TSG_DoubleMutant = set() # 696 tumors
    Gene_Mutation_Tumor_GEQ3_Dict = pickle.load(open("Gene_Mutation_Tumor_GEQ3_Dictionary.p", "rb"))
    ######3
    
    #Gene_
    Doublet_List = []
    
    for ind in df_new2.index:
        Doublet_List.append((df_new2.iloc[ind]["Gene"]+"_"+df_new2.iloc[ind]["Mut1"],df_new2.iloc[ind]["Gene"]+"_"+df_new2.iloc[ind]["Mut2"]))

    #Tumors caarrying these 194 doublets
    Doublet_Tumor_Dict = { } #(gene_mut1,gene_mut2) are keys

    for double in Doublet_List:
        Doublet_Tumor_Dict[double] = Gene_Mutation_Tumor_GEQ3_Dict[double[0]].intersection(Gene_Mutation_Tumor_GEQ3_Dict[double[1]])
#############3

    for double in Doublet_Tumor_Dict.keys():
            m1 = double[0]
            m2 = double[1]
            if ((m1,m2) in OG_Doublets) or ((m2,m1) in OG_Doublets) :
                OG_DoubleMutant = OG_DoubleMutant.union(Doublet_Tumor_Dict[double])
            elif ((m1,m2) in TSG_Doublets) or ((m2,m1) in TSG_Doublets) :
                TSG_DoubleMutant = TSG_DoubleMutant.union(Doublet_Tumor_Dict[double])
    Intersect = OG_DoubleMutant.intersection(TSG_DoubleMutant)  #29 tumors common in both  
    Intersect    
    
    print("OG double-mutant tumors: {}, TSG double-mutant tumors: {}, have OG-TSG {} ".format(len(OG_DoubleMutant), len(TSG_DoubleMutant), len(Intersect )))
        
    ####################
    ###################
    #get the list of all mutations in Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt file
    #for the double mutant tumors
    OG_Mutants_Passenger_Mutation_Info = {tumor1:set() for tumor1 in OG_DoubleMutant} #remove 29 tumors in the intersection
    TSG_Mutants_Passenger_Mutation_Info = {tumor2:set() for tumor2 in TSG_DoubleMutant}
    with open("Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt","r") as infile:
        for line in infile:
            splitted = line.rstrip("\n").split("\t")
            mut = splitted[0]+"_"+splitted[-2]
            if splitted[2] in OG_Mutants_Passenger_Mutation_Info.keys() :
                if mut not in set(All_Drivers):
                    OG_Mutants_Passenger_Mutation_Info[splitted[2]].add(mut)
    with open("Merged_Genie_Tcga_VAFgreater12point5Percent__new.txt","r") as infile:
        for line in infile:
            splitted = line.rstrip("\n").split("\t")
            mut = splitted[0]+"_"+splitted[-2]
            if splitted[2] in TSG_Mutants_Passenger_Mutation_Info.keys() : 
                if mut not in set(All_Drivers):
                    TSG_Mutants_Passenger_Mutation_Info[splitted[2]].add(mut)
                    
    # og_genie_passenger={}
    # for tumor in OG_Mutants_Passenger_Mutation_Info.keys():
    #     if tumor.startswith("GENIE"):
    #         og_genie_passenger[tumor] = OG_Mutants_Passenger_Mutation_Info[tumor]
            
            
    # tsg_genie_passenger={}
    # for tumor in TSG_Mutants_Passenger_Mutation_Info.keys():
    #     if tumor.startswith("GENIE"):
    #         tsg_genie_passenger[tumor] = TSG_Mutants_Passenger_Mutation_Info[tumor] 
    #Write passengercounts to file
    #Write passengercounts to file
    #then draw the figure
    #tumor-og/tsg info-passenger counts
    for tumor in OG_Mutants_Passenger_Mutation_Info.keys():
        tag = "Oncogenes"
        count = str(len(OG_Mutants_Passenger_Mutation_Info[tumor]))
        with open("./FIGURES/Figure1D_PassengerLoad_Input.txt","a") as outfile:
            outfile.write(tumor+"\t"+tag+"\t"+count+"\n")
    for tumor in TSG_Mutants_Passenger_Mutation_Info.keys():
        tag = "Tumor Suppressors"
        count = str(len(TSG_Mutants_Passenger_Mutation_Info[tumor]))
        with open("./FIGURES/Figure1D_PassengerLoad_Input.txt","a") as outfile:
            outfile.write(tumor+"\t"+tag+"\t"+count+"\n")

    # len(TSG_Mutants_Passenger_Mutation_Info["TCGA-AP-A056-01"])
    # len(OG_Mutants_Passenger_Mutation_Info["TCGA-AP-A056-01"])


    '''
    #tumor count per cancer
    #patient counts per cancer
    #mutation counts per cancer: passenger + major + minor olarak stacked bar plot yap   '''

   

    df3 = pd.read_csv("./FIGURES/Figure1D_PassengerLoad_Input.txt",sep="\t",header=None) 
    
    
    df_og = df3[df3[1] == "Oncogenes"]
    df_tsg = df3[df3[1] == "Tumor Suppressors"]
    from scipy.stats import mannwhitneyu
    U1, p = mannwhitneyu(df_og[2].to_list(), df_tsg[2].to_list(),alternative='two-sided')
    p  
    
    sns.color_palette("pastel")  
    ax = sns.boxplot(x=1, y=2, data=df3)
    ax.set_yscale("log")
    ax.set_ylabel("Passenger Load")
    ax.set_xlabel("Passenger Load Comparison p_val= {} ".format(p))
    ax.figure.savefig("./FIGURES/Figure1D_OG_TSG_doubleMutants_PassengerLoad.svg",fontsize=12,family='Arial',fig_size=(3,4))
    ax.figure.savefig("./FIGURES/Figure1D_OG_TSG_doubleMutants_PassengerLoad.pdf",fontsize=12,family='Arial',fig_size=(3,4))
    ##################
    """igure1D = Figure1D_BoxPlot_OG_TSG_PassengerLoad()
    OG double-mutant tumors: 578, TSG double-mutant tumors: 371, have OG-TSG 14 """
  
                
#figure1D = Figure1D_BoxPlot_OG_TSG_PassengerLoad()

########

def Figure1E_Latent_Driver_TumorCount_BoxStripPlot():
    

    import pandas as pd, numpy as np

    df_new1 = pd.read_csv("Rev#3_Significant_doubles_VAF_Nonsense_filtered.txt", sep = "\t")
    
    df_new2 = df_new1[ df_new1['OG_TSG'].isin(["OG","TSG"])]
    df_new2.reset_index(inplace=True)

    ########

    ##### D:known driver d:driver
    Dd_Dict = {} 
    Count_Dict ={}

    for ind in df_new2.index:
        mut1 = df_new2.iloc[ind]["Gene"] + "_" + df_new2.iloc[ind]["Mut1"]
        mut2 = df_new2.iloc[ind]["Gene"] + "_" + df_new2.iloc[ind]["Mut2"]
        count1 = df_new2.iloc[ind]["DoubleMut#"] + df_new2.iloc[ind]["OnlyMut1"]
        count2 = df_new2.iloc[ind]["DoubleMut#"] + df_new2.iloc[ind]["OnlyMut2"]
        type1 = df_new2.iloc[ind]["OncogenicDriver_Mut1"]
        type2 = df_new2.iloc[ind]["OncogenicDriver_Mut2"]
        Dd_Dict[mut1] = type1
        Dd_Dict[mut2] = type2
        Count_Dict[mut1] = count1
        Count_Dict[mut2] = count2
        
    len(set(Dd_Dict.keys())) # whole set 295 known driver and driver constituents
          #251 doublet constituents
    d_count=0 #178 driver[d]
    D_count=0 #117 KnownDriver[D]
    for mut in Dd_Dict.keys():
        if Dd_Dict[mut]=="driver[d]" :
            d_count+=1
        else:
            D_count+=1
            
    #count the number of the constituents on tsg and og 
    ###get the list of all drivers
    #write mtation, mutant tumors, known driver-driver info as input file
    for mut in Count_Dict.keys():
        with open("./FIGURES/Figure1E_MajorMinorPatientCountBoxStripPlotInput.txt","a") as outfile:
            outfile.write(mut+"\t"+str(Count_Dict[mut])+"\t"+Dd_Dict[mut]+"\n")

    #####read the input file
    # check p value with mann whitney u test
    #draw rhe box plot with strip plot 

    df=pd.read_csv("./FIGURES/Figure1E_MajorMinorPatientCountBoxStripPlotInput.txt",sep="\t",header=None)
    df.columns=["Mutation","#Patients","Label"]
    from scipy import stats
    rvs1=df[df["Label"]=="KnownDriver[D]"]["#Patients"]
    rvs2=df[df["Label"]=="driver[d]"]["#Patients"]

    #x=stats.ttest_ind(rvs1,rvs2, equal_var = False)
    #If the p value returned is less than.05, then the null hypothesis is rejected and there is evidence that the data is not from a normally distributed population.
    stats.shapiro(rvs1)
    x=stats.mannwhitneyu(rvs1,rvs2, alternative="two-sided")
    stats.ttest_ind(rvs1,rvs2, equal_var = False)
    p=x.pvalue
    p
    "{:.3f}".format(p)
    #p= float("{:.3f}".format(p))
    #p
    import seaborn as sns
    palette ={"KnownDriver[D]": "blue", "driver[d]": "skyblue"}
    ax = sns.boxplot(x="Label", y="#Patients",data=df,color="white")
    ax = sns.stripplot(x="Label", y="#Patients",data=df,hue="Label",palette=palette)    
    ax.set_yscale("log")
    ax.set_xlabel("MannWhitney-p value={}".format(p))
    ax.set_ylabel("Number of Tumors")

    #117 known driver, 178 driver mutations
    ax.figure.savefig("./FIGURES/Figure1E_MajorMinorPatientCountBoxStripPlot.svg",fontsize=12,family='Arial',bbox_inches="tight",dpi=10000)
    ax.figure.savefig("./FIGURES/Figure1E_MajorMinorPatientCountBoxStripPlot.pdf",fontsize=12,family='Arial',bbox_inches="tight",dpi=10000)
    ##################

#figure1E = Figure1E_Latent_Driver_TumorCount_BoxStripPlot()



## Figure 2A

def Figure2A_Tissue_Bubble_Plot():
    #Load mutation_tumor, gene_tumor dictionaries
    import pickle, pandas as pd
    Gene_Mutation_Tumor_GEQ3_Dict = pickle.load(open("Gene_Mutation_Tumor_GEQ3_Dictionary.p", "rb"))
    Gene_All_Mutants = pickle.load(open("Gene_All_Mutants_Dictionary.p", "rb"))
    Tissue_Tumor_Dict = pickle.load(open("GENIE_TCGA_ClinicalVAF>12_5_Tissue_Tumor_Dictionary.p", "rb"))
    #load filtered double mutation tab delimited file
    

#############3
    df_new1 = pd.read_csv("Rev#3_Significant_doubles_VAF_Nonsense_filtered.txt", sep = "\t")

    df_new1.reset_index(drop = True, inplace = True)

    Genes = set(df_new1["Gene"].to_list()) #59 genescarries at least one doublet observed on 3 or more tumors
    sorted(list(Genes))
    #############
    Doublets = []
    DoubleComponents=set() 
    for ind in df_new1.index:
        m1 = df_new1.iloc[ind]["Gene"]+"_"+df_new1.iloc[ind]["Mut1"]
        m2 = df_new1.iloc[ind]["Gene"]+"_"+df_new1.iloc[ind]["Mut2"]
        Doublets.append((m1,m2))
        DoubleComponents.add(m1)
        DoubleComponents.add(m2)
    df_new1.iloc[0]["Gene"]
    ############
    #Tumors caarrying these 194 doublets
    Doublet_Tumor_Dict = { } #(gene_mut1,gene_mut2) are keys

    for double in Doublets:
        Doublet_Tumor_Dict[double] = Gene_Mutation_Tumor_GEQ3_Dict[double[0]].intersection(Gene_Mutation_Tumor_GEQ3_Dict[double[1]])
#############3
    ############

    #######################
    #for the 59 genes, get the list of double mutant tumors
    Gene_DoubleMutants_Dict = { g:set() for g in Genes}
    DoubleMutantTumors_All = set() #double  mutant all tumors, will be used for filtering tissues
    for d in Doublet_Tumor_Dict.keys():
        d1 = d[0]
        d2 = d[1]
        gen = d1.split("_")[0]
        Gene_DoubleMutants_Dict[gen] = Gene_DoubleMutants_Dict[gen].union(Doublet_Tumor_Dict[d])
        DoubleMutantTumors_All = DoubleMutantTumors_All.union(Doublet_Tumor_Dict[d])
    ##############

    ################Select tissues that carries at least 3 tumors with double muatations
    TissueFiltering =[]
    for tissue in Tissue_Tumor_Dict.keys():
        if len(DoubleMutantTumors_All.intersection(Tissue_Tumor_Dict[tissue])) >= 3:
            TissueFiltering.append(tissue)
            
    #For each gene find the intersection with each tissue for double mutant and all mutant tumors
    Gene_DoubleMutants_Tissue_Dict = {g:{} for g in Gene_DoubleMutants_Dict.keys()}

    for gene in Gene_DoubleMutants_Tissue_Dict.keys():
        for tissue in TissueFiltering:
            Gene_DoubleMutants_Tissue_Dict[gene][tissue] = len(Gene_DoubleMutants_Dict[gene].intersection(Tissue_Tumor_Dict[tissue]))

    #########3
    #Select Tissues where at least 3 tumors carries a double mutation

    ##############
    #For each gene find the intersection with each tissue for double mutant and all mutant tumors
    Gene_AllMutants_Tissue_Dict = {g:{} for g in Gene_All_Mutants.keys()}

    for gene in Gene_All_Mutants.keys():
        for tissue in TissueFiltering:
            Gene_AllMutants_Tissue_Dict[gene][tissue] = len(Gene_All_Mutants[gene].intersection(Tissue_Tumor_Dict[tissue]))

    ################

    Fraction_GeneDouble_All_Mutants_Tissue_Dict = { g:{} for g in Gene_All_Mutants.keys() }

    for gene in sorted(list(Genes)):
        for tissue in TissueFiltering:
            try:
                Fraction_GeneDouble_All_Mutants_Tissue_Dict[gene][tissue] = (Gene_DoubleMutants_Tissue_Dict[gene][tissue] / Gene_AllMutants_Tissue_Dict[gene][tissue])*100
            except ZeroDivisionError:
                print(gene)
                Fraction_GeneDouble_All_Mutants_Tissue_Dict[gene][tissue] = 0
    ############
    #Be careful about the thresholds
    with open("./FIGURES/Figure2A_Tissue_Bubble_Input_Filtered_Tissue_Capped.txt","a") as outfile:
        outfile.write("Gene"+"\t"+"Tissue"+"\t"+"DoubleMutant"+"\t"+"Fraction"+"\t"+"GeneMutant"+"\n")              

    for gene in Fraction_GeneDouble_All_Mutants_Tissue_Dict.keys():
        #if len(Gene_DoubleMutants_Dict[gene]) > = 0: #at least 5 double mutant tumors 
        for tissue in Fraction_GeneDouble_All_Mutants_Tissue_Dict[gene].keys():
            #if len(Tissue_Tumor_Dict[tissue]) > 100 : #use as threshold
            with open("./FIGURES/Figure2A_Tissue_Bubble_Input_Filtered_Tissue_Capped.txt","a") as outfile:
                outfile.write(gene+"\t"+tissue+"\t"+str(Gene_DoubleMutants_Tissue_Dict[gene][tissue])+"\t"+str(Fraction_GeneDouble_All_Mutants_Tissue_Dict[gene][tissue])+"\t"+str(Gene_AllMutants_Tissue_Dict[gene][tissue])+"\n")              





    import pandas as pd
    df = pd.read_csv("./FIGURES/Figure2A_Tissue_Bubble_Input_Filtered_Tissue_Capped.txt", sep="\t")
    df.columns
    df=df[['Gene', 'Tissue', 'DoubleMutant', 'Fraction']]    

    len(df)/21
     

    df = df.sort_values(["Tissue"], ascending=True) \
        .groupby(["Tissue"], sort=False) \
        .apply(lambda x: x.sort_values(['Gene'], ascending=True)) \
        .reset_index(drop=True)


    #plot kısmı

    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set()
    import matplotlib.pyplot as plt
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.cm as cm

    fig = plt.figure()
    fig.set_size_inches(16,9)
    #ax = fig.add_subplot(111)  
    ax = plt.gca()
    ax.set_facecolor('xkcd:white')
    ax.grid(color="black",which='both',axis='both', linestyle='--',linewidth=0.1)
    plt.figure(figsize=(16,9))
        
    import pandas as pd
    df['DoubleMutant'] = df['DoubleMutant'].astype(int)
    df['Fraction'] = df['Fraction'].astype(float)
    dfzero=df[df["DoubleMutant"]==0]
    len(set(dfzero["Tissue"].to_list()))
    df['DoubleMutant1'] = [0]*len(df)
    for ind in df.index:
        if df.loc[ind,'DoubleMutant']>=100:
            df.at[ind,'DoubleMutant1'] = 100
        else:
            df.at[ind,'DoubleMutant1'] = df.loc[ind,'DoubleMutant']
            
    scatter=ax.scatter(x="Gene", y="Tissue",c="Fraction",s=10*(df["DoubleMutant1"]*1),cmap=cm.Blues,edgecolors="black",vmin=0,vmax=20,data=df)
    type(df["DoubleMutant1"])
    ax.tick_params(axis='x', rotation=90)

    ax.set_xlabel('Double Mutant Genes  (n = {})'.format(len(Genes)))
    ax.set_ylabel('Tissues (At least 3 double mutant tumors)')
    ax.spines['bottom'].set_color('0.5')
    ax.spines['top'].set_color('0.5')
    ax.spines['right'].set_color('0.5')
    ax.spines['left'].set_color('0.5')
    ########custom legend
    # rankings, we only want to show 5 of them in the legend.
    legend1 = ax.legend(*scatter.legend_elements(num=10),
                        loc="upper left", title="Percentage", bbox_to_anchor=(1.01, 0.5))
    ax.add_artist(legend1)
    # Produce a legend for the price (sizes).  We use the *func* argument to supply the inverse of the function
    # used to calculate the sizes from above. The *fmt* ensures to show the fraction

    kw = dict(prop="sizes", num=8, color=scatter.cmap(0.7),   #fmt="$ {x:.2f}",
              func=lambda s: (s/1)/10)
    legend2 = ax.legend(*scatter.legend_elements(**kw),
                        loc="lower right", title="#Patients", bbox_to_anchor=(1.40, 0.5))
    plt.xticks(fontname = "Arial",fontsize=13)            
    plt.yticks(fontname = "Arial",fontsize=13)  # This argument will change the font.

    ax.figure.savefig("./FIGURES/Figure2A_Tissue_Bubble_Input_Filtered_Tissue_Capped_TissueBubblePlot_AllGenes_Filtered_Tissue_CappedSize_Capped.svg",bbox_inches='tight',dpi=12000)

    ax.figure.savefig("./FIGURES/Figure2A_Tissue_Bubble_Input_Filtered_Tissue_Capped_TissueBubblePlot_AllGenes_Filtered_Tissue_CappedSize_Capped.pdf",bbox_inches='tight',dpi=12000)


################3

def Figure2B_CancerSubtype_BubblePlot():
    
    #Tissue_Tumor_Info dictionary--->cancer subtype info
    c=0 #62566 tumors have tissue info out of 62567
    Tissue_Tumor_Dict = {}
    with open("GENIE_TCGA_ClinicalVAF>12_5_Tissue_Subtype_new.txt","r") as infile: 
        for line in infile:
            c+=1
            splitted = line.rstrip("\n").split("\t")
            tissue = splitted[-2] #select cancer subtype
            if tissue not in Tissue_Tumor_Dict.keys():
                Tissue_Tumor_Dict[tissue] = set()
                Tissue_Tumor_Dict[tissue].add(splitted[0])
            elif tissue in Tissue_Tumor_Dict.keys():
                Tissue_Tumor_Dict[tissue].add(splitted[0])

   #Load mutation_tumor, gene_tumor dictionaries
    import pickle, pandas as pd
    Gene_Mutation_Tumor_GEQ3_Dict = pickle.load(open("Gene_Mutation_Tumor_GEQ3_Dictionary.p", "rb"))
    Gene_All_Mutants = pickle.load(open("Gene_All_Mutants_Dictionary.p", "rb"))
    #Tissue_Tumor_Dict = pickle.load(open("GENIE_TCGA_ClinicalVAF>12_5_Tissue_Tumor_Dictionary.p", "rb"))
    #load filtered double mutation tab delimited file
    

#############3
    df_new1 = pd.read_csv("Rev#3_Significant_doubles_VAF_Nonsense_filtered.txt", sep = "\t")

    df_new1.reset_index(drop = True, inplace = True)

    Genes = set(df_new1["Gene"].to_list()) #59 genescarries at least one doublet observed on 3 or more tumors
    sorted(list(Genes))
    #############
    Doublets = []
    DoubleComponents=set() 
    for ind in df_new1.index:
        m1 = df_new1.iloc[ind]["Gene"]+"_"+df_new1.iloc[ind]["Mut1"]
        m2 = df_new1.iloc[ind]["Gene"]+"_"+df_new1.iloc[ind]["Mut2"]
        Doublets.append((m1,m2))
        DoubleComponents.add(m1)
        DoubleComponents.add(m2)
    df_new1.iloc[0]["Gene"]
    ############
    #Tumors caarrying these 194 doublets
    Doublet_Tumor_Dict = { } #(gene_mut1,gene_mut2) are keys

    for double in Doublets:
        Doublet_Tumor_Dict[double] = Gene_Mutation_Tumor_GEQ3_Dict[double[0]].intersection(Gene_Mutation_Tumor_GEQ3_Dict[double[1]])
#############3
    ############

    #######################
    #for the 59 genes, get the list of double mutant tumors
    Gene_DoubleMutants_Dict = { g:set() for g in Genes}
    DoubleMutantTumors_All = set() #double  mutant all tumors, will be used for filtering tissues
    for d in Doublet_Tumor_Dict.keys():
        d1 = d[0]
        d2 = d[1]
        gen = d1.split("_")[0]
        Gene_DoubleMutants_Dict[gen] = Gene_DoubleMutants_Dict[gen].union(Doublet_Tumor_Dict[d])
        DoubleMutantTumors_All = DoubleMutantTumors_All.union(Doublet_Tumor_Dict[d])
    ##############

    ################Select tissues that carries at least 3 tumors with double muatations
    TissueFiltering =[]
    for tissue in Tissue_Tumor_Dict.keys():
        if len(DoubleMutantTumors_All.intersection(Tissue_Tumor_Dict[tissue])) >= 3:
            TissueFiltering.append(tissue)
            
    #For each gene find the intersection with each tissue for double mutant and all mutant tumors
    Gene_DoubleMutants_Tissue_Dict = {g:{} for g in Gene_DoubleMutants_Dict.keys()}

    for gene in Gene_DoubleMutants_Tissue_Dict.keys():
        for tissue in TissueFiltering:
            Gene_DoubleMutants_Tissue_Dict[gene][tissue] = len(Gene_DoubleMutants_Dict[gene].intersection(Tissue_Tumor_Dict[tissue]))

    #########3
    #Select Tissues where at least 3 tumors carries a double mutation

    ##############
    #For each gene find the intersection with each tissue for double mutant and all mutant tumors
    Gene_AllMutants_Tissue_Dict = {g:{} for g in Gene_All_Mutants.keys()}

    for gene in Gene_All_Mutants.keys():
        for tissue in TissueFiltering:
            Gene_AllMutants_Tissue_Dict[gene][tissue] = len(Gene_All_Mutants[gene].intersection(Tissue_Tumor_Dict[tissue]))

    ################

    ################

    ################

    Fraction_GeneDouble_All_Mutants_Tissue_Dict = { g:{} for g in Gene_All_Mutants.keys() }

    for gene in sorted(list(Genes)):
        for tissue in TissueFiltering:
            try:
                Fraction_GeneDouble_All_Mutants_Tissue_Dict[gene][tissue] = (Gene_DoubleMutants_Tissue_Dict[gene][tissue] / Gene_AllMutants_Tissue_Dict[gene][tissue])*100
            except ZeroDivisionError:
                print(gene)
                Fraction_GeneDouble_All_Mutants_Tissue_Dict[gene][tissue] = 0
    ############
    ############
    #Be careful about the thresholds
    with open("./FIGURES/Figure2B_CancerSubtype_Bubble_Input_Filtered_Tissue_Capped.txt","a") as outfile:
        outfile.write("Gene"+"\t"+"Tissue"+"\t"+"DoubleMutant"+"\t"+"Fraction"+"\t"+"GeneMutant"+"\n")              

    for gene in Fraction_GeneDouble_All_Mutants_Tissue_Dict.keys():
        #if len(Gene_DoubleMutants_Dict[gene]) > = 0: #at least 5 double mutant tumors 
        for tissue in Fraction_GeneDouble_All_Mutants_Tissue_Dict[gene].keys():
            #if len(Tissue_Tumor_Dict[tissue]) > 100 : #use as threshold
            with open("./FIGURES/Figure2B_CancerSubtype_Bubble_Input_Filtered_Tissue_Capped.txt","a") as outfile:
                outfile.write(gene+"\t"+tissue+"\t"+str(Gene_DoubleMutants_Tissue_Dict[gene][tissue])+"\t"+str(Fraction_GeneDouble_All_Mutants_Tissue_Dict[gene][tissue])+"\t"+str(Gene_AllMutants_Tissue_Dict[gene][tissue])+"\n")              





    import pandas as pd
    df = pd.read_csv("./FIGURES/Figure2B_CancerSubtype_Bubble_Input_Filtered_Tissue_Capped.txt", sep="\t")
    df.columns
    df=df[['Gene', 'Tissue', 'DoubleMutant', 'Fraction']]    
        
     

    df = df.sort_values(["Tissue"], ascending=True) \
        .groupby(["Tissue"], sort=False) \
        .apply(lambda x: x.sort_values(['Gene'], ascending=True)) \
        .reset_index(drop=True)


    #plot kısmı

    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns; sns.set()
    import matplotlib.pyplot as plt
    import matplotlib.pyplot as plt
    import numpy as np
    import matplotlib.cm as cm

    fig = plt.figure()
    fig.set_size_inches(16,11)
    #ax = fig.add_subplot(111)  
    ax = plt.gca()
    ax.set_facecolor('xkcd:white')
    ax.grid(color="black",which='both',axis='both', linestyle='--',linewidth=0.1)
    plt.figure(figsize=(16,11))
        
    import pandas as pd
    df['DoubleMutant'] = df['DoubleMutant'].astype(int)
    df['Fraction'] = df['Fraction'].astype(float)
    dfzero=df[df["DoubleMutant"]==0]
    len(set(dfzero["Tissue"].to_list()))
    df['DoubleMutant1'] = [0]*len(df)
    for ind in df.index:
        if df.loc[ind,'DoubleMutant']>=70:
            df.at[ind,'DoubleMutant1'] = 70
        else:
            df.at[ind,'DoubleMutant1'] = df.loc[ind,'DoubleMutant']
            
    scatter=ax.scatter(x="Gene", y="Tissue",c="Fraction",s=10*(df["DoubleMutant1"]*1),cmap=cm.Blues,edgecolors="black",vmin=0,vmax=20,data=df)
    type(df["DoubleMutant1"])
    ax.tick_params(axis='x', rotation=90)

    ax.set_xlabel('Double Mutant Genes  (n = {})'.format(len(Genes)))
    ax.set_ylabel('Cancer Subtypes (At least 3 double mutant tumors)')
    ax.spines['bottom'].set_color('0.5')
    ax.spines['top'].set_color('0.5')
    ax.spines['right'].set_color('0.5')
    ax.spines['left'].set_color('0.5')
    ########custom legend
    # rankings, we only want to show 5 of them in the legend.
    legend1 = ax.legend(*scatter.legend_elements(num=10),
                        loc="upper left", title="Percentage", bbox_to_anchor=(1.01, 0.5))
    ax.add_artist(legend1)
    # Produce a legend for the price (sizes).  We use the *func* argument to supply the inverse of the function
    # used to calculate the sizes from above. The *fmt* ensures to show the fraction

    kw = dict(prop="sizes", num=8, color=scatter.cmap(0.7),   #fmt="$ {x:.2f}",
              func=lambda s: (s/1)/10)
    legend2 = ax.legend(*scatter.legend_elements(**kw),
                        loc="lower right", title="#Patients", bbox_to_anchor=(1.40, 0.5))
    plt.xticks(fontname = "Arial",fontsize=13)            
    plt.yticks(fontname = "Arial",fontsize=13)  # This argument will change the font.


    ax.figure.savefig("./FIGURES/Figure2B_CancerSubtype_BubblePlot_AllGenes_Filtered_Subtype_Capped.svg",bbox_inches='tight',dpi=12000)

    ax.figure.savefig("./FIGURES/Figure2B_CancerSubtype_BubblePlot_AllGenes_Filtered_Subtype_Capped.pdf",bbox_inches='tight',dpi=12000)


figure2A = Figure2A_Tissue_Bubble_Plot()


figure2B = Figure2B_CancerSubtype_BubblePlot()


def Figure3A():
    Gene_Mutation_Tumor_GEQ3_Dict = pickle.load(open("./Data/Gene_Mutation_Tumor_GEQ3_Dictionary.p", "rb"))
    Gene_Mutation_Tumor_LEQ2_Dict = pickle.load(open("./Data/Gene_Mutation_Tumor_LEQ2_Dictionary.p", "rb"))
    Gene_All_Mutants = pickle.load(open("./Data/DATA_CODES_new/OOP_OUTPUT/Gene_All_Mutants_Dictionary.p", "rb"))
    Tissue_Tumor_Dict = pickle.load(open("./Data/GENIE_TCGA_ClinicalVAF>12_5_Tissue_Tumor_Dictionary.p", "rb"))
    #load filtered double mutation tab delimited file
    import matplotlib.pyplot as plt   
    print(plt.rcParams['font.family'])
    # list of fonts in sans-serif
    plt.rcParams['font.sans-serif']
    # change the default font family
    plt.rcParams.update({'font.family':'Arial'})
    #############3
        
        
    # import matplotlib.pyplot as plt 
    # plt.rcParams['svg.fonttype'] = 'none' 
    import seaborn as sns

    #Read in data & create total column
    import pandas as pd
            

    df_significant= pd.read_cexcel("Supplementary_Table_1.xlsx", sep = "\t")

    df_new1 = df_significant

    df_significant.columns


    #read the file and draw the paired dot plots for different genes


    gene_name = "EGFR"

    df_new1.columns
    data=df_new1[df_new1["Gene"] == gene_name ]

    data.reset_index(drop=True,inplace=True)
    data.columns
    data["Fraction1"].round(decimals = 3)
    data["Fraction2"].round(decimals = 3)
    ####order labels and fraction and double mutation components
    data["Doublet_new"] = [None]*len(data)
    data["Fraction1_new"] = [0]*len(data)
    data["Fraction2_new"] = [0]*len(data)
    data["Class1_new"] = [None]*len(data)
    data["Class2_new"] = [None]*len(data)

    for ind in data.index:
        print(ind)
        if data.iloc[ind]["Fraction1"] < data.loc[ind,"Fraction2"] :
            print(11)
            data.at[ind,"Fraction1_new"] = data.loc[ind,"Fraction1"]
            data.at[ind,"Fraction2_new"] = data.loc[ind,"Fraction2"]
            data.at[ind,"Class1_new"] = data.loc[ind,"Class1"]
            data.at[ind,"Class2_new"] = data.loc[ind,"Class2"]
            data.at[ind,"Doublet_new"] = data.loc[ind,"Mut1"] + "/" + data.loc[ind,"Mut2"]
            
        else:
            print(22)
            data.at[ind,"Fraction1_new"] = data.loc[ind,"Fraction2"]
            data.at[ind,"Fraction2_new"] = data.loc[ind,"Fraction1"]
            data.at[ind,"Class1_new"] = data.loc[ind,"Class2"]
            data.at[ind,"Class2_new"] = data.loc[ind,"Class1"]
            data.at[ind,"Doublet_new"] = data.loc[ind,"Mut2"] + "/" + data.loc[ind,"Mut1"]

    data.loc[3,"Class1"]        


    final_data1 = data.sort_values(by=['Fraction1_new'], ascending=True) 
    final_data = final_data1.sort_values(by=['Activation'], ascending=True)       
    df_new1_1 = final_data.loc[(final_data['Activation'] >= 100)]     

    df_new1_2 = final_data.loc[(final_data['Activation'] >= 75)&(final_data['Activation'] < 100)]        
    df_new1_3 = final_data.loc[(final_data['Activation'] >= 50)&(final_data['Activation'] < 75)]        
    df_new1_4 = final_data.loc[(final_data['Activation'] < 50)]        

    df_new1_1.sort_values(by=['Fraction1_new'], ascending=True,inplace=True)
    df_new1_2.sort_values(by=['Fraction1_new'], ascending=True,inplace=True)
    df_new1_3.sort_values(by=['Fraction1_new'], ascending=True,inplace=True)
    df_new1_4.sort_values(by=['Fraction1_new'], ascending=True,inplace=True)
                              
                             
    final_data=pd.concat([df_new1_4,df_new1_3,df_new1_2,df_new1_1])

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    db = final_data[['Fraction1_new','Fraction2_new',"Class1_new","Class2_new","Doublet_new","DoubleMut#",'GeneMutant#']]
    db.sort_values(by=['Fraction2_new'], ascending=True,inplace=True)
    colors = {'StrongDriver':'purple','WeakDriver':'orchid', 'StrongLatentDriver':'royalblue', 'WeakLatentDriver':'skyblue'}
    plt.figure(figsize=(6,7))
    ax = plt.axes()
    # OR
    ax.set(facecolor = "white")
    ax.spines['left'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['bottom'].set_color('black')

    plt.vlines(x=db['Doublet_new'], ymin=db['Fraction1_new'], ymax=db['Fraction2_new'], color='grey', alpha=0.4,linewidth=0.4*db["DoubleMut#"])
    plt.scatter( db['Doublet_new'],db['Fraction1_new'], c=db['Class1_new'].map(colors), alpha=1,s=100)
    plt.scatter(db['Doublet_new'],db['Fraction2_new'], c=db['Class2_new'].map(colors), alpha=1,s=100 )
    #plt.yscale('log')
    plt.xticks(rotation=90)
    # The following two lines generate custom fake lines that will be used as legend entries:
    markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in colors.values()]
    plt.legend(markers,colors.keys(), numpoints=1)    
    plt.xlabel('Double Mutation Positions')
    plt.ylabel('Fraction (%) among {} Mutants'.format(gene_name))     
    plt.yticks(fontname = "Arial",fontsize=12)  # This argument will change the font.

    plt.xticks(fontname = "Arial",fontsize=12)            
    #plt.title("PIK3CA  Mechanism of Activation Plot")
    plt.xlabel('Double Mutation Constituents')
    plt.savefig("./Data/PairedDotPlot_PIK3CA.svg".format(gene_name),bbox_inches="tight",dpi=600)
    plt.savefig("./Data/PairedDotPlot_PIK3CA_.pdf".format(gene_name),bbox_inches="tight",dpi=600)

