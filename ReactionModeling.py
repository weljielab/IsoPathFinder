
# coding: utf-8

###########################################################################
#This script is designed to read in the list of MetaCyc atom-mapping solutions,
#.. and incorporate an isotope tracer (e.g. 13C6 Glucose) to enumerate
#.. the possible metabolite fates of the labeled carbons
#Inputs - reaction-linksAdd.txt, atom-mappings-smiles-flyAdd.dat, PubChemCpdMatchAdd.csv,  
#Outputs - Dictionary_FromRound_[0-9].pkl
#
#1) All reactions which contain a labeled metabolite are identified
#2) The carbons in the product metabolites are labeled for each reaction
#3) Those newly labeled product metabolites are added to the compound list
#4) The reactants, products, and ECs are built into a dictionary and exported, which can be used to build isotope labeling routes (separate .ipynb)
#5) This process is repeated for as many reactions as desired
#
# Note: since 13C02 is a common product, one can optionally remove CO2
#.. at the end of each reaction round - this step will GREATLY reduce the
#.. number of labeling possibilities
#
#This script is built to perform 13C tracing (or any isotope of carbon), however modifications can be made throughout the script to do e.g. 15N or 18O, etc..
#
#
# Important package versions used in this script (some aren't mentioned here, I suspect different versions of e.g. itertools will not have a substantial impact on the script):
# numpy: 1.11.3
# pandas: 0.19.2
# dill: 0.2.5
# pathos: 0.2.0
# Note this script was executed on an Ubuntu 15.10 Virtual Machine
###########################################################################

import numpy, pandas, os, sys, re, itertools, csv
from itertools import chain
from collections import defaultdict
import dill as pickle
from pathos.helpers import mp

#Import the reaction links of RXNs to ECs
#Note some reactions were manually added, thus an adapted MetaCyc file 'reaction-links.dat' is used here
MetaCycReactionLinks=pandas.read_csv('reaction-linksAdd.txt',sep='\t',skiprows=1)
MetaCycReactionLinks=MetaCycReactionLinks.drop(0)

#Trim reaction links
MetaCycReactionLinks=MetaCycReactionLinks[MetaCycReactionLinks.columns[0:2]]

#Rename columns
MetaCycReactionLinks.columns=['MetaCyc','EC']

#Fix some naming errors
MetaCycReactionLinks.MetaCyc=[re.sub('\?','+-RXN',x) for x in MetaCycReactionLinks.MetaCyc]

####################################################
#Note! Some ECs are manually removed 'ahead of time'
#These include Rubisco or bacterial enzymes, and can
#... be modified by the user
####################################################

###Change this list to your preference
ECDrop=['EC-4.1.2.22','EC-4.1.1.39']
###

MetaCycReactionLinks=MetaCycReactionLinks[~MetaCycReactionLinks['EC'].isin(ECDrop)]
MetaCycReactionLinks=MetaCycReactionLinks.reset_index(drop=True)


#Create a dictionary of MetaCyc RXN IDs with EC numbers as values
MetaCycReactionLinksDict=dict(zip(MetaCycReactionLinks.MetaCyc,MetaCycReactionLinks.EC))

#Read in the Smiles Matrix of reactions to compounds
#Note some mappings were added manually here
#If you want to read all MetaCyc reactions, use atom-mappings-smilesAdd.dat, or -mouseAdd or flyAdd for other organisms
#SmilesMap=pandas.read_csv('atom-mappings-smiles-humanAdd.dat',sep='\t',skiprows=1,header=None,usecols=[0,1])
SmilesMap=pandas.read_csv('atom-mappings-smiles-flyAdd.dat',sep='\t',skiprows=1,header=None,usecols=[0,1])
SmilesMap.columns=['Reaction','Compounds']

#Filter the SmilesMap based on the KEGG EC match
SmilesMap=SmilesMap[SmilesMap['Reaction'].isin(MetaCycReactionLinks.MetaCyc)]
SmilesMap=SmilesMap.reset_index(drop=True)

#Read in list of Smiles->'English' conversions
#This list is derived from the pubchempy python module
PubChemCpds=pandas.read_csv('PubChemCpdMatchAdd.csv',sep=',')

#Compiled expression matches - may shave some time off later given the number of string matches I perform
SmilesStringRe1=re.compile(r':\d?\d]*')
SmilesStringRe2=re.compile(r'\[|\]')
SmilesStringRe3=re.compile(r'^\'|\'$')
MetaCycSplitter=re.compile(r'(>>|\.)')
CarbFinder=re.compile(r'C\*?')
LabelFinder=re.compile(r'C\*:\d+')
SpaceRemover=re.compile(r'^ ')

#Dictionary to pull out 'English name' from the Smiles format
PubChemCpdsDict=dict(zip(PubChemCpds.Smiles,PubChemCpds.Name))

#Function to get a smiles string from the MetaCyc-specific Smiles annotation
#Using compiled expressions
def SmilesFromMetaCycCpd(OneCpdString):
    
    #Get rid of the bracketing
    OneCpdString=SmilesStringRe1.sub(repl='',string=OneCpdString)
    
    OneCpdString=SmilesStringRe2.sub(repl='',string=OneCpdString)
    
    #remove extra apostrophes if they exist
    if OneCpdString.startswith('\'') and OneCpdString.endswith('\''):
        OneCpdString=SmilesStringRe3.sub(repl='',string=OneCpdString)
        
    return(OneCpdString)

############################################################################
#Define the isotope tracer of interest here!

#If you need help finding the right SMILES for your input tracer of interest
#... check the PubChemCpds list
############################################################################

#Pulling any hexopyranose from PubChem
GLUCOSESMILES=['C(O)C1(C(O)C(O)C(O)C(O)O1)','C(C1(C(C(C(C(O1)O)O)O)O))O',
              'C(O)C1(OC(C(C(C1O)O)O)O)','C(C1(OC(C(C(C1O)O)O)O))O',
              'C(O)C1(OC(O)C(O)C(O)C(O)1)','C(C1(C(O)C(O)C(O)C(O)O1))O']

#Defining the initial list of labeled compounds - this list will seed the 1st round of metabolic reactions
LabListInitial=['C*(O)C*1(C*(O)C*(O)C*(O)C*(O)O1)','C*(C*1(C*(C*(C*(C*(O1)O)O)O)O))O',
              'C*(O)C*1(OC*(C*(C*(C*1O)O)O)O)','C*(C*1(OC*(C*(C*(C*1O)O)O)O))O',
              'C*(O)C*1(OC*(O)C*(O)C*(O)C*(O)1)','C*(C*1(C*(O)C*(O)C*(O)C*(O)O1))O']

#Dictionary of MetaCyc Rxn to MetaCyc Compound list split up
SmilesMapDict=dict(zip(SmilesMap.Reaction,SmilesMap.Compounds))

#Dictionary of MetaCyc Rxn to Compound items in SMILES format
SmilesMapSmilesSplit=SmilesMap.Compounds.str.split('[.]|[>>]')
Reducer=lambda x:[SmilesFromMetaCycCpd(y) for y in x]
SmilesMapSmilesSplit=SmilesMapSmilesSplit.apply(Reducer)
SmilesMapSmilesCpdDict=dict(zip(SmilesMap.Reaction,SmilesMapSmilesSplit))

#Master compound list built from SmilesMap
MasterCpdList=list(chain.from_iterable(SmilesMapSmilesSplit))
MasterCpdList=list(set(MasterCpdList))

#Master dictionary of the reactions which contain any one of the 12k compounds - useful for GenerateNewReaction function
CpdToRxnDict={}
for cpd in range(len(MasterCpdList)):
    xlist=[k for (k,v) in SmilesMapSmilesCpdDict.items() if MasterCpdList[cpd] in v] #56s for cpd list of 12k
    CpdToRxnDict.setdefault(MasterCpdList[cpd],[])
    CpdToRxnDict[MasterCpdList[cpd]].extend(xlist)

#Get the index in a reaction that is already split by its .'s and >>'s of the >> (i.e. what splits the reactants and products)
def GetTheReactantsLength(OneSplitSMILESReaction):
    for i in range(0,len(OneSplitSMILESReaction)):
        if len(OneSplitSMILESReaction[i])==0:
            return i

#Dictionary of the point of equilibration sign based on str.split('[.]|[>>]'), with keys as the MetaCyc RXN IDs
ReactantLengthDict=dict(zip(SmilesMap.Reaction,SmilesMapSmilesSplit.apply(GetTheReactantsLength)))

#Insert a labeled compound into a MetaCyc reaction by finding its isotope-stripped version
def InsertLabeledMetabIntoString(BigString,LabMetab,StripMetab):
    
    #Find the carbons
    Carbs=CarbFinder.findall(LabMetab)
    
    #Iterative list of the carbon locations with a label
    match=[j for j in range(len(Carbs)) if Carbs[j].__contains__('C*')]

    #Split the MetaCyc reaction compound string into reactants and products
    TempRxn=MetaCycSplitter.split(str(BigString))
    
    #Convert each item to Smiles format
    TempRxn1=[SmilesFromMetaCycCpd(x) for x in TempRxn]
    
    for metab in range(len(TempRxn)):
        
        #Find the carbon locations
        CIndicies=[m.start() for m in re.finditer('C',TempRxn[metab])]
        
        #Listify
        ListMetab=list(TempRxn[metab])
        
        #If any of the items from the MetaCyc reactants/products list match the metabolite which contains the isotope..
        if TempRxn1[metab]==StripMetab:
            
            for location in range(len(match)):
                
                #Insert the isotope at the right location
                ListMetab[CIndicies[match[location]]]='C*'
                
                #Rejoin the list
                ConvertedMetab=''.join(ListMetab)
                
                #Replace old metabolite with newly labeled form
                TempRxn[metab]=ConvertedMetab
            
        #Convert back to MetaCyc format
        BigString=''.join(item for item in TempRxn)
        
    return(BigString)

#Take the list of labeled compounds, and find all the Metacyc reactions that contain the compounds, build a matrix
def BuildReactionMatrixFromCpdList(ListOfLabeledCpds):
    
    MasterRxnList=[]
    MasterCompoundList=[]
    
    #Strip all labeled metabolites of isotopes for matching purposes
    StrippedList=[x.replace('C*','C') for x in ListOfLabeledCpds]
    
    for labcpd,nolabcpd in zip(ListOfLabeledCpds,StrippedList):
        
        try:
            
            #Try using the dictionary to generate a list of reactions for each given metabolite
            RxnList=CpdToRxnDict[nolabcpd]
            
            #Get the Metacyc compounds from each reaction
            RxnCompounds=[SmilesMapDict[x] for x in RxnList]
            
            #Insert the labeled compound into the Metacyc compounds list
            RxnCompounds=[InsertLabeledMetabIntoString(x,labcpd,nolabcpd) for x in RxnCompounds]
            
        except:
            RxnList=[]
            RxnCompounds=[]
        
        #Build out the lists
        MasterRxnList.extend(RxnList)
        MasterCompoundList.extend(RxnCompounds)

    #Turn into Series' for later apply functions and indexing
    PrepMatrix=pandas.Series([MasterRxnList,MasterCompoundList])
    
    #Build DataFrame
    PrepMatrix=BuildReactionMatrixExport(PrepMatrix)

    return(PrepMatrix)

#Function to build the reaction matrix export
def BuildReactionMatrixExport(PreLabeledReactionMatrix):
    
        #Build new data frame with same rows as ReactionMatrix
        DataFrameExport=pandas.DataFrame(index=numpy.arange(len(PreLabeledReactionMatrix[0])),columns=['Reaction','Compounds','Reactants','Products'])
        
        #Fill the reaction column from the lists of reactions and compounds with isotope addition from 'BuildReactionMatrix' function
        DataFrameExport['Reaction']=PreLabeledReactionMatrix[0]
        DataFrameExport['Compounds']=PreLabeledReactionMatrix[1]
        
        del PreLabeledReactionMatrix
        
        return(DataFrameExport)

#Take list of carbon numbers, find them in the metabolite of interest, and replace Cs with C*s (i.e. with isotope)
def MiniInsertLabel(CarbonRegexList,Metabolite):
    for num in CarbonRegexList:
        Metabolite=re.sub(pattern='\\b'+'C:'+num+'\\b',repl='C*:'+num,string=Metabolite)
    return(Metabolite)

#Take labeled reactants, find carbon numbers, map to products and label those carbons
def GenericLabelingReactionMatrixApply(UnlabeledReactionMatrix):
    
    #Split into metabolite items
    UnlabeledReactionMatrix=MetaCycSplitter.split(UnlabeledReactionMatrix)
    
    #Get all reactants with label
    ReactantList=[item for item in UnlabeledReactionMatrix if '*' in item]
    
    #find carbom numbers to get labeled
    FindThese=LabelFinder.findall(str(ReactantList))
    
    #Return only the numbers that are assigned to the carbons
    FindThese=[x.replace('C*:','') for x in FindThese]
    
    #Get cpds not currently labeled
    HoldList=[item for item in UnlabeledReactionMatrix if not '*' in item]
    
    #Return the smiles format of reactants labeled cpds
    ReactantList=[SmilesFromMetaCycCpd(x) for x in ReactantList]
    
    #Label any metabolite that contains a carbon number from FindThese
    HoldList=[MiniInsertLabel(FindThese,x) for x in HoldList]
    
    #Pull out only metabolites which obtained a labeled (i.e. the Products)
    ProductList=[item for item in HoldList if '*' in item]
    
    #Return Smiles
    ProductList=[SmilesFromMetaCycCpd(x) for x in ProductList]
    
    #Build tuple of reactants and products
    UnlabeledReactionMatrix=(ReactantList,ProductList)
    
    del ReactantList,ProductList,HoldList

    return(UnlabeledReactionMatrix)

#Unpack the tuple of reactants and products, fill in the appropriate columns of the dataframe
def UnpackAndFill(ReactedMatrix):
    Reactants,Products=zip(*ReactedMatrix.Compounds)
    ReactedMatrix.Reactants=Reactants
    ReactedMatrix.Products=Products
    del Reactants,Products
    ReactedMatrix.drop('Compounds',axis=1,inplace=True)
    return(ReactedMatrix)

#New function, take products from Labeled reaction matrix, make new labeled list
def NewLabeledCpdList(ReactedMatrix):
    
    #Build list of cpds from the products column
    NewList=list(chain.from_iterable(ReactedMatrix.NewCpds))

    #Some metabolites had additional spaces which messed up the matching
    NewList=[re.split(' ',x) for x in NewList]
    NewList=list(chain.from_iterable(NewList))
    NewList=[x for x in NewList if '*' in x]
    
    ######################################################
    #Important!
    #The user can add metabolites (by their smiles format)
    #... to remove at each round
    #
    #Adding too many metabolites here will probably slow
    #... the program significantly, fyi..
    ######################################################
    
    #Remove 13CO2 at the end of each round! Otherwise a large number of the reaction possibilities contain a CO2 which is probably not aligned with experimental reality 
    NewList=[x for x in NewList if re.sub('\*','',x)!='C(=O)=O']
    
    #Remove labeled sorboses, doesn't seem to have physiological relevance
    NewList=[x for x in NewList if re.sub('\*','',x)!='C(O)C(=O)C(O)C(O)C(O)CO']
    
    #Keep unique elements of list
    NewList=list(set(NewList))
    
    return(NewList)

#Use MetaCycReactionLinksDict for reaction conversion from MetaCyc's RXN to EC
def ReactionConvert(RunReactionMatrix):
    if str(MetaCycReactionLinksDict[RunReactionMatrix])!='nan':
        RunReactionMatrix=MetaCycReactionLinksDict[RunReactionMatrix]
    return(RunReactionMatrix)

#Name converter from Smiles to 'English'
def ReactantConvertOne(x):
    try:
        x=PubChemCpdsDict[x]
    except:
        x=''
    return(x)

#Convert a compound string so that it can be matched to an 'English' name
#Compiled patterns
ChargeRemover=re.compile(r'\+|\-')
ExtraPostasRemover=re.compile('^\'|\'$')

def LabelStrippersOne(x):
    
    #Remove labels #Not needed when keeping smiles
    #x=x.replace('*','')
    
    #Remove extra Rs? #Not needed when keeping smiles
    #x=x.replace('R','')
    
    #Remove charges
    x=ChargeRemover.sub('',string=str(x))
    
    #Remove any rextra apostrophes
    x=ExtraPostasRemover.sub('',string=str(x))
    
    return(x)

#Convert a labeled metabolite string into 'English isotopomer name (e.g. Glucose M+6)
def CombinedConvertMetabs(OneLabeledMatrixCpd):
    
    #Get number of label counts # Not needed when keeping smiles
    #OneLabeledMatrixCpdCarbs=OneLabeledMatrixCpd.count('*')
    
    #Strip labels and charges
    OneLabeledMatrixCpdsEnglish=LabelStrippersOne(OneLabeledMatrixCpd)
    
    #Convert to English using metabolite dictionary # Not needed when keeping smiles format
    #OneLabeledMatrixCpdsEnglish=ReactantConvertOne(OneLabeledMatrixCpdsEnglish)
    
    #Format string to use both english name and label number # Not needed when keeping smiles format
    #OneLabeledMatrixCpd=str('{0} M+{1}').format(OneLabeledMatrixCpdsEnglish,OneLabeledMatrixCpdCarbs)
    
    #Changed from return(OneLabeledMatrixCpd) to keep smiles names
    return(OneLabeledMatrixCpdsEnglish)

#Function to call the conversion of reactants and products to english
def UpdatedMetabConvert(LabeledMatrixCpds):
    LabeledMatrixCpds=[CombinedConvertMetabs(x) for x in LabeledMatrixCpds]
    return(LabeledMatrixCpds)

#Trim reaction matrix to drop nonmatching metabolites (i.e. they can't be named)
def TrimExportReactionMatrix(ExportReactionMatrix):

    #Drop the nans
    ExportReactionMatrix=ExportReactionMatrix.dropna(subset=['Reactants'])
    ExportReactionMatrix=ExportReactionMatrix.dropna(subset=['Products'])
    
    #drop row indicies that have a blank Labeled Reactants or Products cell
    ExportReactionMatrix=ExportReactionMatrix.drop(ExportReactionMatrix[ExportReactionMatrix['Reactants'].map(len)==0].index,axis=0)
    ExportReactionMatrix=ExportReactionMatrix.drop(ExportReactionMatrix[ExportReactionMatrix['Products'].map(len)==0].index,axis=0)

    return(ExportReactionMatrix)

#If metabolite didn't get an english name match, drop it
def DropMissingMetabs(ConvertedMatrixColumnCpd):

    ConvertedMatrixColumnCpd=[item for item in ConvertedMatrixColumnCpd if not re.search('^ M+[0-9]*',item)]

    return(ConvertedMatrixColumnCpd)

#convert to english, build the dictionary for round 1
def ConvertMatrixToEnglish(LabeledReactionMatrix):
    
    #Convert RXN to EC
    LabeledReactionMatrix.Reaction=LabeledReactionMatrix.Reaction.apply(ReactionConvert)
   
    ##############
    #Update - don't convert names, do that at the end of the IsoPathFinder Step, and instead keep the SMILES
    #.. format instead
    
    #Convert reactants and products to english M+x isotopomer names #Added back to include the charge strip
    LabeledReactionMatrix.Reactants=LabeledReactionMatrix.Reactants.apply(UpdatedMetabConvert)    
    LabeledReactionMatrix.Products=LabeledReactionMatrix.Products.apply(UpdatedMetabConvert)

    #Remove unnamed metabolites
    ##LabeledReactionMatrix.Reactants=LabeledReactionMatrix.Reactants.apply(DropMissingMetabs) #622ms, 25% faster than before
    ##LabeledReactionMatrix.Products=LabeledReactionMatrix.Products.apply(DropMissingMetabs) #622ms, 25% faster than before

    #Drop empty fields
    ##LabeledReactionMatrix=TrimExportReactionMatrix(LabeledReactionMatrix)
    ##LabeledReactionMatrix=LabeledReactionMatrix.reset_index(drop=True)
    ############
    
    return(LabeledReactionMatrix)

#Build dictiontary/hashable results from the labeling matrix
def MakingDict(LabeledRxnMatrix):
    
    Dictionary=defaultdict(lambda: defaultdict(list))

    for rxns in range(len(LabeledRxnMatrix)):
        
        #Remove any extra spaces
        Reactants=[SpaceRemover.sub('',x) for x in LabeledRxnMatrix.Reactants[rxns]]

        Products=[SpaceRemover.sub('',x) for x in LabeledRxnMatrix.Products[rxns]]
        
        #Add values to the Product keys, both the reactants as well as the reactants' as keys to the enzymes
        if len(Products)==1:
            if LabeledRxnMatrix.Reaction[rxns] in Dictionary[str(Products).strip('[]|\'')][Reactants[0]]:
                pass
            else:
                Dictionary[str(Products).strip('[]|\'')][Reactants[0]].append(LabeledRxnMatrix.Reaction[rxns])
           
        if len(Products)>1:
            for cpd in range(len(Products)):
                if LabeledRxnMatrix.Reaction[rxns] in Dictionary[str(Products[cpd]).strip('[]|\'')][Reactants[0]]:
                    pass
                else:
                    Dictionary[str(Products[cpd]).strip('[]|\'')][Reactants[0]].append(LabeledRxnMatrix.Reaction[rxns])

    
    return(Dictionary)


#What I had found was that for the first couple reaction rounds, the dataframes
#.. were generally small and only took a few seconds
#.. However, in later rounds I needed to use parallelization (especially true if running with the entire MetaCyc model)
#.. which is the code below (with specific parallel functions)

NextRoundReaction=BuildReactionMatrixFromCpdList(LabListInitial)
NextRoundReaction.Compounds=NextRoundReaction.Compounds.apply(GenericLabelingReactionMatrixApply)
NextRoundReaction=UnpackAndFill(NextRoundReaction)
NextRoundReaction['NewCpds']=NextRoundReaction.Products
NextRoundCpdList=NewLabeledCpdList(NextRoundReaction)
NextRoundReaction=ConvertMatrixToEnglish(NextRoundReaction)
DictionaryRound=MakingDict(NextRoundReaction)
output=open('Dictionary_FromRound_1.pkl','wb', -1) #Change
pickle.dump(DictionaryRound,output,protocol=4)
print('Round 1 Done') #Change


#Build and run the reaction, also append the new labeled cpd list
def PipePandaProcess(CpdList):
    
    #build reaction matrix
    ReactionMatrix=BuildReactionMatrixFromCpdList(CpdList)
    
    #Run labeling reaction
    ReactionMatrix.Compounds=ReactionMatrix.Compounds.apply(GenericLabelingReactionMatrixApply)
    
    #Prep step
    ReactionMatrix=UnpackAndFill(ReactionMatrix)
    
    #Add new column for labeled list export later
    ReactionMatrix['NewCpds']=ReactionMatrix.Products
    
    #Convert matrix to English
    ReactionMatrix=ConvertMatrixToEnglish(ReactionMatrix)
        
    return(ReactionMatrix)

#Successful write out in parallel function
#Note this function will write out a temporary dataframe and read it back in pieces to build the reaction dictionary, which alleviates memory concerns by writing to disk
def ParallelRxnBuildWriteOut(LabeledList,BuildFunction,RoundNum):
    
    #Split into - using 16 cores, break up evenly
    SplitList=numpy.array_split(LabeledList,16*32)
    
    #Fill new list of labeled cpds
    NewCpdList=[]
    
    #Build multiprocessing pool object with 16 processors
    pooler=mp.Pool(16)

    #Set up the new file
    with open('TempPandaDF.csv','w') as fp:
        
        #For each result in the Build function using an item from SplitList
        for result in pooler.imap(BuildFunction,SplitList):
            
            #Building new cpd list
            NewCpds=NewLabeledCpdList(result)
            NewCpdList.extend(NewCpds)
            
            #Each result is a Pandas Object, so write it to csv
            result.to_csv(fp,index=False,header=False)

    
    pooler.close()
    pooler.join()

    #Unique items
    NewCpdList=list(set(NewCpdList))
    
    #Write File, optional
    #CpdFile='LabeledCpds_FromRound_{0}.csv'.format(RoundNum) #Change
    #with open(CpdFile,'w') as output:
    #    writer=csv.writer(output,lineterminator='\n')
    #    for val in NewCpdList:
    #        writer.writerow([val])
    
    #Return CpdList for next round
    return(NewCpdList)

#Successful write out in parallel function
def ParallelRxnBuildWriteOutRd2(LabeledList,BuildFunction,RoundNum):
    
    #No chunks for round 2
    SplitList=numpy.array_split(LabeledList,16)
    
    #Fill new list of labeled cpds
    NewCpdList=[]
    
    #Build multiprocessing pool object with num of cpus
    pooler=mp.Pool(16)

    #Set up the new file
    with open('TempPandaDF.csv','w') as fp:
        
        #For each result in the Build function using an item from SplitList
        for result in pooler.imap(BuildFunction,SplitList):
            
            #Building new cpd list
            NewCpds=NewLabeledCpdList(result)
            NewCpdList.extend(NewCpds)
            
            #Each result is a Pandas Object, so write it to csv
            result.to_csv(fp,index=False,header=False)

    
    pooler.close()
    pooler.join()

    #Unique items
    NewCpdList=list(set(NewCpdList))
    
    #Write File, optional
    #CpdFile='LabeledCpds_FromRound_{0}.csv'.format(RoundNum) #Change
    #with open(CpdFile,'w') as output:
    #    writer=csv.writer(output,lineterminator='\n')
    #    for val in NewCpdList:
    #        writer.writerow([val])
    
    #Return CpdList for next round
    return(NewCpdList)

def ConvertBackForDict(PandaDF):
    #Convert back to original pandas formatting
    PandaDF.Reactants=PandaDF.Reactants.apply(ConvertBackForDictSub)
    PandaDF.Products=PandaDF.Products.apply(ConvertBackForDictSub)
    return(PandaDF)

def ConvertBackForDictSub(PandaDFColumn):
    PandaDFColumn=re.sub('\'|\[|\]','',PandaDFColumn)
    PandaDFColumn=re.split(', ',PandaDFColumn)
    return(PandaDFColumn)

def AddingDict(ExistingDictionary,LabeledRxnMatrix):
    
    LabeledRxnMatrix=LabeledRxnMatrix.reset_index(drop=True)
    
    for rxns in range(len(LabeledRxnMatrix)):
        
        #Remove any extra spaces
        Reactants=[SpaceRemover.sub('',x) for x in LabeledRxnMatrix.Reactants[rxns]]

        Products=[SpaceRemover.sub('',x) for x in LabeledRxnMatrix.Products[rxns]]
        
        #Add values to the Product keys, both the reactants as well as the reactants' as keys to the enzymes
        if len(Products)==1:
            if LabeledRxnMatrix.Reaction[rxns] in ExistingDictionary[str(Products).strip('[]|\'')][Reactants[0]]:
                pass
            else:
                ExistingDictionary[str(Products).strip('[]|\'')][Reactants[0]].append(LabeledRxnMatrix.Reaction[rxns])
           
        if len(Products)>1:
            for cpd in range(len(Products)):
                if LabeledRxnMatrix.Reaction[rxns] in ExistingDictionary[str(Products[cpd]).strip('[]|\'')][Reactants[0]]:
                    pass
                else:
                    ExistingDictionary[str(Products[cpd]).strip('[]|\'')][Reactants[0]].append(LabeledRxnMatrix.Reaction[rxns])

    return(ExistingDictionary)

#Since the data frames got too big, I opted to write out to file then read back in the results in chunks to build the compiled dictionaries of
#.. the 13C-labeled metabolic network


def BuildDictionaryFromCsv(RoundNum):
    #start with blank Dict
    RoundDictionary=defaultdict(lambda: defaultdict(list))
    
    #Read in the Pandas DF by 10000 row chunks - when I build the loops, the TempPandaDF will be
    #a dataframe written to disk as a temporary hold for all the reactions run in a given round
    
    for chunk in pandas.read_csv('TempPandaDF.csv',sep=',',index_col=False,skip_blank_lines=True,names=['Reaction','Reactants','Products','NewCpds'],chunksize=10000):
        #Convert the chunk into a format to be built into dictinoary (i.e. original pandas format)
        chunk=ConvertBackForDict(chunk)
        #Use the chunk to add into the exisitng dictionary
        RoundDictionary=AddingDict(RoundDictionary,chunk)
    
    #Write out dictionary
    output=open('Dictionary_FromRound_{0}.pkl'.format(RoundNum),'wb',-1)
    pickle.dump(RoundDictionary,output,protocol=4)
    print('Round {0} Done'.format(RoundNum))
    
    #Remove the TempPandaDF so it doesn't conflict with the next round
    os.remove('TempPandaDF.csv')

##Run parallel reactions
#Note I have these separated by round, one could write a simple loop as well, I liked knowing how long each round took when I was speed-optimizing the code

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,2)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(2)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,3)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(3)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,4)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(4)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,5)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(5)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,6)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(6)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,7)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(7)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,8)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(8)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,9)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(9)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,10)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(10)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,11)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(11)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,12)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(12)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,13)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(13)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,14)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(14)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,15)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(15)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,16)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(16)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,17)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(17)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,18)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(18)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,19)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(19)

NextRoundCpdList=ParallelRxnBuildWriteOutRd2(NextRoundCpdList,PipePandaProcess,20)
#Build the dictioanry for a given round, using the TempPandaDF.csv file on disk
BuildDictionaryFromCsv(20)

#Keep going if you want, however the path search tended to get overwhelming by 15 rounds, so more calculations may be overkill here

