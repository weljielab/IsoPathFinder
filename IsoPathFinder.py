
# coding: utf-8

###############################################################################################################################
#Tracing Paths - this script is designed for glucose tracing
#Inputs - Dictionary_FromRound_[0-9].pkl, CpdConvertDict.csv, PubChemCpdMatchAdd.csv
#Outputs - Serine M+3_Paths_9Rxns.csv, Glutamine M+2_Paths_14Rxns.csv (examples)

#To use this on your own, make sure you've run the ReactionModeling script first to compile the tracer results, before you can
#.. query the network to construct paths

#This script will use the dictionaries from the local directory, so if you build e.g. 13C6 glucose tracer results, and want to
#... query that network, then have those dictionaries here

#This needs to be run in Linux, and multiple cores (e.g. 4 or more) with 16GB+ of RAM is HIGHLY recommended
#Note you will need to install pathos and dill packages, and pandas needs to be v0.19 or higher
###############################################################################################################################

import numpy, pandas, os, sys, re, itertools, csv, gc

from itertools import chain, repeat, combinations, permutations, cycle

from collections import defaultdict

from pathos.helpers import mp

from contextlib import closing

import dill as pickle

#Load all the dictionaries
ResultsList=[]
#You will likey not want to test beyond 15 steps anyway given the number of possible paths, but you can change the number of 
#.. dictionaries to read in here
for MyDict in range(16):
    ResultsList.append(pickle.load(open('Dictionary_FromRound_{0}.pkl'.format(MyDict+1),'rb')))

##Read in file for name conversion, make dictionaries
CpdConvert=pandas.read_csv('CpdConvertDict.csv')
CpdConvert=CpdConvert.dropna()
MetaToCommon=dict(zip(CpdConvert.Name,CpdConvert.Common))
CommonToMeta=dict(zip(CpdConvert.Common,CpdConvert.Name))

#Read in file for pubchem smiles to common name
PubChemCpds=pandas.read_csv('PubChemCpdMatchAdd.csv',sep=',')
PubChemCpdsDict=dict(zip(PubChemCpds.Smiles,PubChemCpds.Name))

#This form of NameToPubChemSmiles builds lists of Common name:SMILES in case 
#.. there are multiple entries for the same metabolite (e.g. 'DL-Serine')
NameToPubChemSmiles=defaultdict(list)
for cpd in range(len(PubChemCpds)):
    NameToPubChemSmiles[PubChemCpds.Name[cpd]].append(PubChemCpds.Smiles[cpd])

#Specify input
tracerinput = 'alpha-D-Glucose M+6'

#Old function - get the smiles from a pubchem name ID
def SmilesFromMetaCycCpd(OneCpdString):
    
    #Get rid of the bracketing
    OneCpdString=SmilesStringRe1.sub(repl='',string=OneCpdString)
    
    OneCpdString=SmilesStringRe2.sub(repl='',string=OneCpdString)
    
    #remove extra apostrophes if they exist
    if OneCpdString.startswith('\'') and OneCpdString.endswith('\''):
        OneCpdString=SmilesStringRe3.sub(repl='',string=OneCpdString)
        
    return(OneCpdString)

#Convert user input to Smiles - return list of labeled SMILES options
def SmilesOptionsFromUserInput(userinput):
    try:
        #Convert the user input to the MetaCyc-compatible name
        Conversion1=CommonToMeta[re.sub(' M\+\d','',userinput)]
        Conversion2=NameToPubChemSmiles[Conversion1]
        return(Conversion2)
    except:
        print('Oops, no match of the user input to Pubchem Smiles, check your spelling')


#Possibilities level off for metabolites, e.g. Serine starts duplicating after round 12, so only need to calculate paths until then
def FinalRoundMetabLength(ResultsDictionaries,MetabName):
    ResultDFs=[]
    
    #Build the dataframes for each round
    for resultdict in range(len(ResultsDictionaries)):
        ResultDFs.append(BuildOneRoundPath(ResultsList[resultdict],MetabName))
    
    #Check if consecutive rounds are equal
    for rxnround in reversed(range(len(ResultDFs))):
        if all(ResultDFs[rxnround][2].isin(ResultDFs[rxnround-1][2]))==False and len(ResultDFs[rxnround-1]>0):
            Length=rxnround
            break
    return(Length)

#Combinations of isotopomers may exist for the same compound 
#    e.g. metabolite X M+2 could be C*(C(C*)OOO), or C*(C*(C)OOO), so we need to run a combinations function to get all possible isotopologues
#If you want to run a single isotopomer (say you know positional 13C labeling), then you can skip this part and run the BuildPathsCoreParallel
#... function manually


#Need number of carbons, number of user input labeled carbons
def IsotopomerCombinations(CpdSmiles,NumOfLabels):
    
    MasterIsotopomerList=[]
    
    for cpd in CpdSmiles:
    
        #Find all the carbons in the original Smiles input
        CarbonList=re.findall('C',cpd)
    
        #Set up the first permutation of labels by adding the specified number of 13Cs
        LabelCarbonList=[re.sub('C','C*',x) for x in CarbonList[0:NumOfLabels]]
    
        #Replace the carbons with labeled carbons
        CarbonList[0:NumOfLabels]=LabelCarbonList
    
        #Make permutations of isotopomers
        CarbonPerms=[x for x in permutations(CarbonList)]
    
        #Build into list
        CarbonPerms=[list(x) for x in list(set(CarbonPerms))]
    
        #Find the indicies of carbons in original Smiles, this will get replaced by 13C
        ReplList=[m.start() for m in re.finditer('C',cpd)]
    
        #Listify the string for positional replacement
        ListifySmiles=[x for x in cpd]
    
        #List to be built of all permutations of isotopomers
        NewListOfLabeledSmiles=[]
    
        #For each isotopomer permutation
        for Perms in range(len(CarbonPerms)):
        
            #Ordered list of the given isotopomer
            LabList=CarbonPerms[Perms]
        
            #For each index to be replaced
            for position in range(len(ReplList)):
            
                #Replace with C* if necessary
                ListifySmiles[ReplList[position]]=LabList[position]
            
            #Rebuild into string
            cpd=''.join(ListifySmiles)
        
            #Add string to all possible isotopomers
            NewListOfLabeledSmiles.append(cpd)
        
        MasterIsotopomerList.extend(NewListOfLabeledSmiles)
        
    return(MasterIsotopomerList)


#Calculate longest path for each isotopomer - cap at 15 due to number of possibilities
def CalculateIsotopomerLength(ResultsDictionaries,IsotopomerComboList):
    Lengths=[]
    for isotope in IsotopomerComboList:
        try:
            Lengths.append(FinalRoundMetabLength(ResultsDictionaries,isotope))
        except:
            pass
    if max(Lengths)>15: 
        return(15)
    else:
        return(max(Lengths))

#Function to remove charges and unwanted characters in smiles strings
ChargeRemover=re.compile(r'\+|\-')
ExtraPostasRemover=re.compile('^\'|\'$')

def LabelStrippersOne(x):
    
    #Remove labels
    x=str(x).replace('*','')
    
    #Remove extra Rs?
    x=x.replace('R','')
    
    #Remove charges
    x=ChargeRemover.sub('',string=str(x))
    
    #Remove any rextra apostrophes
    x=ExtraPostasRemover.sub('',string=str(x))
    
    return(x)

#Name converter from Smiles to 'English'
def ReactantConvertOne(x):
    try:
        x=PubChemCpdsDict[x]
    except:
        x='nomatch'
    return(x)


def CombinedConvertMetabs(OneLabeledMatrixCpd):
    
    #Get number of label counts
    OneLabeledMatrixCpdCarbs=str(OneLabeledMatrixCpd).count('*')
    
    #Strip labels and charges
    OneLabeledMatrixCpdsEnglish=LabelStrippersOne(OneLabeledMatrixCpd)
    
    #Convert to English using metabolite dictionary
    OneLabeledMatrixCpdsEnglish=ReactantConvertOne(OneLabeledMatrixCpdsEnglish)
    
    #Format string to use both english name and label number
    OneLabeledMatrixCpd=str('{0} M+{1}').format(OneLabeledMatrixCpdsEnglish,OneLabeledMatrixCpdCarbs)
    
    return(OneLabeledMatrixCpd)

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

#Taken from  itertools guide
def islice(iterable, *args):
    s = slice(*args)
    it = iter(range(s.start or 0, s.stop or sys.maxsize, s.step or 1))
    try:
        nexti = next(it)
    except StopIteration:
        return
    for i, element in enumerate(iterable):
        if i == nexti:
            yield element
            nexti = next(it)

#Taken from itertools guide
def roundrobin(*iterables):
    #"roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    # Recipe credited to George Sakkis
    pending = len(iterables)
    nexts = cycle(iter(it).__next__ for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))

#Convert each isotopomer in a row to the common name with the M+x in tact
def IsotopomerToCommonName(OneRowPathMatrix):
    ConvertList=[CombinedConvertMetabs(x) for x in OneRowPathMatrix.values[0]if not type(x)==list] #good
    NonConvertList=[x for x in OneRowPathMatrix.values[0] if type(x)==list]
    OneRowPathMatrix=pandas.DataFrame(list(roundrobin(ConvertList,NonConvertList))).T
    return(OneRowPathMatrix)

#Remove metabolites which cannot be name-converted from smiles
def RemoveTheUnnamed(OnePathMatrix):
    DropIndices=[]
    for rxn in range(len(OnePathMatrix)):
        if any([x for x in OnePathMatrix[rxn:rxn+1].values[0] if 'nomatch' in x]):
            DropIndices.append(rxn)
    OnePathMatrix=OnePathMatrix.drop(DropIndices)
    OnePathMatrix=OnePathMatrix.reset_index(drop=True)
    return(OnePathMatrix)

#Convert path matrix to common name
def ConvertPathMatrixToEnglish(OnePathMatrix):
    NewMatrix=pandas.DataFrame()
    for row in range(len(OnePathMatrix)):
        NewMatrix=NewMatrix.append(IsotopomerToCommonName(OnePathMatrix[row:row+1]))
    NewMatrix=NewMatrix.reset_index(drop=True)
    NewMatrix=RemoveTheUnnamed(NewMatrix)
    return(NewMatrix)

#With a given input isotopologue name, calculate the possible isotopomer combinations and run each through the BuildPathsCoreParallel function
def BuildPathsSetupIsotopomers(Dictionaries,userinputname,PathLength):
    
    #Get the smiles of interest
    CpdMatch=SmilesOptionsFromUserInput(userinputname)
    
    #Get M+x
    IsotopologueNum=re.findall('M\+\d',userinputname)[0]
    
    #Get number of labeled carbons
    NumOfLabels=int(re.findall('\d',IsotopologueNum)[0])
    
    #Generate cominbations of isotopomers to search for
    IsotopomerList=IsotopomerCombinations(CpdMatch,NumOfLabels)
    
    if PathLength=='':
    
        #Get the max path length - note this is currently capped at 15
        StopDict=CalculateIsotopomerLength(Dictionaries,IsotopomerList) #subtract one to correctly index ReusltList - oct 2017
    
    else:
        
        StopDict=PathLength #same as above for the 1 subtraction - oct 2017
    
    #Hard set the length
    #StopDict=11

    for Isotopomer in IsotopomerList:
        try:
            for SubLength in list(reversed(range(StopDict+1))):
                print('Calculating',Isotopomer,'Paths of Length',SubLength) #adding one to index the correct round - oct 2017
                BuildPathsCoreParallel(SubLength,Isotopomer,userinputname,StopDict)
        
        except:
            print('Isotopomer failed')
            continue
                    
    #Given the parallel form of this script, the DF should be read back in
    #.. and duplicate rows should be dropped
    try:
        rxnlength=int(re.findall('(\d+)Rxns', '{0}_Paths_{1}Rxns.csv'.format(userinputname,PathLength))[0])
        TempDF=pandas.read_csv('{0}_Paths_{1}Rxns.csv'.format(userinputname,PathLength),header=None,error_bad_lines=False, names=list(range(rxnlength*2+1)))
        TempDF=TempDF.drop_duplicates()
        TempDF=TempDF.reset_index(drop=True)
        TempDF=DropDuplicatedMetabNamesPostCalc(TempDF)
        TempDF=RemoveDeadEnds(TempDF, rxnlength) #If paths don't start with glucose
        with open('{0}_Paths_{1}Rxns.csv'.format(userinputname,PathLength),'w') as fp: #'w' to overwrite the file
            TempDF.to_csv(fp,index=False,header=False)
        fp.close()
    except:
        pass

#Drop any duplicate name within the same row - this function is to remove futile cycling in metabolic path searches
def DropDuplicatedMetabNamesPostCalc(SeedPathMatrix):
    DropIndices=[]
    for rxn in range(len(SeedPathMatrix)):
        if len(list(set([x for x in SeedPathMatrix.iloc[rxn].values if 'M+' in str(x)]))) != len(list([x for x in SeedPathMatrix.iloc[rxn].values if 'M+' in str(x)])):
            DropIndices.append(rxn)
    SeedPathMatrix=SeedPathMatrix.drop(DropIndices)
    SeedPathMatrix=SeedPathMatrix.reset_index(drop=True)
    return(SeedPathMatrix)

#Paths may have dead ends on the way back to the input tracer, these should be dropped
def RemoveDeadEnds(LabelDF, rxnlength):
    colnames = []
    #rename columns
    for rxn in reversed(range(1, rxnlength+2)):
        colnames.append(['Reaction ' + str(rxn), 'Metabolite ' + str(rxn)])
    colnames = list(chain.from_iterable(colnames))
    colnames = colnames[1:]
    LabelDF.columns = colnames
    #Drop dead end paths/rows
    DropList = []
    MetabOnly = LabelDF[list(x for x in LabelDF.columns if 'Metabolite' in x)] #only check metabolite columns
    for row in range(len(MetabOnly)):
        for metab in reversed(range(len(MetabOnly.iloc[row].values))):
            #Note this only works for the glucose search
            if str(MetabOnly.iloc[row].values[metab]) == 'nan' and MetabOnly.iloc[row].values[metab-1] != tracerinput and str(MetabOnly.iloc[row].values[metab-1]) != 'nan':
                DropList.append(row)
    DropList=list(set(DropList))
    LabelDF=LabelDF.drop(DropList)
    LabelDF=LabelDF.reset_index(drop=True)
    return(LabelDF)

#Build a DF that contains the product of interest and the res
def BuildOneRoundPath(ResultDictionary,MetabName):
    
    #Get around key error problem
    try:
        #Get the reactants and enzymes key'ed to the metabolite of interest
        ReactantList=list(ResultDictionary[MetabName].keys())
        EnzymeList=list(ResultDictionary[MetabName].values())
    except:
        ReactantList=['']
        EnzymeList=['']
    
    #Repeated list of original MetabName to match into a pandas DF
    MetabNameFlat=list(repeat(MetabName,len(ReactantList)))
    
    #Build DF
    Output=pandas.DataFrame([MetabNameFlat,EnzymeList,ReactantList]).T
    return(Output)

#ProductDictionaryResult is the result of one round of BuildOneRoundPath
#Prior round will be the x-1 round of the x'th round which built the ProductDictionaryResult
def AddOnPath(ProductDictionaryResult,PriorRoundDictionary):

    #Make sure the paths exist first, otherwise just exit the function
    if ProductDictionaryResult is None:
        return
    
    else:
    
        #Build a list of new matrices - using last column
        TempDF=[BuildOneRoundPath(PriorRoundDictionary,x) for x in ProductDictionaryResult[ProductDictionaryResult.columns[-1]]]
    
        #Faster here than below?
        #Remove the paths which cannot be connected (which have blank fields)
        TempDF=[x for x in TempDF if not '' in x.values]
    
        #Consider a situaion where the product has no match in the ResultsList Dictionary
        #... we'll decide to either return nothing or continue filling the DF
    
        if len(TempDF)==0:
            Output=None
    
        else:
    
            #Extend the prior rounds' rows to match the new TempDF so column-binding can be performed
            ExtendedProductDictionaryResult=[list(repeat(ProductDictionaryResult.iloc[x][0:len(ProductDictionaryResult.columns)-1],len(TempDF[x]))) for x in range(len(TempDF))]
            ExtendedProductDictionaryResult=list(chain.from_iterable(ExtendedProductDictionaryResult))
        
            if len(ExtendedProductDictionaryResult)==1:
                ExtendedProductConcat=pandas.DataFrame(pandas.concat([x for x in ExtendedProductDictionaryResult])).T
            else:
                ExtendedProductConcat=pandas.concat([x for x in ExtendedProductDictionaryResult],axis=1).T
 
            if len(TempDF)>1:
                TempDFConcat=pandas.concat([x for x in TempDF])
            else:
                TempDFConcat=TempDF[0]
            
            ExtendedProductConcat=ExtendedProductConcat.reset_index(drop=True)
            TempDFConcat=TempDFConcat.reset_index(drop=True)

            Output=pandas.concat([ExtendedProductConcat,TempDFConcat],axis=1)
    
            #Drop rows with blank fields
            #Output=Output[Output[Output.columns[-1]].values!='']
    
            #Sorbose problem.. rename it
            Output=Output.replace(to_replace='Sorbose',value='Sorbose, L- M+6')
    
            return(Output)

#Write a function to convert any metabolite names in an output file to more common names
#Its ugly but it works, lets assume dataframes are not millions of rows
def ConvertNames(PathMatrix):
    NewMatrix=pandas.DataFrame()
    for rxns in range(len(PathMatrix)):
        OneRow=PathMatrix[rxns:rxns+1]
        OneRow=OneRow.reset_index(drop=True)
        for x in OneRow:
            try: 
                carbs=re.findall('M\+\d',OneRow.loc[0][x])[0]
                convertname=MetaToCommon[re.sub(' M\+\d','',OneRow.loc[0][x])]
                metabname=convertname+' '+carbs
                OneRow.loc[0][x]=metabname
            except:
                pass
        NewMatrix=NewMatrix.append(OneRow,ignore_index=True)
    return(NewMatrix)

#Drop any duplicate name within the same row - this function is to remove futile cycling in metabolic path searches
def DropDuplicatedMetabNamesAny(SeedPathMatrix):
    DropIndices=[]
    for rxn in range(len(SeedPathMatrix)):
        if len(list(set([x for x in SeedPathMatrix[rxn:rxn+1].values[0] if '*' in x])))!=len(list([x for x in SeedPathMatrix[rxn:rxn+1].values[0] if '*' in x])): 
            DropIndices.append(rxn)
    SeedPathMatrix=SeedPathMatrix.drop(DropIndices)
    SeedPathMatrix=SeedPathMatrix.reset_index(drop=True)
    return(SeedPathMatrix)

#Start with a function to build the paths from a chosen ending point, then this will be embedded in a larger function that
#starts from each possible ending point (e.g. Round 1-20)

def BuildPathsCoreNonParallel(StopDictNum,SeedPathMatrix,MetabName):
        
    StopDictNum = StopDictNum + 1
    #Build paths going backwards, back to the input tracer
    for rxnround in list(reversed(range(StopDictNum))):
        
        #Add on one round to the existing routes
        SeedPathMatrix=AddOnPath(SeedPathMatrix,ResultsList[rxnround])
        
    #Reset column names to a numbered order to append shorter metabolite routes later
    SeedPathMatrix.columns=range(len(SeedPathMatrix.columns))  
    
    SeedPathMatrix=SeedPathMatrix.reset_index(drop=True)
    
    #Drop any duplicate metab name - added Jan29
    SeedPathMatrix=DropDuplicatedMetabNamesAny(SeedPathMatrix)

    #Start from one reaction earlier - loop from 1 to StopDict to build progressively shorter paths
    for subrxn in range(1,StopDictNum,1):
        StartingSubSet=BuildOneRoundPath(ResultsList[StopDictNum-subrxn],MetabName)
        if len(StartingSubSet)>0:
            SeedPathSub=AddOnPath(StartingSubSet,ResultsList[StopDictNum-subrxn-1])
            for rxnround in list(reversed(range(StopDictNum-subrxn-1))):
                SeedPathSub=AddOnPath(SeedPathSub,ResultsList[rxnround])
                
                #If an extension cannot be made in the short path search
                if len(SeedPathSub)==0:
                    del(SeedPathSub)
                
                else:
                    SeedPathSub=SeedPathSub.reset_index(drop=True)
                    #Drop any duplicate metab name - added Jan29
                    SeedPathSub=DropDuplicatedMetabNamesAny(SeedPathSub)
                
        if 'SeedPathSub' in locals():
            SeedPathSub.columns=range(len(SeedPathSub.columns))
            SeedPathMatrix=SeedPathMatrix.append(SeedPathSub,ignore_index=True)

    #Convert metab names from SMILES
    SeedPathMatrix=ConvertPathMatrixToEnglish(SeedPathMatrix)

    #Convert ugly metabolite names        
    SeedPathMatrix=ConvertNames(SeedPathMatrix)
    
    #Write to file
    #SeedPath.to_csv('{0}_Paths.csv'.format(userinputname))
    return(SeedPathMatrix)

def AddOnPathParallel3(SplitMatrix,StopDictNum):
    #Build paths going backwards, back to input tracer
    #for rxnround in list(reversed(range(StopDictNum-1))):
    for rxnround in list(reversed(range(StopDictNum))):

        #Add on one round to the existing routes
        SplitMatrix=AddOnPath(SplitMatrix,ResultsList[rxnround])
        
        #Function to remove any rows which contain multiple instances of MetabName, avoid metabolic cycles/loops
        SplitMatrix=DropDuplicatedMetabNamesAny(SplitMatrix) #Replaced Jan29

    #Reset column names to a numbered order to append shorter metabolite routes later
    SplitMatrix.columns=range(len(SplitMatrix.columns))
    
    SplitMatrix=SplitMatrix.reset_index(drop=True)
    
    SplitMatrix=DropDuplicatedMetabNamesAny(SplitMatrix)
    
    if len(SplitMatrix)==0:
            
        return(SplitMatrix)
            
    else:
        SplitMatrix.columns=range(len(SplitMatrix.columns))
        
        #Convert from SMILES
        SplitMatrix=ConvertPathMatrixToEnglish(SplitMatrix)
        #Convert ugly names            
        SplitMatrix=ConvertNames(SplitMatrix)
                
        return(SplitMatrix)
    
    #return(SplitMatrix)

#Current
#Input with the particular dictionary round of interest
def BuildPathsCoreParallel(ResultDictionaryNum,metabname,userinputname,StopDictNum):
    
    #The starting set of reactants that led to the product isotopomer of interest
    #StartingSet=BuildOneRoundPath(ResultsList[ResultDictionaryNum],metabname)
    StartingSet=BuildOneRoundPath(ResultsList[ResultDictionaryNum-1],metabname)

    if len(StartingSet)>0:
    
        #Starter path matrix
        #SeedPath=AddOnPath(StartingSet,ResultsList[ResultDictionaryNum-1])
        SeedPath=AddOnPath(StartingSet,ResultsList[ResultDictionaryNum-2])

        if SeedPath is not None:
            try:
                #Go one more path matrix
                #SeedPath=AddOnPath(SeedPath,ResultsList[ResultDictionaryNum-2])
                SeedPath=AddOnPath(SeedPath,ResultsList[ResultDictionaryNum-3])

                #if len(SeedPath)>48:
    
                SeedPathSplit=list(numpy.array_split(SeedPath,32))

                #StopDictRepeat=list(repeat(ResultDictionaryNum-2,32))
                StopDictRepeat=list(repeat(ResultDictionaryNum-3,32))

                pooler=mp.Pool(16)
                        
                try:
                        
                    with open('{0}_Paths_{1}Rxns.csv'.format(userinputname,StopDictNum),'a') as fp: #originally 'w' - but append for looping through shorter path lengths
        
                        for result in pooler.starmap(AddOnPathParallel3,zip(SeedPathSplit,StopDictRepeat)):
                            result.to_csv(fp,index=False,header=False)
                    pooler.close()
                    pooler.join()
                    gc.collect()
                    
                except:
                        
                    try:
                                
                        SeedPathSplit=list(numpy.array_split(SeedPath,16))
                        #StopDictRepeat=list(repeat(ResultDictionaryNum-2,16))
                        StopDictRepeat=list(repeat(ResultDictionaryNum-3,16))

                            #pooler=mp.Pool(8)
                            #with closing(mp.Pool(8)) as pooler:
                        with open('{0}_Paths_{1}Rxns.csv'.format(userinputname,StopDictNum),'a') as fp: #originally 'w' - but append for looping through shorter path lengths
                            for result in pooler.starmap(AddOnPathParallel3,zip(SeedPathSplit,StopDictRepeat)):
                                result.to_csv(fp,index=False,header=False)

                        pooler.close()
                        pooler.join()
                        gc.collect()
                            
                    except:

                        #print('Error 1- Some of these paths cannot be connected, try another isotopologue if no results are written to csv [meaning no paths could be connected]')
                        
                        try:
                                #Try smaller splits
                            SeedPathSplit=list(numpy.array_split(SeedPath,8))
                            #StopDictRepeat=list(repeat(ResultDictionaryNum-2,8))
                            StopDictRepeat=list(repeat(ResultDictionaryNum-3,8))

                                #pooler=mp.Pool(4)
                            with open('{0}_Paths_{1}Rxns.csv'.format(userinputname,StopDictNum),'a') as fp: #originally 'w' - but append for looping through shorter path lengths
                                for result in pooler.starmap(AddOnPathParallel3,zip(SeedPathSplit,StopDictRepeat)):
                                    result.to_csv(fp,index=False,header=False)
                    
                            pooler.close()
                            pooler.join()
                            gc.collect()
                        
                        except:

                            #print('Error 2- Some of these paths cannot be connected, try another isotopologue if no results are written to csv [meaning no paths could be connected]')
                        
                            try:
                                #Try smaller splits
                                SeedPathSplit=list(numpy.array_split(SeedPath,4))
                                #StopDictRepeat=list(repeat(ResultDictionaryNum-2,4))
                                StopDictRepeat=list(repeat(ResultDictionaryNum-3,4))

                                #pooler=mp.Pool(4)
                                with open('{0}_Paths_{1}Rxns.csv'.format(userinputname,StopDictNum),'a') as fp: #originally 'w' - but append for looping through shorter path lengths
                                    for result in pooler.starmap(AddOnPathParallel3,zip(SeedPathSplit,StopDictRepeat)):
                                        result.to_csv(fp,index=False,header=False)
                    
                                pooler.close()
                                pooler.join()
                                gc.collect()
            
                            except:
                
                                try:
                                #Try smaller split
                                    SeedPathSplit=list(numpy.array_split(SeedPath,2))
                                    #StopDictRepeat=list(repeat(ResultDictionaryNum-2,2))
                                    StopDictRepeat=list(repeat(ResultDictionaryNum-3,2))

                            #pooler=mp.Pool(2)
                                    with open('{0}_Paths_{1}Rxns.csv'.format(userinputname,StopDictNum),'a') as fp: #originally 'w' - but append for looping through shorter path lengths
                                        for result in pooler.starmap(AddOnPathParallel3,zip(SeedPathSplit,StopDictRepeat)):
                                            result.to_csv(fp,index=False,header=False)
                                    pooler.close()
                                    pooler.join()
                                    gc.collect()
                    
                                except:
 
                                    #print('Error 3 - Some of these paths cannot be connected, try another isotopologue if no results are written to csv [meaning no paths could be connected]')
                  
                                    try:
                                        #print('No Parallel')
                                        #Output=BuildPathsCoreNonParallel(ResultDictionaryNum-2,SeedPath,metabname)
                                        
                                        ##Output=BuildPathsCoreNonParallel(ResultDictionaryNum-3,SeedPath,metabname)
                                        SeedPathSplit=list(numpy.array_split(SeedPath,1))
                                        StopDictRepeat=list(repeat(ResultDictionaryNum-3,1))

                                        ##with open('{0}_Paths_{1}Rxns.csv'.format(userinputname,StopDictNum),'a') as fp: #originally 'w' - but append for looping through shorter path lengths
                                        ##    Output.to_csv(fp,index=False,header=False)
                                        with open('{0}_Paths_{1}Rxns.csv'.format(userinputname,StopDictNum),'a') as fp: #originally 'w' - but append for looping through shorter path lengths
                                            for result in pooler.starmap(AddOnPathParallel3,zip(SeedPathSplit,StopDictRepeat)):
                                                result.to_csv(fp,index=False,header=False)
                                            
                                        pooler.close()
                                        pooler.join()
                                        gc.collect()
                    
                                    except:
                        
                                        print('Paths cannot be built')
                                        pooler.close()
                                        pooler.join()
                                        gc.collect()
                                    
                                        return
                        
            except:
                print('Isotopomer Failed')
                pooler.close()
                pooler.join()
                gc.collect()
                return

#def BuildPaths(): #User input form
def BuildPaths(userinput,pathlength):
    #userinput=input('Enter isotopologue, or type "help" for lists of common isotopologues (then call BuildPaths again with specified isotopologue of interest): ') #User input form
    
    #pathlength=input('Enter path length to search for (please enter a number)- if nothing is entered, length will be automatically calculated and a max length will be set to 15: ') #User input form
    
    if pathlength!='':
    
        pathlength=int(pathlength)
    
        if pathlength>15:
            pathlength=15
    
    if userinput=='help':
        print(CpdConvert.Common)
    if userinput!='help' and userinput.startswith('Acetyl-CoA'):
        print('The number of possible Acetyl-CoA isotopomers is too great to calculate (e.g. 23 choose 2 for M+2 with a 23-carbon compound), and we recommend manually going into this script and specifying your particular isotpomer of interest')
    if userinput!='help' and userinput.startswith('Acetyl-CoA')==False:
        try:
            
            print('Calculating results... this may take a few minutes')
            
            BuildPathsSetupIsotopomers(ResultsList,userinput,pathlength)
                         
            print('Calculations finished, if paths could be found, they will be stored as a csv in the local directory')
            #quit()
        
        except:
            pass
    

#If you want to run a query yourself, an example would be to call BuildPaths('Serine M+3',12)
#BuildPaths('Serine M+3',9)
#BuildPaths('Glutamine M+2', 14)

