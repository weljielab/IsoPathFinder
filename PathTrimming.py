
# coding: utf-8

################################################################################
#Using calculated isotopomer paths, two main approaches of trimming will be used
#Inputs - Xu2011Gill2015_ECs_Cyclers.csv, Tepper2013_GibbsECs.csv, Serine M+3_Paths_9Rxns.csv, Glutamine M+2_Paths_14Rxns.csv
#Outputs - SerM3_Paths_9RxnsFilter.csv, SerM3_Paths_9RxnsNoFilter.csv, GlnM2_Paths_14RxnsNoFilter.csv, GlnM2_Paths_14RxnsNoFilter.csv

#1) Using circadian EC hits from Xu, 2011 and Gill, 2015, to pull out any rows
#...which contain any ECs, and the number of ECs in a row given will also be added
#
#2) Gibbs energies have been estimated in Tepper, 2013, which will be used to 
#...predict path thermodynamics
#
#Note this is an optional script, and an example of how one may go about trimming
#.. excess results through other means, such as previously published 'omic data
#.. or reaction thermodynamics
#
#This script should be used as a guide, as the trimming procedures can be very
#.. dependent on context and experimental objectives
################################################################################

import numpy as np
import pandas as pd
import os, sys, re, itertools, csv
from itertools import chain, repeat
from collections import defaultdict
import dill as pickle
from pathos.helpers import mp

#Read in a list of significant circadian ECs from Xu et al, 2011 and Gill et al, 2015
SigECList=pd.read_csv('Xu2011Gill2015_ECs_Cyclers.csv',header=None)
ECList=list(SigECList[0])

#Edit names to match the format of the PathSearch results
ECList=['EC-'+x for x in ECList]

#Read in the list of estimated Gibbs energies for all annotated enzymes, from Tepper et al., 2013
GibbsList=pd.read_csv('Tepper2013_GibbsECs.csv',skiprows=1)
GibbsECList=['EC-'+x for x in list(GibbsList['ec number'])]
#Note this will take the last delta G value in a list of duplicated ECs

GibbsListDict=dict(zip(GibbsECList,GibbsList['standart Gibbs energy (kJ/mol)']))

#From the fly studies, drop any paths which don't have circadian enzymes
#NOTE!# This step is optional, it is entirely possible that pathway activities may be cyclic without a circadian enzyme (for more detail, see Thurley et al, PNAS, 2017)

def KeggECTrimmer(RowItem,RxnNumber):
    #Convert the stringy EC item into a list of one/multiple ECs
    ListyECs=[re.sub('\[|\'|\]','',x) for x in re.split('\, ',RowItem)]
    
    #if any item across the row is found in the circadian enzyme list
    if any([y for y in ListyECs if y in KeggECList]):
        RowItem=[y for y in ListyECs if y in KeggECList]
        return(RowItem)
    else:
        return(RxnNumber)

#Find ECs for reporting - if there isn't an EC (say only a Metacyc RXN is left), it can be dropped - of course the 'RXN' format in MetaCyc may represent real reactions, but if it isn't connection to an EC the interpretation is more difficult
def MatchFlyKeggECs(PathMatrix):
    DropIndices=[]
    for rxn in range(len(PathMatrix)):
        ECMatchList=[KeggECTrimmer(z,rxn) for z in PathMatrix[rxn:rxn+1].values[0] if type(z) is str and 'EC' in z]
        if any([x==rxn for x in ECMatchList]):
            DropIndices.append(rxn)
            
    #Store indicies where there isn't a RXN -> EC conversion
    PathMatrix=PathMatrix.drop(DropIndices)
    PathMatrix=PathMatrix.reset_index(drop=True)
    return(PathMatrix)

#Sum up the gibbs energies for a given metabolic path
def SumGibbsAcrossRow(PathMatrix):
    MiniGibbList=[]
    for items in PathMatrix:

        try:
            if any([GibbsListDict[x] for x in GibbsECList if x in re.split('\'|\,',items)]):
                
                MiniGibbList.append(min([GibbsListDict[x] for x in GibbsECList if x in re.split('\'|\,',items)]))
        except:
            pass
        
    PathMatrix['Gibbs Sum']=sum(MiniGibbList)
    return(PathMatrix)

#Match ECs in isotopomer paths to significant EC hits from Xu,2011 and Gill,2015
def CircadianECCounter(PathMatrix):
    
    TempECList=[]
    
    for items in PathMatrix:
        
        try:
            
            if any([x for x in ECList if x in re.split('\'|\,',items)])==True:
                    
                    #Add the EC to the list
                TempECList.extend([x for x in ECList if x in re.split('\'|\,',items)])
                    
        except:
            pass
    
    #Only unique items in the ECList
    TempECList=list(set(TempECList))
    
    PathMatrix['Number of Circadian ECs']=len(TempECList)
    PathMatrix['Circadian ECs']=TempECList
    
    return(PathMatrix)

#Function to find ECs in the PathMatrix which match to Xu,2011 and Gill,2015 circadian hits
#This function is currently very slow for large data frames
#Converted to apply functions to take out loops and still slow, some optimization may be required
def MatchCircECHitsAndGetGibbs(PathMatrix, Call='Filter'):
    
    PathMatrix=PathMatrix.reset_index(drop=True)
    
    #Drop duplicates
    PathMatrix=PathMatrix.drop_duplicates()
    
    PathMatrix=PathMatrix.reset_index(drop=True)
    
    #ECs added
    PathMatrixECs=PathMatrix.apply(CircadianECCounter,axis=1)
    
    #Caculate Gibbs across each path
    PathMatrix=PathMatrix.apply(SumGibbsAcrossRow,axis=1)

    
    PathMatrix['Number of Circadian ECs']=PathMatrixECs['Number of Circadian ECs']
    PathMatrix['Circadian ECs']=PathMatrixECs['Circadian ECs']
    
    #Drop rows which don't match any ECs (optional)
    if Call=='Filter':
        PathMatrixDrop=PathMatrix[PathMatrix['Number of Circadian ECs']!=0]
        PathMatrixDrop=PathMatrixDrop[PathMatrixDrop['Gibbs Sum']<0]
    if Call=='No Filter':
        PathMatrixDrop=PathMatrix
    
    PathMatrixDrop=PathMatrixDrop.reset_index(drop=True)
    
    ##Sort by number of matches
    PathMatrixDrop=PathMatrixDrop.sort_values('Gibbs Sum',ascending=True)
    PathMatrixDrop=PathMatrixDrop.reset_index(drop=True)
    
    #Drop certain unwanted items
    PathMatrixDrop=DropSpecifics(PathMatrixDrop)
    
    PathMatrixDrop=PathMatrixDrop.reset_index(drop=True)
    
    return(PathMatrixDrop)

#Drop certain metabolites that aren't relevant or don't match beyond a PubChem ID
#More features can be added here depending on the organism, experimenters' needs etc..
#This is very slow, may need a new way to do this if many more metabolites need to be added

def DropSpecifics(PathMatrix):
    
    #Regex for certain ECs
    EC1=re.compile('\[\'EC-4.1.2.22\'\]')
    EC2=re.compile('\[\'EC-4.1.2.9\'\]')
    EC3=re.compile('\'EC-4.1.1.39\'')
    
    DropIndex=[]
    for rxn in range(len(PathMatrix)):
        if any([x for x in PathMatrix.iloc[rxn].values if type(x) is str and 'Sorbose' in x]):
            DropIndex.append(rxn)
        if any([x for x in PathMatrix.iloc[rxn].values if type(x) is str and 'AC1L1ARE' in x]):
            DropIndex.append(rxn)
        if any([x for x in PathMatrix.iloc[rxn].values if type(x) is str and 'AC1L19OT' in x]):
            DropIndex.append(rxn)
        if any([x for x in PathMatrix.iloc[rxn].values if type(x) is str and '498-23-7' in x]):
            DropIndex.append(rxn)
        if any([x for x in PathMatrix.iloc[rxn].values if type(x) is str and 'AAAHB' in x]):
            DropIndex.append(rxn)
        if any([x for x in PathMatrix.iloc[rxn].values if type(x) is str and 'ETHYLENE' in x]):
            DropIndex.append(rxn)
        if any([x for x in PathMatrix.iloc[rxn].values if type(x) is str and EC1.search(x)]):
            DropIndex.append(rxn)
        if any([x for x in PathMatrix.iloc[rxn].values if type(x) is str and EC2.search(x)]):
            DropIndex.append(rxn)
        if any([x for x in PathMatrix.iloc[rxn].values if type(x) is str and EC3.search(x)]):
            DropIndex.append(rxn)
    DropIndex=list(set(DropIndex))
    PathMatrix=PathMatrix.drop(DropIndex)
    PathMatrix=PathMatrix.reset_index(drop=True)
    
    return(PathMatrix)

#Cleanup and remove duplicate rows
def CleanAndRewrite(PathDF, rxnlength):
    PathDF=PathDF.sort_values(PathDF.columns[-3],ascending=True)
    PathDF=PathDF.drop_duplicates()
    PathDF=PathDF.reset_index(drop=True)
    #Goal - combine rows
    colnames = []
    
    #rename columns
    for rxn in reversed(range(1, rxnlength+2)):
        colnames.append(['Reaction ' + str(rxn), 'Metabolite ' + str(rxn)])
    colnames = list(chain.from_iterable(colnames)) + ['Gibbs'] + ['Circadian Enzyme Number'] + ['Circadian Enzymes']
    colnames = colnames[1:]
    PathDF.columns = colnames
    
    #Bring it all together: replace that row of the SerDF with the combined ECs
    #Without the reaction columns, check if rows are duplicated
    SubNames = [col for col in PathDF.columns if 'Metabolite' in col or 'Gibbs' in col]
    SubSer = PathDF[SubNames]
    SubSer = SubSer.drop_duplicates()
    #if they are, then you need to combine the ECs in those rows
    dupindex = [x for x in list(PathDF.index) if x not in list(SubSer.index)]
    for dup in dupindex:
        #find max value of an index from SerDF which is less than dupindex
        serindex = (max([x for x in list(PathDF.index) if x < dup]))
        SerMini = PathDF.loc[serindex:dup]
        fillvals = []
        
        for col in SerMini.columns:
            rxnlist = [x for x in SerMini[col]]
            rxnlist = [re.sub('\[|\]|\"', '', str(x)) for x in rxnlist]
            rxnlist = list(set(rxnlist))
            fillvals.append(rxnlist)
            
        temp = pd.DataFrame([list((set(x))) for x in fillvals]).T
        for val in range(len(temp.iloc[1])):
            if temp.iloc[1][val] is not None:
                temp.iloc[0][val] = temp.iloc[0][val] + ", " + temp.iloc[1][val]
        
        temp.columns = PathDF.columns
        SubNames = [col for col in temp.columns if 'Reaction' in col or 'Circadian Enzymes' in col]
        for col in SubNames:
            temp.iloc[0][col] = [x for x in list(temp[col]) if x is not None]
        PathDF.iloc[serindex] = temp.iloc[0]
    
    PathDF = PathDF.drop(dupindex)
    PathDF = PathDF.reset_index(drop = True)
    
    for number in range(len(PathDF)):
        if type(PathDF.iloc[number]['Circadian Enzyme Number'])==str:
            PathDF.iloc[number]['Circadian Enzyme Number'] = max(PathDF.iloc[number]['Circadian Enzyme Number'])
    
    #Clean hideous lists
    for enz in range(len(PathDF)):
        PathDF.iloc[enz]['Circadian Enzymes'] = list(set(re.split('\,', re.sub('\\\|\'|\[|\]|\"| ', '', str(PathDF['Circadian Enzymes'][enz])))))
    
    return(PathDF)

#Read in the results of the isotopologue file of interest
rxnlength=int(re.findall('(\d+)Rxns', 'Serine M+3_Paths_9Rxns.csv')[0])
MetabDF=pd.read_csv('Serine M+3_Paths_9Rxns.csv',header=None,error_bad_lines=False, names=list(range(rxnlength*2+1)))

#Clean up serine
ParaNum = 2
PathMatrixSplit=list(np.array_split(MetabDF,ParaNum)) #this could give an error if the number of rows is <16
    
pooler=mp.Pool(ParaNum)

Call = list(np.repeat('Filter', ParaNum))
#In case it exists already
if 'Trimmed_Paths.csv' in os.listdir():
    os.remove('Trimmed_Paths.csv')

with open('Trimmed_Paths.csv','a') as fp: #originally 'w' - but append for looping through shorter path lengths
    #written out to file as its being built
        #for result in pooler.imap(MatchCircECHitsAndGetGibbsFilter,PathMatrixSplit):
        for result in pooler.starmap(MatchCircECHitsAndGetGibbs,zip(PathMatrixSplit, Call)):
            #Each result is a Pandas Object, so write it to csv
            result.to_csv(fp,index=False,header=False)

pooler.close()
pooler.join()

##Read back in to sort
##Sort by number of matches
PathDF=pd.read_csv('Trimmed_Paths.csv',header=None)
Out=CleanAndRewrite(PathDF, rxnlength)
Out=Out.T
Out.columns=['Path ' + str(x) for x in list(Out.columns+1)]
Out.to_csv('SerM3_Paths_9RxnsFilter.csv')

#No filter, what was lost?
rxnlength=int(re.findall('(\d+)Rxns', 'Serine M+3_Paths_9Rxns.csv')[0])
MetabDF=pd.read_csv('Serine M+3_Paths_9Rxns.csv',header=None,error_bad_lines=False, names=list(range(rxnlength*2+1)))

#Clean up serine
ParaNum = 2
PathMatrixSplit=list(np.array_split(MetabDF,ParaNum)) #this could give an error if the number of rows is <16
    
pooler=mp.Pool(ParaNum)

Call = list(np.repeat('No Filter', ParaNum))
#In case it exists already
if 'Trimmed_Paths.csv' in os.listdir():
    os.remove('Trimmed_Paths.csv')

with open('Trimmed_Paths.csv','a') as fp: #originally 'w' - but append for looping through shorter path lengths
    #written out to file as its being built
        #for result in pooler.imap(MatchCircECHitsAndGetGibbsFilter,PathMatrixSplit):
        for result in pooler.starmap(MatchCircECHitsAndGetGibbs,zip(PathMatrixSplit, Call)):
            #Each result is a Pandas Object, so write it to csv
            result.to_csv(fp,index=False,header=False)

pooler.close()
pooler.join()

##Read back in to sort
##Sort by number of matches
PathDF=pd.read_csv('Trimmed_Paths.csv',header=None)
Out=CleanAndRewrite(PathDF, rxnlength)
Out=Out.T
Out.columns=['Path ' + str(x) for x in list(Out.columns+1)]
Out.to_csv('SerM3_Paths_9RxnsNoFilter.csv')

#Read in the results of the isotopologue file of interest
#MetabDF=pd.read_csv('Glutamine M+2_Paths_14Rxns.csv',header=None,error_bad_lines=False)
rxnlength=int(re.findall('(\d+)Rxns', 'Glutamine M+2_Paths_14Rxns.csv')[0])
MetabDF=pd.read_csv('Glutamine M+2_Paths_14Rxns.csv', header=None, error_bad_lines=False, names=list(range(rxnlength*2+1)))

#Clean up gln
ParaNum = 4
PathMatrixSplit=list(np.array_split(MetabDF,ParaNum)) #this could give an error if the number of rows is <16
    
pooler=mp.Pool(ParaNum)

Call = list(np.repeat('Filter', ParaNum))
#In case it exists already
if 'Trimmed_Paths.csv' in os.listdir():
    os.remove('Trimmed_Paths.csv')

with open('Trimmed_Paths.csv','a') as fp: #originally 'w' - but append for looping through shorter path lengths
    #written out to file as its being built
        #for result in pooler.imap(MatchCircECHitsAndGetGibbsFilter,PathMatrixSplit):
        for result in pooler.starmap(MatchCircECHitsAndGetGibbs,zip(PathMatrixSplit, Call)):
            #Each result is a Pandas Object, so write it to csv
            result.to_csv(fp,index=False,header=False)

pooler.close()
pooler.join()

##Read back in to sort
##Sort by number of matches
PathDF=pd.read_csv('Trimmed_Paths.csv',header=None)
Out=CleanAndRewrite(PathDF, rxnlength)
Out=Out.T
Out.columns=['Path ' + str(x) for x in list(Out.columns+1)]
Out.to_csv('GlnM2_Paths_14RxnsFilter.csv')

#Read in the results of the isotopologue file of interest
#MetabDF=pd.read_csv('Glutamine M+2_Paths_14Rxns.csv',header=None,error_bad_lines=False)
rxnlength=int(re.findall('(\d+)Rxns', 'Glutamine M+2_Paths_14Rxns.csv')[0])
MetabDF=pd.read_csv('Glutamine M+2_Paths_14Rxns.csv', header=None, error_bad_lines=False, names=list(range(rxnlength*2+1)))

#Clean up gln
ParaNum = 4
PathMatrixSplit=list(np.array_split(MetabDF,ParaNum)) #this could give an error if the number of rows is <16
    
pooler=mp.Pool(ParaNum)

Call = list(np.repeat('No Filter', ParaNum))
#In case it exists already
os.remove('Trimmed_Paths.csv')

with open('Trimmed_Paths.csv','a') as fp: #originally 'w' - but append for looping through shorter path lengths
    #written out to file as its being built
        #for result in pooler.imap(MatchCircECHitsAndGetGibbsFilter,PathMatrixSplit):
        for result in pooler.starmap(MatchCircECHitsAndGetGibbs,zip(PathMatrixSplit, Call)):
            #Each result is a Pandas Object, so write it to csv
            result.to_csv(fp,index=False,header=False)

pooler.close()
pooler.join()

##Read back in to sort
##Sort by number of matches
PathDF=pd.read_csv('Trimmed_Paths.csv',header=None)
Out=CleanAndRewrite(PathDF, rxnlength)
Out=Out.T
Out.columns=['Path ' + str(x) for x in list(Out.columns+1)]
Out.to_csv('GlnM2_Paths_14RxnsNoFilter.csv')

