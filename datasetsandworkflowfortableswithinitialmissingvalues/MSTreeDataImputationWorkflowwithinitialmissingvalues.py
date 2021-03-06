#!/usr/bin/python2.7

import pickle
import os
import networkx as nx
import copy
import random as rand
import time
import matplotlib.pyplot as plt
import excelutilities
import FisherExact
from sklearn.linear_model import LinearRegression
import math
import operator
import itertools
import openpyxl as px #############################################################################################################

#note type(2) or type(2).__name are used to make this analysis automaticlally use the default int like object for python

def readinsaveddatastructs(currentdataset,originaldatafolder="init_datasets",pklfilename="saveddatastructs.pkl"):
	"""use this function to read in saved datadic for ms dataset"""
	os.chdir(currentdataset);os.chdir(originaldatafolder);
	savedfile=open(pklfilename);saveddatadic=pickle.load(savedfile);savedfile.close();
	os.chdir("..");os.chdir("..");
	return saveddatadic

#make file to read excelfile
def readexcelandreturntable(filename,sheettoread):
    # remove later
    #move to shared excel utilities
    wkbk=px.load_workbook(filename)# need to load excel package form shared packages
    currentsheet=wkbk.get_sheet_by_name(sheettoread)
    exceltable=[]
    for row in currentsheet.iter_rows():
        excelrow=[]
        for cell in row:
            #print cell, cell.value, type(cell.value),type(cell.value)==type(None)
            if type(cell.value)==type(None):
                break
            elif cell.value==u"X":
                excelrow.append("X")
            elif type(cell.value)!=type(u"2"):
                excelrow.append(int(cell.value))
            else:
                excelrow.append(cell.value)
        exceltable.append(excelrow)
    return exceltable
    
def trimtable(table):
     trimmedtable=[]
     for row in xrange(1,len(table)):
        trimmedtable.append(table[row][1:])
     return trimmedtable

def alltypesthesame(table,typename):
    return sum((type(value)!=typename) for row in table for value in row)==0

def flatten(twodlist):
     return [value for sublist in twodlist for value in sublist]

def countothertypes(tablerow,typename):
    return sum(1 for item in tablerow if (type(item)!=typename))# equivalent to above but faster

def pickfullcol(table):
     pickrows=[]
     savedrows=[]
     ttable=zip(*table)
     for colpos in xrange(len(ttable)):
          if countothertypes(ttable[colpos],type(0).__name__)==0:#note that type(2).__name__ will return int or type of integer for installed python
               pickrows.append(colpos)
               savedrows.append(ttable[colpos])
     return [ list(row) for row in zip(*savedrows)]

def numdiff(l1,l2,nummissing,numcells,nummarkers):
    assert(len(l1)==len(l2))
    return sum(l1[pos]!=l2[pos] for pos in xrange(len(l1)))
     
def magnitudeofdiff(l1,l2,nummissing,numcells,nummarkers):
    assert(len(l1)==len(l2))
    return sum(abs(l1[pos]-l2[pos]) for pos in xrange(len(l1)))

def ratioofdiff(l1,l2,nummissing,numcells,nummarkers):
    assert(len(l1)==len(l2))
    return sum(l1[pos]!=l2[pos] for pos in xrange(len(l1)))/float(nummissing)

def ratioofmagnitudediff(l1,l2,nummissing,numcells,nummarkers):
    assert(len(l1)==len(l2))
    return sum(abs(l1[pos]!=l2[pos]) for pos in xrange(len(l1)))/float(nummissing)
    
def readinexceltablesasdict(saveddatadic,shouldimpute=True):
    imputedic={}
    anytablehasmissingentries=False
    for treepos in xrange(len(saveddatadic["onlyuniquecellsintree"])):
        datasetname=((saveddatadic["treefolders"][treepos]).split("/"))[2]
        imputedic[datasetname]={}
        excelname=saveddatadic["excelname"][treepos]
        sheetname=saveddatadic["excelsheetname"][treepos]
        imputedic[datasetname]["unimputedoriginaltable"]=readexcelandreturntable(excelname,sheetname)
        if not alltypesthesame(trimtable(imputedic[datasetname]["unimputedoriginaltable"]),type(2).__name__):anytablehasmissingentries=True#note that type(2).__name__ will return int or type of integer for installed python
    imputedic["anytablehasmissingentries"]=anytablehasmissingentries
    imputedic["toimpute"]=shouldimpute
    return imputedic

def getlistofmissingpos(impdic):
    dicforlists={}
    for datasetname in impdic:
        if type(impdic[datasetname])==dict:
            dicforlists[datasetname]=[]
            for rowpos in xrange(1,len(impdic[datasetname]["unimputedoriginaltable"])):
                for colpos in xrange(1,len(impdic[datasetname]["unimputedoriginaltable"][0])):
                    if type(impdic[datasetname]["unimputedoriginaltable"][rowpos][colpos])!=type(2):# aka it is not an int
                        dicforlists[datasetname].append((rowpos,colpos))
    return dicforlists

def filtertableformissingpositions(table,missingpositions):
    filteredtable=[]
    for rowpos in xrange(len(table)):
        nextrow=[]
        for colpos in xrange(len(table[0])):
            if (rowpos,colpos) not in missingpositions:
                nextrow.append(table[rowpos][colpos])
        filteredtable.append(nextrow)
    return filteredtable

    
def removeatuniformvaluefrommatrix(table,numrowpos,numcolpos,possiblepositions):
    """ function that takes a 2d table with labeled rows and columns and one at a time blanking out an entry before returning it
according to the paper, Souverein, O. W., Zwinderman, A. H. and Tanck, M. W. T. (2006), Multiple Imputation of Missing Genotype
Data for Unrelated Individuals. Annals of Human Genetics, 70: 372 381. doi: 10.1111/j.1529-8817.2005.00236.x, this method can be
refered to as missing completely at random, which means each entry that is blanked out independently as a subset of the data """
    rowpos,colpos=rand.choice(possiblepositions)
    while(sum(1 for pos in possiblepositions if pos[1] == colpos)<2):rowpos,colpos=rand.choice(possiblepositions)# this ensures that atleast one position is not blanked out in each row allowig 1nn to work fine later
    possiblepositions.remove((rowpos,colpos))# pick a previously unpicked entry to blank out and remove entry as an entry to pick
    table[rowpos][colpos]="X"

def removeatuniformvaluefrommatrixlimitingblanksinanyroworcol(table,numrowpos,numcolpos,possiblepositions):
    """ function that takes a 2d table with labeled rows and columns and one at a time blanking out an entry before returning it
according to the paper, Souverein, O. W., Zwinderman, A. H. and Tanck, M. W. T. (2006), Multiple Imputation of Missing Genotype
Data for Unrelated Individuals. Annals of Human Genetics, 70: 372 381. doi: 10.1111/j.1529-8817.2005.00236.x, this method can be
refered to as not missing at random, which means each entry that is blanked out based on what other entries are already blanked
out, in this case if a row or a col has more than a certain ratio of blanks, no more blanked out in that row or col will be
introduced """
    rowratio,colratio=.9,.9#rowratio,colratio=.5,.5
    rowpos,colpos=rand.choice(possiblepositions)
    rowmissingcount=countothertypes(table[rowpos][1:],type(2))#row col title may be included throws ratio off a bit #note that type(2).__name__ will return int or type of integer for installed python
    colmissingcount=countothertypes(zip(*table)[colpos][1:],type(2))
    if colpos==18:
        print "check",rowpos,colpos,table[rowpos][colpos], table[rowpos][0],table[0][colpos]
        print (rowmissingcount>=int(round(rowratio*numrowpos))), (colmissingcount>=int(round(colratio*numcolpos))),(sum(1 for pos in possiblepositions if pos[1] == colpos)<2)
        print rowmissingcount,rowratio,numrowpos,"    ",colmissingcount,colratio,numcolpos,"    ",sum(1 for pos in possiblepositions if pos[1] == colpos)
        print "out loop"
    while((rowmissingcount>=int(round(rowratio*numrowpos)))or(colmissingcount>=int(round(colratio*numcolpos))) or (sum(1 for pos in possiblepositions if pos[1] == colpos)<2)):
        if colpos==18:print "enter loop:"
        rowpos,colpos=rand.choice(possiblepositions)
        rowmissingcount=countothertypes(table[rowpos][1:],type(2))# put row in a list to make it a 2d list 
        colmissingcount=countothertypes(zip(*table)[colpos][1:],type(2))
        if colpos==18:
            print "check",rowpos,colpos,table[rowpos][colpos], table[rowpos][0],table[0][colpos]
            print (rowmissingcount>=int(round(rowratio*numrowpos))), (colmissingcount>=int(round(colratio*numcolpos))),(sum(1 for pos in possiblepositions if pos[1] == colpos)<2)
            print rowmissingcount,rowratio,numrowpos,"    ",colmissingcount,colratio,numcolpos,"    ",sum(1 for pos in possiblepositions if pos[1] == colpos)
    if colpos==18:print "exit loop",rowpos,colpos,table[rowpos][colpos], table[rowpos][0],table[0][colpos]
    possiblepositions.remove((rowpos,colpos))
    table[rowpos][colpos]="X"

def rauvfmlbiaroc(table,numrow,numcol,otherdata):return removeatuniformvaluefrommatrixlimitingblanksinanyroworcol(table,numrow,numcol,otherdata)

def rauvfm(table,numrow,numcol,otherdata):return removeatuniformvaluefrommatrix(table,numrow,numcol,otherdata)


def producetableswithmissingentries(imputedic,numtables,nummissingentries, picknewposandremovefunc):
    assert(hasattr(nummissingentries,"__iter__"))
    for entry in imputedic:
        if type(imputedic[entry])==dict:
            print("doing entry %s for the method %s\n"%(entry,picknewposandremovefunc.__name__))
            imputedic[entry]["missingentriestablesdata_madeby_%s"%(picknewposandremovefunc.__name__)]={}
            numrow=len(imputedic[entry]["unimputedoriginaltable"])-1 #adjust for fact that first entry in row or col must be a str name
            numcol=len(imputedic[entry]["unimputedoriginaltable"][0])-1
            numrowpos,numcolpos=numcol,numrow
            for row in imputedic[entry]["unimputedoriginaltable"]: print row
            print imputedic[entry]["unimputedoriginaltable"][1][0], numrow, "row", zip(*imputedic[entry]["unimputedoriginaltable"])[0]
            print imputedic[entry]["unimputedoriginaltable"][0][1], numcol, "col", imputedic[entry]["unimputedoriginaltable"][0] 
            initpossiblepositions=[(rpos,colpos) for rpos in xrange(1,numrow) for colpos in xrange(1,numcol)]
            allpossiblepositions=[]
            for pos in initpossiblepositions:
                if type(imputedic[entry]["unimputedoriginaltable"][pos[0]][pos[1]])==int:allpossiblepositions.append(pos)
            print len(allpossiblepositions)
            for nme in nummissingentries:
                imputedic[entry]["missingentriestablesdata_madeby_%s"%(picknewposandremovefunc.__name__)]["tableswith %d additional missingentries"%(nme)]={}
                imputedic[entry]["missingentriestablesdata_madeby_%s"%(picknewposandremovefunc.__name__)]["tableswith %d additional missingentries"%(nme)]["tables"]=[]
                for tablenum in xrange(numtables):
                    nexttable=copy.deepcopy(imputedic[entry]["unimputedoriginaltable"])
                    possiblepositions=copy.deepcopy(allpossiblepositions)
                    for posnum in xrange(nme):picknewposandremovefunc(nexttable,numrowpos,numcolpos,possiblepositions)
                    imputedic[entry]["missingentriestablesdata_madeby_%s"%(picknewposandremovefunc.__name__)]["tableswith %d additional missingentries"%(nme)]["tables"].append(nexttable)
                    for row in zip(*imputedic[entry]["missingentriestablesdata_madeby_%s"%(picknewposandremovefunc.__name__)]["tableswith %d additional missingentries"%(nme)]["tables"][-1])[1:]:
                        #print row, sum(1 for pos in row[1:] if type(pos) == int)
                        #for item in row[1:]: print type(item),item
                        try:assert(sum(1 for pos in row[1:] if type(pos) == int)>0)
                        except:
                                print zip(*imputedic[entry]["missingentriestablesdata_madeby_%s"%(picknewposandremovefunc.__name__)]["tableswith %d additional missingentries"%(nme)]["tables"][-1]).index(row)
                                print "\n"
                                print row, sum(1 for pos in row[1:] if type(pos) == int),picknewposandremovefunc.__name__,initpossiblepositions
                                print possiblepositions,len(possiblepositions),len(initpossiblepositions)
                                for pos in initpossiblepositions:
                                        if pos not in possiblepositions:
                                                print type(imputedic[entry]["missingentriestablesdata_madeby_%s"%(picknewposandremovefunc.__name__)]["tableswith %d additional missingentries"%(nme)]["tables"][-1][pos[0]][pos[1]]),imputedic[entry]["missingentriestablesdata_madeby_%s"%(picknewposandremovefunc.__name__)]["tableswith %d additional missingentries"%(nme)]["tables"][-1][pos[0]][pos[1]],pos


                                print zip(*imputedic[entry]["missingentriestablesdata_madeby_%s"%(picknewposandremovefunc.__name__)]["tableswith %d additional missingentries"%(nme)]["tables"][-1])
                                print zip(*imputedic[entry]["unimputedoriginaltable"])
                                print "\n"
                                print "\n"
                    #print "\n\n"


###Subfunction for knn imputation function#####

def missingvalue(data):
    missing =False
    count=0
    for row in xrange(1,len(data)):
        for col in xrange(1,len(data[0])):
            if type(data[row][col])!=type(2):
                missing =True
                count+=1
    return missing, count

def almostEqual(d1, d2):# dont need almost equal
    #from 15112 cmu website
    #off by .01 percent
    epsilon1 = d1*.01
    epsilon2 = d2*.01
    return (abs(d2 - d1) < min(epsilon1,epsilon2))

def compareintrow(r1,r2):
    assert(len(r1)==len(r2))
    countclosematches=0
    for pos in xrange(len(r1)):
        if (type(r1[pos])==type(2)) and (type(r2[pos])==type(2)):# mabe replaced type.name== with type==
            if r1[pos]==r2[pos]:countclosematches+=1#remove almostequal?
    return countclosematches

def showmissingentryrows(cdata):
    for row in xrange(1,len(cdata)):
        show=False
        for pos in xrange(3,len(cdata[row])):
            if (type(cdata[row][pos])!=type(2)):show=True
        if show: 
            print(cdata[row])
            print("\n"*2)

def addmeanimpvaluesotmatrix(table):
     initimptable=copy.deepcopy(table)
     tinitimptable=[ list(row) for row in zip(*initimptable)]
     collen=len(tinitimptable[0])
     for colpos in xrange(1, len(tinitimptable)):
          data=tinitimptable[colpos][1:]#only get data for row
          filtereddata=filter(lambda a: a!="X",data)#remove non numerical values from row
          average=sum(filtereddata)/float(len(filtereddata))
          for rowpos in xrange(1,collen):
               if tinitimptable[colpos][rowpos]=="X":tinitimptable[colpos][rowpos]=average
     return[ list(row) for row in zip(*tinitimptable)]

def findoptimalsubsetofeaturesusingAIC(table,colpos,filledinblanks,withblanks,numberofentriestoremove=0):
     table=trimtable(table)
     coldata=withblanks#coldata=ttable.pop(colpos)
     usedpositions=range(len(table[0]))#To allow another transpose
     regr= LinearRegression()
     for nextremovedentry in xrange(numberofentriestoremove):
          bestvalue,entrypos=0,None
          for possibleentrytoremove in usedpositions:
               tempsubsetofdata=[list(table[rowpos]) for rowpos in usedpositions if possibleentrytoremove!=rowpos]
               mcoltable=[list(row) for row in zip(*tempsubsetofdata)]
               AICdata= AICscorelinregmodel(len(mcoltable[0]),len(mcoltable),regr,
                                            mcoltable,coldata,filledinblanks)
               if AICdata<bestvalue:bestvalue,entrypos=AICdata,possibleentrytoremove
          usedpositions.remove(entrypos)
     return usedpositions


def AICscorelinregmodel(numparams,numobservation,linregrmod,datacols,truevalues,noblankvalues):
     linregrmod.fit(datacols,noblankvalues[1:])
     impvalues=linregrmod.predict(datacols)
     snoblankvalues=noblankvalues[1:]
     struevalues=truevalues[1:]
     usedtrue,usedimputed=[],[]
     for pos in xrange(len(struevalues)):
          if type(struevalues[pos])!=str:
               usedtrue.append(snoblankvalues[pos])
               usedimputed.append(impvalues[pos])
     assert(len(usedtrue)==len(usedimputed))
     RSS=sum((usedtrue[pos]-usedimputed[pos])**2 for pos in xrange(len(usedtrue)))
     if (RSS==0.0): RSS=0.00000000000001 # prevetn doing log(0)
     numusedforRSS=len(impvalues)
     assert(numusedforRSS>0)
     return 2*(numparams+1) + numusedforRSS*math.log(RSS/float(numusedforRSS))


def makedisttable(table):
     labels=[row[0] for row in table[1:]]
     trimedtable=trimtable(table)
     assert(len(labels)==len(trimedtable))
     disttable=[]
     for firstcellpos in xrange(len(trimedtable)):
          nextrow=[]
          for secondcellpos in xrange(len(trimedtable)):
               fcelldata=trimedtable[firstcellpos]
               scelldata=trimedtable[secondcellpos]
               assert(len(fcelldata)==len(scelldata))
               nextrow.append(sum( abs(fcelldata[pos]-scelldata[pos]) for pos in xrange(len(fcelldata)) if ((type(fcelldata[pos])!=str) and (type(scelldata[pos])!=str))))                         
          disttable.append(nextrow)
     for pos in xrange(len(trimedtable)): assert(disttable[pos][pos]==0)
     return disttable, labels

def doneighborjoiningtreebuilding(initmat,initlabels,nodedist):
    #turn matrix into a dictionary of distances
    D={}
    for speciespos1 in xrange(len(initlabels)):
        D[initlabels[speciespos1]]={}
        for speciespos2 in xrange(len(initlabels)):
            D[initlabels[speciespos1]][initlabels[speciespos2]]=initmat[speciespos1][speciespos2]
    nextlabel=initlabels
    while(len(nextlabel)>1):
        numitemssummed=len(nextlabel)# allother items, 'n'
        u={}
        #step 1
        for speciespos1 in xrange(len(initlabels)):
            distsum=0
            for speciespos2 in xrange(len(initlabels)):
                if speciespos1==speciespos2:continue
                distsum+=D[initlabels[speciespos1]][initlabels[speciespos2]]
            u[initlabels[speciespos1]]=distsum/(numitemssummed-2)

        #step2
        minval=10**5
        minpair=None
        for speciespos1 in xrange(len(initlabels)):
            for speciespos2 in xrange(len(initlabels)):
                if speciespos1==speciespos2:continue
                metric=D[initlabels[speciespos1]][initlabels[speciespos2]]-u[initlabels[speciespos1]]-u[initlabels[speciespos2]]
                if metric<minval:
                    minval=metric
                    minpair=(speciespos1,speciespos2)

        #part of step 5 - delete from lable list
        nextlabel=list(initlabels)
        nextlabel.remove(initlabels[minpair[0]])
        nextlabel.remove(initlabels[minpair[1]])
        nextlabel.append("-".join([initlabels[minpair[0]],initlabels[minpair[1]]]))

        #step 3 join 2 nodes
        nodedist[nextlabel[-1]]={}
        nodedist[nextlabel[-1]][initlabels[minpair[0]]]=(D[initlabels[minpair[0]]][initlabels[minpair[1]]]+u[initlabels[minpair[0]]]-u[initlabels[minpair[1]]])/2
        #nodedist[initlabels[minpair[0]]][nextlabel[-1]]=(D[initlabels[minpair[0]]][initlabels[minpair[1]]]+u[initlabels[minpair[0]]]-u[initlabels[minpair[1]]])/2
        nodedist[nextlabel[-1]][initlabels[minpair[1]]]=(D[initlabels[minpair[0]]][initlabels[minpair[1]]]+u[initlabels[minpair[1]]]-u[initlabels[minpair[0]]])/2        
        #nodedist[initlabels[minpair[1]]][nextlabel[-1]]=(D[initlabels[minpair[0]]][initlabels[minpair[1]]]+u[initlabels[minpair[1]]]-u[initlabels[minpair[0]]])/2        

        #step 4 make distances to new node from all other nodes
        D[nextlabel[-1]]={}
        for speciespos1 in xrange(len(initlabels)):
            if speciespos1 in minpair: continue
            D[nextlabel[-1]][initlabels[speciespos1]]=(D[initlabels[minpair[0]]][initlabels[speciespos1]]+D[initlabels[minpair[1]]][initlabels[speciespos1]]-D[initlabels[minpair[0]]][initlabels[minpair[1]]])/2
            D[initlabels[speciespos1]][nextlabel[-1]]=D[nextlabel[-1]][initlabels[speciespos1]]
            
        #other part of step5- delete from distance table
        #del with dict entry

        for oldsp in minpair:
            del D[initlabels[oldsp]]
            for currentspecies in nextlabel:
                if currentspecies in D:
                    if initlabels[oldsp] in D[currentspecies]:
                        del D[currentspecies][initlabels[oldsp]]
        initlabels=nextlabel


        if len(nextlabel)==2:
            #step 6
            secondlastnode=nextlabel[0]
            lastnode=nextlabel[1]
            nextlabel=["-".join([initlabels[minpair[0]],initlabels[minpair[1]]])]
            nodedist[nextlabel[-1]]={}
            nodedist[nextlabel[-1]][secondlastnode]=D[secondlastnode][lastnode]
            #nodedist[secondlastnode][nextlabel[-1]]=D[secondlastnode][lastnode]

            nodedist[nextlabel[-1]][lastnode]=D[secondlastnode][lastnode]
            #nodedist[lastnode][nextlabel[-1]]=D[secondlastnode][lastnode]
            #nodedist[secondlastnode][lastnode]=D[secondlastnode][lastnode]
            break


    # turn edge data into tree
    searchabletree=nx.Graph()
    for snode in nodedist:
         if snode not in searchabletree.nodes():
              searchabletree.add_node(snode)
         for enode in nodedist[snode]:
              if enode not in searchabletree.nodes():
                   searchabletree.add_node(enode)
              #print nodedist[snode][enode]
              searchabletree.add_edge(snode,enode,weight=nodedist[snode][enode])
    #print searchabletree.nodes(data=True)
    nx.set_node_attributes(searchabletree,"value",{node:"X" for node in searchabletree.nodes()})
    #print searchabletree.nodes(data=True)
    nx.set_node_attributes(searchabletree,"possiblevalue",{node:{} for node in searchabletree.nodes()})
    #print searchabletree.nodes(data=True)
    return nodedist,searchabletree


def searchnodedata(tree,node):
     for data in tree.nodes(data=True):
          if node==data[0]: return data[1]
          
def assignleavevalues(tree,leaves):
     for leaf in leaves:
          value=searchnodedata(tree,leaf)["value"]
          if type(value)==str:
               for nearnode in nx.bfs_edges(tree,leaf):
                    if type(searchnodedata(tree,nearnode[1])["value"])!=str:
                         searchnodedata(tree,leaf)["value"]=searchnodedata(tree,nearnode[1])["value"]
                         break

def recursivenodelabelling(tree,startingnode,parentnode,allposvalue="all"):#,values
     if len(tree.neighbors(startingnode))==1:
          posvaldic=searchnodedata(tree,startingnode)["possiblevalue"]
          value=searchnodedata(tree,startingnode)["value"]
          for val in allposvalue:
               if val!=value:posvaldic[val]=10**10
               else:posvaldic[val]=0
          return startingnode,posvaldic,0 #return the opt c for this leave which by definition is 0
     nodedata=[]
     posvaldic=searchnodedata(tree,startingnode)["possiblevalue"]
     for newnode in tree.neighbors(startingnode):
          if newnode!=parentnode:
               listofvaluesforsubtree=recursivenodelabelling(tree,newnode,startingnode,allposvalue)
               nodedata.append(listofvaluesforsubtree)
     minC=10**10
     for posval in allposvalue:
          Ccostaccumulator=0
          for othernode in nodedata:
               Ccostaccumulator+=min(othernode[2]+1,othernode[1][posval])
          if Ccostaccumulator<minC:minC=Ccostaccumulator
          posvaldic[posval]=Ccostaccumulator
     return startingnode,posvaldic,minC

def assignnodesviatraceback(tree,startingnode,parentnode,parentchar=None):
     currentcharacter=searchnodedata(tree,startingnode)["value"]
     if parentchar==None:
          if type(currentcharacter)==str:
               posvaldic=searchnodedata(tree,startingnode)["possiblevalue"]
               bestval=min(posvaldic,key=posvaldic.get)
               searchnodedata(tree,startingnode)["value"]=bestval
          else:bestval=currentcharacter
     else:
          if type(currentcharacter)==str:
               posvaldic=searchnodedata(tree,startingnode)["possiblevalue"]
               lowestC=min(posvaldic.values())
               bestval=min(posvaldic,key=posvaldic.get)
               if (lowestC+1)>(posvaldic[parentchar]):searchnodedata(tree,startingnode)["value"]=parentchar
               else:searchnodedata(tree,startingnode)["value"]=bestval
          else:bestval=currentcharacter
     for newnode in tree.neighbors(startingnode):
          if newnode!=parentnode:
               assignnodesviatraceback(tree,newnode,startingnode,bestval)

def runFitch(tree,allposvalues):
     startedge=rand.choice(tree.edges())
     #print "Fitch: phase1"
     listofvaluesforsubtree=recursivenodelabelling(tree,startedge[0],startedge[1],allposvalues)#,values)
     #print "phase2"
     listofvaluesforsubtree=recursivenodelabelling(tree,startedge[1],startedge[0],allposvalues)#,values)
     #print "phase3"
     assignnodesviatraceback(tree,startedge[0],startedge[1])
     #print "phase4"
     assignnodesviatraceback(tree,startedge[1],startedge[0])


######################################3333
def useFitchforimputation(data,*arg):
     disttable,labels=makedisttable(data)
     nodedist,searchabletree=doneighborjoiningtreebuilding(disttable,labels,{})
     return assignmissingvaluesusingmodifiedfitch(data,searchabletree,labels)
     
def assignmissingvaluesusingmodifiedfitch_old(table,treegraph,labels): 
     ttable=[list(row) for row in zip(*table)]#for i in ttable: print i
     impttable=[ttable[0]]
     leaves=[x for x in treegraph.nodes_iter() if treegraph.degree(x)==1]
     for rowpos in xrange(1,len(ttable)):
          row=ttable[rowpos]
          if countothertypes(row[1:],type(2))>0:
               copytree=treegraph.copy()#here is where tree for each msl is made
               maxvalueposvalue=0
               minvalueposvalue=10**10
               for pos in xrange(len(labels)):
                    if (type((row[1:])[pos])!=str):
                         if (row[1:])[pos]> maxvalueposvalue:maxvalueposvalue=(row[1:])[pos]
                         if (row[1:])[pos]< minvalueposvalue:minvalueposvalue=(row[1:])[pos]
               assert(len(row)==(len(labels)+1))
               nx.set_node_attributes(copytree,"value",{labels[pos]:row[pos+1] for pos in xrange(len(labels))})
               allposvalues=range(minvalueposvalue,maxvalueposvalue+1)
               assignleavevalues(copytree,leaves)
               runFitch(copytree,allposvalues)
               newrow=[row[0]]
               for pos in xrange(1,len(row)):
                    if type(row[pos])!=str:newrow.append(row[pos])
                    else: newrow.append(searchnodedata(copytree,labels[pos-1])["value"])
          else:newrow=list(row)
          impttable.append(newrow)
     return [list(finrow) for finrow in zip(*impttable)]


def findleafswithmissingvalues(tree,leaves,dic):
     disofposleavevals=dic
     for leaf in leaves:
          value=searchnodedata(tree,leaf)["value"]
          #print leaf
          if type(value)==str:
               disofposleavevals[leaf]=[]
               oldedgedist=None
               for nearnode in nx.bfs_edges(tree,leaf):
                    #print tree.edges(data=True)
                    #print nearnode,nearnode[0],nearnode[1],tree[nearnode[0]][nearnode[1]]["weight"]
                    if type(searchnodedata(tree,nearnode[1])["value"])!=str:
                         
                         if (type(oldedgedist)!=type(None)) and (tree[nearnode[0]][nearnode[1]]["weight"]> oldedgedist):break# note wieght of edge is distance vie edge
                         #print searchnodedata(tree,nearnode[1])["value"]
                         disofposleavevals[leaf].append(searchnodedata(tree,nearnode[1])["value"])
                         if (oldedgedist!=-1): oldedgedist=tree[nearnode[0]][nearnode[1]]["weight"]
                         #print oldedgedist
                         #searchnodedata(tree,leaf)["value"]=searchnodedata(tree,nearnode[1])["value"]
     namesofleaves=disofposleavevals.keys()
     valuesforkeys=[disofposleavevals[key] for key in namesofleaves]
     return namesofleaves

def getdistancestoallotheredges(tree,leaf):
     dicofallpaths=nx.shortest_path(tree,source=leaf)
     lenofallpaths={}
     for othernode in dicofallpaths:
          lenofallpaths[othernode]=0
          for pairposition in xrange(len(dicofallpaths[othernode])-1):
               weightofedge=tree[dicofallpaths[othernode][pairposition]][dicofallpaths[othernode][pairposition+1]]["weight"]
               #if weightofedge<0:lenofallpaths[othernode]+=0# adjustment neccessary if an edge is negative making further away nodes seem closer
               #else:
               lenofallpaths[othernode]+=weightofedge
     orderedlistofpaths=sorted([(lenofallpaths[othern],othern,searchnodedata(tree,othern)["value"]) for othern in lenofallpaths])# each element has len of path name of other node and that nodes value
     return lenofallpaths,orderedlistofpaths

#dtable,lab=makedisttable(exam3);nodedist,searchabletree=doneighborjoiningtreebuilding(dtable,lab,{});copytree=assignmissingvaluesusingmodifiedfitch_old(exam3,searchabletree,lab);a,b,c=findposleavevalues(copytree,lab,{});dicofdist,sorteddistlist=getdistancestoallotheredges(copytree,"A|B|C|D|E|F");findnearestneighborvalues(sorteddistlist);

def findnearestneighborvalues(orderedlistofpaths):
     dists,othernodes,values=zip(*orderedlistofpaths)
     uniquedists=sorted(list(set(dists)))
     #print dists,othernodes,values,uniquedists
     foundnn=False
     while( not foundnn):
          currentlistofvalues=[]
          currentdist=uniquedists.pop(0)
          for item in orderedlistofpaths:
               if item[0]==currentdist:
                    if type(item[2])==int:
                         currentlistofvalues.append(item[2])
          if len(currentlistofvalues)>0:foundnn=True
          #print currentlistofvalues
     return currentdist, currentlistofvalues

def totalweightoftree(tree):
	weight=0
	for edg in tree.edges(data=True):
		#print edg
		dist=edg[2]["weight"]
		#if dist<0:
		#	weight+=0
		#else:
		weight+=edg[2]["weight"]
	return weight

def assignmissingvaluesusingmodifiedfitch(table,treegraph,labels): 
     ttable=[list(row) for row in zip(*table)]#for i in ttable: print i
     impttable=[ttable[0]]
     leaves=[x for x in treegraph.nodes_iter() if treegraph.degree(x)==1]
     #print leaves
     for rowpos in xrange(1,len(ttable)):
          row=ttable[rowpos]
          if countothertypes(row[1:],type(2))>0:
               copytree=treegraph.copy()#here is where tree for each msl is made
               maxvalueposvalue=0
               minvalueposvalue=10**10
               for pos in xrange(len(labels)):
                    if (type((row[1:])[pos])!=str):
                         if (row[1:])[pos]> maxvalueposvalue:maxvalueposvalue=(row[1:])[pos]
                         if (row[1:])[pos]< minvalueposvalue:minvalueposvalue=(row[1:])[pos]
               assert(len(row)==(len(labels)+1))
               nx.set_node_attributes(copytree,"value",{labels[pos]:row[pos+1] for pos in xrange(len(labels))})
               allposvalues=range(minvalueposvalue,maxvalueposvalue+1)
               #print allposvalues
               #print copytree.edges(data=True)
               #print copytree.nodes(data=True)
               #start loop
               nameofleaves=findleafswithmissingvalues(copytree,labels,{})
               listofcorrespondingposvalues=[]
               for leafpos in xrange(len(nameofleaves)):
                    dicofdist,sorteddistlist=getdistancestoallotheredges(copytree,nameofleaves[leafpos])
                    nndist,nnvalues=findnearestneighborvalues(sorteddistlist)# where values from neareast neighbors are found # assumes there are no negative distances
                    listofcorrespondingposvalues.append(nnvalues)
                    #print leafpos,nameofleaves[leafpos]
                    #print sorteddistlist
                    #print nnvalues
                    #print listofcorrespondingposvalues
               lowestcosttree=10**10
               besttree=None
               for posvalueset in itertools.product(*listofcorrespondingposvalues):# get all combinations of possible set of values for missing leafs-cartesian product
                    postree=copytree.copy()
                    nx.set_node_attributes(postree,"value",{nameofleaves[pos]:posvalueset[pos] for pos in xrange(len(nameofleaves))})
                    runFitch(postree,allposvalues)
                    if totalweightoftree(postree)<lowestcosttree:
                            lowestcosttree=totalweightoftree(postree)
                            besttree=postree
               copytree=besttree
               #
               #assert(False)

               newrow=[row[0]]
               for pos in xrange(1,len(row)):
                    if type(row[pos])!=str:newrow.append(row[pos])
                    else: newrow.append(searchnodedata(copytree,labels[pos-1])["value"])
          else:newrow=list(row)
          impttable.append(newrow)
     return [list(finrow) for finrow in zip(*impttable)]


def knnimputation(data,*arg):#k=1,weights=[1]):
    if len(arg)>0:
        assert(type(arg[0])==type(2))
        assert(arg[0]>0)
        assert(type(arg[1])==list)
        for item in arg[1]:assert(type(item)==float)
        assert(sum(arg[1])==1)
        k,weights=arg[0],arg[1]
        assert(len(weights)==k)
    else:
        k,weights=1,[1.0]    
    # make a ordered list of best hits list
    copydata=copy.deepcopy(data)
    # assumes a certain format
    for row in xrange(1,len(copydata)):
        tdata=copy.deepcopy(copydata)
        tdata.pop(0)
        tdata.remove(copydata[row])
        orderedlist=[]
        while(len(orderedlist)<(len(copydata)-2)):
            bestmatchpos=None
            nummatches=0
            for otherrow in xrange(len(tdata)):
                if nummatches<compareintrow(copydata[row][1:],tdata[otherrow][1:]):
                    nummatches=compareintrow(copydata[row][1:],tdata[otherrow][1:])
                    bestmatchpos=otherrow
            if bestmatchpos==None:bestmatchpos=rand.choice(xrange(len(tdata)))# if no entry is picked pick on at unifrom random to minimize bias
            orderedlist.append(tdata.pop(bestmatchpos))
        assert(len(orderedlist)==(len(copydata)-2))
        for pos in xrange(1,len(copydata[row])):
            if (type(copydata[row][pos])!=type(2)):
                closestentries=[]
                rowcounter=0
                for r in xrange(len(orderedlist)):
                    #print "before", orderedlist[r], orderedlist[r][pos]
                    if (type(orderedlist[r][pos])==type(2)):
                        #print "after", orderedlist[r], orderedlist[r][pos]
                        closestentries.append(orderedlist[r][pos])
                        rowcounter+=1
                        #print closestentries,rowcounter
                        if rowcounter==k: break
                newvalue=0
                #print closestentries, weights, copydata[row]
                assert(len(closestentries)==len(weights))
                for i in xrange(k): newvalue+=closestentries[i]*weights[i]
                copydata[row][pos]=(type(2))(round(newvalue))
    showmissingentryrows(copydata)
    return copydata#imputatedata

def linregimputationscikitlearn_iter(table,numberofentriestoremove=0,cutoffforchange=.01,maxiter=100,*arg):
     tablecopy=copy.deepcopy(table)
     initimptable=addmeanimpvaluesotmatrix(tablecopy)
     olditeration=initimptable
     oldttable=zip(*olditeration)
     nextiterationimp=copy.deepcopy(tablecopy)
     ttable=zip(*nextiterationimp)
     maxdiff=10**10
     counter=0
     while(maxdiff>cutoffforchange):
          if counter>=maxiter:break
          counter+=1
          for coltoimputepos in xrange(1,len(ttable)):# rows for nextieration are transposed to cols in ttables
               coltoimpute=ttable[coltoimputepos]#ignore first entry only text
               if countothertypes(coltoimpute[1:],type(2))>0:
                    fullyspecifiedentries=[list(oldttable[row]) for row in xrange(len(oldttable)) if coltoimputepos!=row]#start here
                    fullyspecifiedtarget=oldttable[coltoimputepos]
                    missingvaluesentries=ttable[coltoimputepos]
                    optset=range(len(fullyspecifiedentries)-1)
                    filteredfullyspecifiedentries=[fullyspecifiedentries[0]]
                    for pos in optset:
                         filteredfullyspecifiedentries.append(fullyspecifiedentries[pos+1])
                    filteredfullyspecifiedentries=[list(row) for row in zip(*filteredfullyspecifiedentries)]
                    regr= LinearRegression()
                    regr.fit(trimtable(filteredfullyspecifiedentries),[[item] for item in (fullyspecifiedtarget[1:])])
                    missingvaluestarget=regr.predict(trimtable(filteredfullyspecifiedentries))
                    assert(len(missingvaluestarget)==len(ttable[coltoimputepos][1:]))# make sure that data obtained and the data recorded  are same sixze without text
                    for rowpos in xrange(1,len(nextiterationimp)):
                         if type(nextiterationimp[rowpos][coltoimputepos])==str:
                              nextiterationimp[rowpos][coltoimputepos]=missingvaluestarget[rowpos-1]
          maxdiff=max(map(operator.sub,flatten(trimtable(olditeration)), flatten(trimtable(nextiterationimp))))#finds the maxdifference between matrices
          olditeration=nextiterationimp
          oldttable=zip(*olditeration)
          nextiterationimp=copy.deepcopy(tablecopy)
          ttable=zip(*nextiterationimp)
     finaltable=olditeration#rename final table
     for improwpos in xrange(1,len(finaltable)):
          for impcolpos in xrange(1,len(finaltable[improwpos])):
               if type(finaltable[improwpos][impcolpos])!=type(2):
                    finaltable[improwpos][impcolpos]=(type(2))(round(finaltable[improwpos][impcolpos]))
     return finaltable

def onennimputationuniform(data,Q=None,*arg):
        result=knnimputation(data,1,[1.0])
        if Q!=None:Q.put(result)
        else:return result
def threennimputationuniform(data,Q=None,*arg):
        result= knnimputation(data,3,[1.0/3]*3)
        if Q!=None:Q.put(result)
        else:return result
def fivennimputationuniform(data,Q=None,*arg):
        result= knnimputation(data,5,[1.0/5]*5)
        if Q!=None:Q.put(result)
        else:return result
def fivennimputationincrease(data,Q=None,*arg):
        result= knnimputation(data,5,[.05,.15,.2,.25,.35])
        if Q!=None:Q.put(result)
        else:return result
def fivennimputationdecrease(data,Q=None,*arg):
        result= knnimputation(data,5,[.35,.25,.2,.15,.05])
        if Q!=None:Q.put(result)
        else:return result
def fivennimputationincreasedegenerate(data,Q=None,*arg):
        result= knnimputation(data,5,[0.0,0.0,0.0,0.0,1.0])
        if Q!=None:Q.put(result)
        else:return result

def linearregressionendswithcutoff(data,Q=None,*arg):
        result= linregimputationscikitlearn_iter(data)
        if Q!=None:Q.put(result)
        else:return result
        
def Fitchforimputation(data,Q=None,*arg):
        result=useFitchforimputation(data)
        if Q!=None:Q.put(result)
        else:return result


def OneNN(data,Q=None,*arg):
        return onennimputationuniform(data,Q)

def LinearReg(data,Q=None,*arg):
        return linearregressionendswithcutoff(data,Q)

def FitchImp(data,Q=None,*arg):
        return Fitchforimputation(data,Q)

def imputeandsavetableusingspecifiedimpmethod(imputedic,impfunc,args=None):          
    for entry in imputedic:
        if type(imputedic[entry])==dict:
            print("doing entry %s for the imputation method %s\n"%(entry,impfunc.__name__))
            for numentry in imputedic[entry]:
                if type(imputedic[entry][numentry])==dict:
                    for subentry in imputedic[entry][numentry]:
                        if type(imputedic[entry][numentry][subentry])==dict:
                            imputedic[entry][numentry][subentry][("%s_imputedtables")%(impfunc.__name__)]={}
                            imputedic[entry][numentry][subentry][("%s_imputedtables")%(impfunc.__name__)]["tables"]=[]
                            if (args==None):imputedic[entry][numentry][subentry][("%s_imputedtables")%(impfunc.__name__)]["imputedoriginaltable"]=impfunc(imputedic[entry]["unimputedoriginaltable"])
                            else:imputedic[entry][numentry][subentry][("%s_imputedtables")%(impfunc.__name__)]["imputedoriginaltable"]=impfunc(imputedic[entry]["unimputedoriginaltable"],*args)
                            for table in imputedic[entry][numentry][subentry]["tables"]:
                                #for row in zip(*table):print row
                                #print "\n\n"
                                if (args==None):imputedic[entry][numentry][subentry][("%s_imputedtables")%(impfunc.__name__)]["tables"].append(impfunc(table))
                                else:imputedic[entry][numentry][subentry][("%s_imputedtables")%(impfunc.__name__)]["tables"].append(impfunc(table,*args))


def saveimputedata(currentdataset,imputedic,statdic=None):
    timestr=time.asctime(time.localtime());date= timestr.split(" ");date=" ".join(date[:3]+date[4:]);
    os.chdir(currentdataset)
    workingfolder=("imputed_datasets"+os.sep)
    if workingfolder[:-1] not in os.listdir("."):os.mkdir(workingfolder)
    os.chdir(workingfolder)
    readm="README.txt";readmfile=open(readm,"w");
    readmfile.write(" Generated Blanks and Imputation began on %s\n"%(timestr));
    for oritab in imputedic:# oritab is for each original table
        print oritab, type(imputedic[oritab])
        if type(imputedic[oritab])==dict:
            os.system("sudo chmod 777 .")
            if oritab not in os.listdir("."):os.mkdir((oritab+os.sep))
            os.system("sudo chmod 777 %s"%oritab+os.sep)
            os.chdir((oritab+os.sep))
            excelutilities.turntable2excelfile(imputedic[oritab]["unimputedoriginaltable"],"%s_originaltable.xlsx"%(oritab),"MSdataset")
            for entry in imputedic[oritab]:
                if type(imputedic[oritab][entry])==dict:
                    readmfile.write("Used the following method to pick entries in table to blank out for %s: %s\n"%(oritab,entry))
                    os.system("sudo chmod 777 .")
                    if entry not in os.listdir("."):os.mkdir((entry+os.sep))
                    os.system("sudo chmod 777 %s"%entry+os.sep)
                    os.chdir((entry+os.sep))  
                    for numentry in imputedic[oritab][entry]:
                        if type(imputedic[oritab][entry][numentry])==dict:
                            readmfile.write("for this iteration did %s.\n"%(numentry))
                            os.system("sudo chmod 777 .")
                            if numentry not in os.listdir("."):os.mkdir((numentry+os.sep))
                            #os.system("sudo chmod 777 %s"%numentry+os.sep)
                            os.chdir((numentry+os.sep)) 
                            for missingvaluetablepos in xrange(len(imputedic[oritab][entry][numentry]["tables"])):excelutilities.turntable2excelfile(imputedic[oritab][entry][numentry]["tables"][missingvaluetablepos],"%s_%d.xlsx"%(entry,missingvaluetablepos),"MSdataset")
                            for subentry in imputedic[oritab][entry][numentry]:
                                if type(imputedic[oritab][entry][numentry][subentry])==dict:
                                    readmfile.write("Used the following method to impute back the entries in table for %s: %s\n"%(oritab,subentry))
                                    os.system("sudo chmod 777 .")
                                    if subentry not in os.listdir("."):os.mkdir((subentry+os.sep))
                                    os.system("sudo chmod 777 %s"%subentry+os.sep)
                                    os.chdir((subentry+os.sep))    
                                    for imputedtablepos in xrange(len(imputedic[oritab][entry][numentry][subentry]["tables"])):excelutilities.turntable2excelfile(imputedic[oritab][entry][numentry][subentry]["tables"][imputedtablepos],"%s_imputedby_%s_%d.xlsx"%(entry,(subentry.split("_"))[0],imputedtablepos),"MSdataset")
                                    os.chdir("..")
                                    if "imputedoriginaltable" not in os.listdir("."):os.mkdir(("imputedoriginaltable"+os.sep))
                                    os.system("sudo chmod 777 %s"%"imputedoriginaltable" +os.sep)
                                    os.chdir(("imputedoriginaltable"+os.sep)) 
                                    excelutilities.turntable2excelfile(imputedic[oritab][entry][numentry][subentry]["imputedoriginaltable"],"%s_imputedby_%s_%s.xlsx"%(entry,(subentry.split("_"))[0],"imputedoriginaltable"),"MSdataset")
                                    os.chdir("..")
                            os.chdir("..")
                    os.chdir("..")
            os.chdir("..")
    readmfile.close()
    saveddatastructfile=open("savedimputedatastructs.pkl","w");pickle.dump(imputedic,saveddatastructfile);
    if statdic!=None:pickle.dump(statdic,saveddatastructfile)
    saveddatastructfile.close()
    os.chdir("..");os.chdir("..");


def makedicofallstats(imputedic,listofmetrics,missposdic):
    reftables={}
    for oritab in imputedic:
        if type(imputedic[oritab])==dict:
            #missinfposfortree=missposdic[oritab]
            #originaltable=trimtable(imputedic[oritab]['originaltable'])
            #addnummissing=countothertypes(flatten(trimtable(imputedic[oritab]["unimputedoriginaltable"])),int)#table starts off with missing values get this number to correct total missing count later
            for entry in imputedic[oritab]:
                if type(imputedic[oritab][entry])==dict:
                    for numentry in imputedic[oritab][entry]:
                        if type(imputedic[oritab][entry][numentry])==dict:    
                            for subentry in imputedic[oritab][entry][numentry]:
                                if type(imputedic[oritab][entry][numentry][subentry])==dict:
                                    reftables[subentry]=trimtable(imputedic[oritab][entry][numentry][subentry]["imputedoriginaltable"])
    #print imputedic[oritab][entry],oritab,entry,numentry,subentry
    #addnummissing=countothertypes(flatten(trimtable(imputedic[oritab][entry][numentry][subentry]["imputedoriginaltable"])),int)#table starts off with missing values get this number to correct total missing count later
    statdic={}
    statdic["listofmetrics"]=listofmetrics
    statdic["collecteddata"]={}
    for metric in listofmetrics:
        currentmetric=metric.__name__
        statdic["collecteddata"][currentmetric]={}
        for oritab in imputedic:# oritab is for each original table
            if type(imputedic[oritab])==dict:
                missinfposfortreetable=missposdic[oritab]
                filteredref=filtertableformissingpositions(imputedic[oritab][entry][numentry][subentry]["imputedoriginaltable"],missinfposfortreetable)
                for entry in imputedic[oritab]:
                    if type(imputedic[oritab][entry])==dict:
                        for numentry in imputedic[oritab][entry]:
                            if type(imputedic[oritab][entry][numentry])==dict:    
                                for subentry in imputedic[oritab][entry][numentry]:
                                    if type(imputedic[oritab][entry][numentry][subentry])==dict:
                                        for reftable in reftables:
                                            if (oritab,entry,numentry,subentry,(reftable.split("_")[0]+"_imputedoriginal")) not in statdic["collecteddata"][currentmetric]:statdic["collecteddata"][currentmetric][(oritab,entry,numentry,subentry,reftable.split("_")[0]+"_imputedoriginal")]=[]
                                            for table in imputedic[oritab][entry][numentry][subentry]["tables"]:
                                                filteredtab=filtertableformissingpositions(table,missinfposfortreetable)
                                                statdic["collecteddata"][currentmetric][(oritab,entry,numentry,subentry,reftable.split("_")[0]+"_imputedoriginal")].append(metric(flatten(filteredref),flatten(filteredtab),type(2)((numentry.split(" "))[1]),len(filteredref),len(filteredref[0])))                                            
    return statdic


def turnstatdicsintoplots(statisticdic,reftabletypes,colorder=None,*args):
        metricnames=[met.__name__ for met in statisticdic["listofmetrics"]]
        rawentries=[entry for entry in statisticdic["collecteddata"][metricnames[0]]]
        listofinputargforimputation=[set(inputargs) for inputargs in zip(*rawentries)]
        datasetnames=list(listofinputargforimputation[0])
        blankoutstrs=[("Blanking Method: "+(entry.split("_")[2])) for entry in list(listofinputargforimputation[1])]#center text how?
        blankoutnames=list(listofinputargforimputation[1])
        imputationnames=[("Imp for Data(%s)"%(entry.split("_"))[0]) for entry in list(listofinputargforimputation[3])]####################################################
        #print reftabletypes[0].split("_")[0], reftabletypes[1].split("_")[0],reftabletypes[2].split("_")[0]
        imputationnames=[(imputationnames[pos]+"and Original(%s)"%(reftabletypes[pos].split("_")[0])) for pos in xrange(len(imputationnames))]######################################################
        #print imputationnames
        imputationentrynames=list(listofinputargforimputation[3])
        if type(colorder)==list:
                imputationentrynames=colorder
                print "entered colorder if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                imputationnames=[("Imp for Data(%s) and Original(%s)"%(imputationentrynames[entry].split("_")[0],reftabletypes[entry].split("_")[0])) for entry in xrange(len(imputationentrynames))]####################
                #print imputationnames
        print(imputationentrynames),"imputationentrynames"
        numberofblankoutnames=list(listofinputargforimputation[2])
        numberofblankedoutentries=sorted([(type(2))(strpart) for entry in numberofblankoutnames for strpart in entry.split() if strpart.isdigit()])
        temp=[]
        for item in numberofblankedoutentries:
                for posstring in numberofblankoutnames:
                        if (" "+str(item)+" ") in posstring:temp.append(posstring)
        numberofblankoutnames=temp
        print(numberofblankoutnames)
        print(numberofblankedoutentries)
        rows=blankoutstrs
        cols=imputationnames
        pad = 10 # in points
        allfigs=[]
        #################################
        #selectedmetric=metricnames[0]

        #################################

        for selectedmetric in metricnames:
                maxy,miny=0,10000000
                maxx,minx=0,10000000

                fig, axes = plt.subplots(nrows=len(rows), ncols=len(cols))#(,figsize=(len(cols)*2,len(rows)*4))# how to make bar and whisker plot in subplots

                for ax, col in zip(axes[0], cols):
                    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                                xycoords='axes fraction', textcoords='offset points',
                                size='large', ha='center', va='baseline')#,size=fontsize)

                for ax, row in zip(axes[:,0], rows):
                    ax.annotate(row, xy=(-.22, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                                xycoords=ax.yaxis.label, textcoords='offset points',
                                size='large', ha='center', va='center',rotation=90)#,size=fontsize)

                plt.setp(axes.flat, xlabel='Number of Missing Entries', ylabel='Metric: %s'%(selectedmetric))

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                print("working on col %d and row %d"%(colpos,rowpos))
                                sorteddata=[]
                                currentimputationmethod=imputationentrynames[colpos]
                                currentcurrentblankpickmethod=blankoutnames[rowpos]
                                for nmepos in xrange(len(numberofblankedoutentries)):
                                        print("begin for this nme %d"%(numberofblankedoutentries[nmepos]))
                                        datafornmeimpandblankout=[]
                                        for accumalteddataname in statisticdic["collecteddata"][selectedmetric]:
                                                if (imputationentrynames[colpos] in accumalteddataname) and (blankoutnames[rowpos] in accumalteddataname) and (numberofblankoutnames[nmepos] in accumalteddataname) and (reftabletypes[colpos] in accumalteddataname):
                                                        print "enter",imputationentrynames[colpos],blankoutnames[rowpos],accumalteddataname
                                                        for metricvalue in statisticdic["collecteddata"][selectedmetric][accumalteddataname]: datafornmeimpandblankout.append(metricvalue)
                                        sorteddata.append(datafornmeimpandblankout)
                                        print("data for %s and %s and %d: \n"%(currentimputationmethod,currentcurrentblankpickmethod,numberofblankedoutentries[nmepos]))
                                        print(datafornmeimpandblankout)
                                boxplotdic=axes[rowpos][colpos].boxplot(sorteddata)
                                axes[rowpos][colpos].set_xticklabels(numberofblankedoutentries)
                                for cap in boxplotdic["caps"]:cap.set_linewidth(4)
                                rangey=axes[rowpos][colpos].get_ylim()
                                rangex=axes[rowpos][colpos].get_xlim()
                                if rangey[1]>maxy:maxy=rangey[1]
                                if rangey[0]<miny:miny=rangey[0]
                                if rangex[1]>maxx:maxx=rangex[1]
                                if rangex[0]<minx:minx=rangex[0]

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                axes[rowpos][colpos].set_ylim(miny,maxy)
                                axes[rowpos][colpos].set_xlim(minx,maxx)
                                
                fig.tight_layout(pad=1,h_pad=.5,w_pad=.5)#         
                allfigs.append((fig, axes))

                fig.suptitle("\nSubplots for Combinations of Parameters to select which and how many Entries to Blank, and Imputation method for %s Metric"%selectedmetric, size =16)
                fig.subplots_adjust(top=.8)

                plt.show()

                print(numberofblankoutnames)
                print(sorteddata)
        return allfigs

def turnstatdicsintoplots_strippedaxis(statisticdic,reftabletypes,colorder=None,coltitle=None,*args):
        metricnames=[met.__name__ for met in statisticdic["listofmetrics"]]
        rawentries=[entry for entry in statisticdic["collecteddata"][metricnames[0]]]
        listofinputargforimputation=[set(inputargs) for inputargs in zip(*rawentries)]
        datasetnames=list(listofinputargforimputation[0])
        blankoutstrs=[("Blanking Method: "+(entry.split("_")[2])) for entry in list(listofinputargforimputation[1])]#center text how?
        blankoutnames=list(listofinputargforimputation[1])
        imputationnames=[("Imp for Data(%s)"%(entry.split("_"))[0]) for entry in list(listofinputargforimputation[3])]
        imputationnames=[(imputationnames[pos]+"and Original(%s)"%(reftabletypes[pos].split("_")[0])) for pos in xrange(len(imputationnames))]
        imputationentrynames=list(listofinputargforimputation[3])
        if type(colorder)==list:
                imputationentrynames=colorder
                print "entered colorder if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
                imputationnames=[("Imp for Data(%s) and Original(%s)"%(imputationentrynames[entry].split("_")[0],reftabletypes[entry].split("_")[0])) for entry in xrange(len(imputationentrynames))]####################
        if type(coltitle)==list:
                imputationnames=coltitle
                print "entered coltitle if !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        print(imputationentrynames),"imputationentrynames"
        numberofblankoutnames=list(listofinputargforimputation[2])
        numberofblankedoutentries=sorted([(type(2))(strpart) for entry in numberofblankoutnames for strpart in entry.split() if strpart.isdigit()])
        temp=[]
        for item in numberofblankedoutentries:
                for posstring in numberofblankoutnames:
                        if (" "+str(item)+" ") in posstring:temp.append(posstring)
        numberofblankoutnames=temp
        print(numberofblankoutnames)
        print(numberofblankedoutentries)

        rows=blankoutstrs
        cols=imputationnames

        pad = 10 # in points

        allfigs=[]

        #################################
        #selectedmetric=metricnames[0]

        #################################

        for selectedmetric in metricnames:
                maxy,miny=0,10000000
                maxx,minx=0,10000000

                fig, axes = plt.subplots(nrows=len(rows), ncols=len(cols))#(,figsize=(len(cols)*2,len(rows)*4))# how to make bar and whisker plot in subplots

                for ax, col in zip(axes[0], cols):
                    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                                xycoords='axes fraction', textcoords='offset points',
                                size='large', ha='center', va='baseline')#,size=fontsize)

                plt.setp(axes.flat, xlabel='Number of Missing Entries', ylabel='Metric: %s'%(selectedmetric))

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                print("working on col %d and row %d"%(colpos,rowpos))
                                sorteddata=[]
                                currentimputationmethod=imputationentrynames[colpos]
                                currentcurrentblankpickmethod=blankoutnames[rowpos]
                                for nmepos in xrange(len(numberofblankedoutentries)):
                                        print("begin for this nme %d"%(numberofblankedoutentries[nmepos]))
                                        datafornmeimpandblankout=[]
                                        for accumalteddataname in statisticdic["collecteddata"][selectedmetric]:
                                                if (imputationentrynames[colpos] in accumalteddataname) and (blankoutnames[rowpos] in accumalteddataname) and (numberofblankoutnames[nmepos] in accumalteddataname) and (reftabletypes[colpos] in accumalteddataname):
                                                        print "enter",imputationentrynames[colpos],blankoutnames[rowpos],accumalteddataname
                                                        for metricvalue in statisticdic["collecteddata"][selectedmetric][accumalteddataname]: datafornmeimpandblankout.append(metricvalue)
                                        sorteddata.append(datafornmeimpandblankout)
                                        print("data for %s and %s and %d: \n"%(currentimputationmethod,currentcurrentblankpickmethod,numberofblankedoutentries[nmepos]))
                                        print(datafornmeimpandblankout)
                                #print sorteddata
                                #print"done with group"
                                boxplotdic=axes[rowpos][colpos].boxplot(sorteddata)
                                axes[rowpos][colpos].set_xticklabels(numberofblankedoutentries)
                                for cap in boxplotdic["caps"]:cap.set_linewidth(4)
                                rangey=axes[rowpos][colpos].get_ylim()
                                rangex=axes[rowpos][colpos].get_xlim()
                                if rangey[1]>maxy:maxy=rangey[1]
                                if rangey[0]<miny:miny=rangey[0]
                                if rangex[1]>maxx:maxx=rangex[1]
                                if rangex[0]<minx:minx=rangex[0]

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                axes[rowpos][colpos].set_ylim(miny,maxy)
                                axes[rowpos][colpos].set_xlim(minx,maxx)
                                
                fig.tight_layout(pad=1,h_pad=.5,w_pad=.5)#         
                allfigs.append((fig, axes))

                plt.show()
                print(numberofblankoutnames)
                print(sorteddata)
        return allfigs


###########################################filteringmarkers#################


def readinsavedimpdatastructs(currentdataset,originaldatafolder="imputed_datasets",pklfilename="savedimputedatastructs.pkl"):
	"""use this function to read in saved impute datadic for ms dataset"""
	os.chdir(currentdataset);os.chdir(originaldatafolder);
	savedfile=open(pklfilename);saveddatadic=pickle.load(savedfile);statdic=pickle.load(savedfile);savedfile.close();
	os.chdir("..");os.chdir("..");
	return saveddatadic,statdic

def defaultfm(table):
        return list(table),[1]*(len(table[0])-1)

def basicfilter(table):
        ttable=[list(col) for col in zip(*table)]
        filtered=[ttable[0]]
        allseenvaluesets=[]
        seencounts=[]
        for col in xrange(1,len(ttable)):
                valueset=tuple(ttable[col][1:])
                if valueset in allseenvaluesets:
                        pos=allseenvaluesets.index(valueset)
                        seencounts[pos]=seencounts[pos]+1
                        continue
                else:
                      allseenvaluesets.append(valueset)
                      seencounts.append(1)
                      filtered.append(ttable[col])
        return [list(col) for col in zip(*filtered)],seencounts

def nullfilter(table):
        return [list(row) for row in table],[1]*(len(table[0])-1)

def makefiltereddata(imputedic,filtermethod):
        colcounts={}
        filtereddic={}
        for entry in imputedic:
                if type(imputedic[entry])==dict:
                    filtereddic[entry]={}
                    colcounts[entry]={}   
                    print("doing entry %s for the filter method %s\n"%(entry,filtermethod.__name__))
                    for numentry in imputedic[entry]:
                        if type(imputedic[entry][numentry])==dict:
                            filtereddic[entry][numentry]={}
                            colcounts[entry][numentry]={}
                            for subentry in imputedic[entry][numentry]:
                                if type(imputedic[entry][numentry][subentry])==dict:
                                    filtereddic[entry][numentry][subentry]={}
                                    colcounts[entry][numentry][subentry]={}
                                    for imputetablemethod in imputedic[entry][numentry][subentry]:
                                        if type(imputedic[entry][numentry][subentry][imputetablemethod])==dict:
                                            filtereddic[entry][numentry][subentry][imputetablemethod]={}
                                            colcounts[entry][numentry][subentry][imputetablemethod]={}
                                            assert(type(imputedic[entry][numentry][subentry][imputetablemethod]["tables"])==list)
                                            alldata=[filtermethod(table) for table in imputedic[entry][numentry][subentry][imputetablemethod]["tables"]]
                                            seperatedata=[list(row) for row in zip(*alldata)]
                                            filtereddic[entry][numentry][subentry][imputetablemethod]["tables"]=seperatedata[0]
                                            colcounts[entry][numentry][subentry][imputetablemethod]["weights"]=seperatedata[1]
                                            seperatedataforimputeorignaldata=filtermethod(imputedic[entry][numentry][subentry][imputetablemethod]["imputedoriginaltable"])
                                            filtereddic[entry][numentry][subentry][imputetablemethod]["imputedoriginaltable"]=seperatedataforimputeorignaldata[0]
                                            colcounts[entry][numentry][subentry][imputetablemethod]["imputedoriginaltablecounts"]=seperatedataforimputeorignaldata[1]
                                            pass
                                        elif type(imputedic[entry][numentry][subentry][imputetablemethod])==list:
                                            alldata=[filtermethod(table) for table in imputedic[entry][numentry][subentry][imputetablemethod]]
                                            seperatedata=[list(row) for row in zip(*alldata)]
                                            assert(imputetablemethod=="tables")
                                            filtereddic[entry][numentry][subentry][imputetablemethod]=seperatedata[0]
                                            colcounts[entry][numentry][subentry]["weights"]=seperatedata[1]
                                        else:
                                            assert(False) 
                                            pass
                                else:
                                    assert(False)
                        elif type(imputedic[entry][numentry])==list:
                                assert(numentry=="unimputedoriginaltable")
                                seperatedata=filtermethod(imputedic[entry][numentry])
                                filtereddic[entry][numentry]=seperatedata[0]
                                colcounts[entry]["weights"]=seperatedata[1]
                        else:
                                filtereddic[entry][numentry]=imputedic[entry][numentry]# need no count for this none dict or list
                else:
                        filtereddic[entry]=imputedic[entry]
        return filtereddic,colcounts
        

def savefilterdata(currentdataset,filterdic,colcountsdic):
    timestr=time.asctime(time.localtime());date= timestr.split(" ");date=" ".join(date[:3]+date[4:]);
    os.chdir(currentdataset)
    workingfolder=("filtered_datasets"+os.sep)
    if workingfolder[:-1] not in os.listdir("."):os.mkdir(workingfolder)
    os.chdir(workingfolder)
    readm="README.txt";readmfile=open(readm,"w");
    readmfile.write(" Generated filtering began on %s\n"%(timestr));
    for oritab in filterdic:# oritab is for each original table
        if type(filterdic[oritab])==dict:
            os.system("sudo chmod 777 .")
            if oritab not in os.listdir("."):os.mkdir((oritab+os.sep))
            os.system("sudo chmod 777 %s"%oritab+os.sep)
            os.chdir((oritab+os.sep))
            excelutilities.turntable2excelfile(filterdic[oritab]["unimputedoriginaltable"],"%s_unimputedoriginaltable.xlsx"%(oritab),"MSdataset")
            for entry in filterdic[oritab]:
                if type(filterdic[oritab][entry])==dict:
                    readmfile.write("Used the following method to pick entries in table to blank out for %s: %s\n"%(oritab,entry))
                    os.system("sudo chmod 777 .")
                    if entry not in os.listdir("."):os.mkdir((entry+os.sep))
                    os.system("sudo chmod 777 %s"%entry+os.sep)
                    os.chdir((entry+os.sep))  
                    for numentry in filterdic[oritab][entry]:
                        if type(filterdic[oritab][entry][numentry])==dict:
                            readmfile.write("for this iteration did %s.\n"%(numentry))
                            os.system("sudo chmod 777 .")
                            if numentry not in os.listdir("."):os.mkdir((numentry+os.sep))
                            os.chdir((numentry+os.sep)) 
                            for missingvaluetablepos in xrange(len(filterdic[oritab][entry][numentry]["tables"])):excelutilities.turntable2excelfile(filterdic[oritab][entry][numentry]["tables"][missingvaluetablepos],"%s_%d.xlsx"%(entry,missingvaluetablepos),"MSdataset")
                            for subentry in filterdic[oritab][entry][numentry]:
                                if type(filterdic[oritab][entry][numentry][subentry])==dict:
                                    readmfile.write("Used the following method to impute back the entries in table for %s: %s\n"%(oritab,subentry))
                                    os.system("sudo chmod 777 .")
                                    if subentry not in os.listdir("."):os.mkdir((subentry+os.sep))
                                    os.system("sudo chmod 777 %s"%subentry+os.sep)
                                    os.chdir((subentry+os.sep))    
                                    for imputedtablepos in xrange(len(filterdic[oritab][entry][numentry][subentry]["tables"])):excelutilities.turntable2excelfile(filterdic[oritab][entry][numentry][subentry]["tables"][imputedtablepos],"%s_imputedby_%s_%d.xlsx"%(entry,(subentry.split("_"))[0],imputedtablepos),"MSdataset")
                                    os.chdir("..")
                                    
                                    if "imputedoriginaltable" not in os.listdir("."):os.mkdir(("imputedoriginaltable"+os.sep))
                                    os.system("sudo chmod 777 %s"%"imputedoriginaltable" +os.sep)
                                    os.chdir(("imputedoriginaltable"+os.sep)) 
                                    excelutilities.turntable2excelfile(filterdic[oritab][entry][numentry][subentry]["imputedoriginaltable"],"%s_imputedby_%s_%s.xlsx"%(entry,(subentry.split("_"))[0],"imputedoriginaltable"),"MSdataset")
                                    os.chdir("..")
                            os.chdir("..")
                    os.chdir("..")
            os.chdir("..")
    readmfile.close()
    saveddatastructfile=open("savedfilterdatastructs.pkl","w");pickle.dump(filterdic,saveddatastructfile);pickle.dump(colcountsdic,saveddatastructfile);saveddatastructfile.close()
    os.chdir("..");os.chdir("..");


############################################################################

def readinsavedfilterdatastructs(currentdataset,originaldatafolder="filtered_datasets",pklfilename="savedfilterdatastructs.pkl"):
        """use this function to read in saved impute datadic for ms dataset"""
        os.chdir(currentdataset);os.chdir(originaldatafolder);
        savedfile=open(pklfilename);saveddatadic=pickle.load(savedfile);coldic=pickle.load(savedfile);savedfile.close();
        os.chdir("..");os.chdir("..");
        return saveddatadic,coldic


def getacentralnode(evotree):
        return rand.choice([node for node in evotree if len(evotree.neighbors(node))>1])

def findnewickrepresentationofatreewrapperwithlabelednodes(tree,pickstartnodefunc=getacentralnode):
        """ function that takes in a networkx tree of the coalescent tree model and returns a newick string representation of tree using the recursive findnewickrepresentationofasubtree """
        if not nx.is_tree(tree):return None
        LCA=pickstartnodefunc(tree)#lastcommonancestor
        dicofsedges=nx.dfs_successors(tree,LCA)
        listofsuccessorstrings=[]
        for successor in dicofsedges[LCA]:
                listofsuccessorstrings.append(findnewickrepresentationofasubtreewithlabelednodes(tree,successor,LCA))
        return "("+ ",".join(listofsuccessorstrings) +")"+LCA+";"

def findnewickrepresentationofasubtreewithlabelednodes(tree,currentroot,parentofroot):
        """ function that takes a network x subtree of the coalescent tree model and return a newick string of a subtree """
        if len(nx.neighbors(tree,currentroot))==1:return currentroot
        subtreesuccessors=nx.dfs_successors(tree,currentroot)
        if parentofroot in subtreesuccessors:subtreesuccessors.pop(parentofroot)
        for node in subtreesuccessors:
                if parentofroot in subtreesuccessors[node]:subtreesuccessors[node].remove(parentofroot)
        listofsuccessorstrings=[]
        for successor in subtreesuccessors[currentroot]:
                listofsuccessorstrings.append(findnewickrepresentationofasubtreewithlabelednodes(tree,successor,currentroot))                
        return "("+ ",".join(listofsuccessorstrings) +")"+currentroot

def findnewickrepresentationofatreewrapperwithlabelednodes_leavesonly(tree,pickstartnodefunc=getacentralnode):
        """ function that takes in a networkx tree of the coalescent tree model and returns a newick string representation of tree using the recursive findnewickrepresentationofasubtree """
        if not nx.is_tree(tree):return None
        LCA=pickstartnodefunc(tree)
        dicofsedges=nx.dfs_successors(tree,LCA)
        listofsuccessorstrings=[]
        for successor in dicofsedges[LCA]:
                listofsuccessorstrings.append(findnewickrepresentationofasubtreewithlabelednodes_leavesonly(tree,successor,LCA))
        if len(listofsuccessorstrings)==1:returnstring= listofsuccessorstrings[0]+";"
        else:returnstring= "("+ ",".join(listofsuccessorstrings) +");"
        return returnstring

def findnewickrepresentationofasubtreewithlabelednodes_leavesonly(tree,currentroot,parentofroot):
        """ function that takes a network x subtree of the coalescent tree model and return a newick string of a subtree """
        if len(nx.neighbors(tree,currentroot))==1:return currentroot
        subtreesuccessors=nx.dfs_successors(tree,currentroot)
        if parentofroot in subtreesuccessors:subtreesuccessors.pop(parentofroot)
        for node in subtreesuccessors:
                if parentofroot in subtreesuccessors[node]:subtreesuccessors[node].remove(parentofroot)
        listofsuccessorstrings=[]
        for successor in subtreesuccessors[currentroot]:
                listofsuccessorstrings.append(findnewickrepresentationofasubtreewithlabelednodes_leavesonly(tree,successor,currentroot))                
        if len(listofsuccessorstrings)==1:returnstring= listofsuccessorstrings[0]
        else:returnstring= "("+ ",".join(listofsuccessorstrings) +")"
        return returnstring


def makeweightdisttable(table,weights=None):
     if weights==None: assert False
     labels=[row[0] for row in table[1:]]
     trimedtable=trimtable(table)
     assert(len(labels)==len(trimedtable))
     disttable=[]
     for firstcellpos in xrange(len(trimedtable)):
          nextrow=[]
          for secondcellpos in xrange(len(trimedtable)):
               fcelldata=trimedtable[firstcellpos]
               scelldata=trimedtable[secondcellpos]
               assert(len(fcelldata)==len(scelldata))
               nextrow.append(sum( weights[pos]*abs(fcelldata[pos]-scelldata[pos]) for pos in xrange(len(fcelldata)) if ((type(fcelldata[pos])!=str) and (type(scelldata[pos])!=str))))                         
          disttable.append(nextrow)
     for pos in xrange(len(trimedtable)): assert(disttable[pos][pos]==0)
     return disttable, labels


def runallnjfortree(imptreetable,weights):
    imptabledist=makedisttable(imptreetable,weights)
    imptreedata=doneighborjoiningtreebuilding(imptabledist[0],imptabledist[1],{})
    return imptreedata[1]

def runallweightnjfortree(imptreetable,weights):
    imptabledist=makeweightdisttable(imptreetable,weights)
    imptreedata=doneighborjoiningtreebuilding(imptabledist[0],imptabledist[1],{})
    return imptreedata[1]

def  readinandsafetreecmps(filterdata,weightdic,reftableimptypes,datatabletypes,treebuilder=runallnjfortree,workdirname="Temp_Work"):
    if workdirname not in os.listdir("."):os.mkdir(workdirname)
    else: assert(False)
    refname=workdirname+"/referencetree.txt"
    imptreesname=workdirname+"/imptrees.txt"
    treecmpoutputname=workdirname+"/tcoutput.txt"
    readintreecmpdata={}
    reftables={}
    """for each original table"""
    for entry in filterdata:
        if type(filterdata[entry])==dict:
            readintreecmpdata[entry]={}
            oldoriginaltable=filterdata[entry]["unimputedoriginaltable"]
            oldoriginalweights=weightdic[entry]["weights"]
            #oldoriginaltreedata=treebuilder(oldoriginaltable,oldoriginalweights)
            #readintreecmpdata[entry]["unimputedoriginaltable"]=oldoriginaltreedata
            #reffile=open(refname,"w")
            #reffile.write(findnewickrepresentationofatreewrapperwithlabelednodes_leavesonly(originaltreedata)+"\n")
            #reffile.close()
            for numentry in filterdata[entry]:
                if type(filterdata[entry][numentry])==dict:
                    readintreecmpdata[entry][numentry]={}
                    for subentry in filterdata[entry][numentry]:
                        if type(filterdata[entry][numentry][subentry])==dict:
                            readintreecmpdata[entry][numentry][subentry]={}
                            #originaltable=filterdata[entry]["unimputedoriginaltable"]
                            for imputetablemethod in filterdata[entry][numentry][subentry]:
                                for reftable in reftableimptypes:
                                        reftables[reftable]={}
                                        for imp in filterdata[entry][numentry][subentry]:
                                                print imp
                                                if reftable.split("_")[0] in imp:
                                                        reftables[reftable]["table"]=filterdata[entry][numentry][subentry][imp]["imputedoriginaltable"]
                                                        reftables[reftable]["counts"]=weightdic[entry][numentry][subentry][imp]["imputedoriginaltablecounts"]
                                        print reftables


                                
                                if type(filterdata[entry][numentry][subentry][imputetablemethod])==dict:
                                    currentreftopick=imputetablemethod.split("_")[0];print currentreftopick;
                                    for curpos in xrange(len(datatabletypes)):
                                            if currentreftopick in datatabletypes[curpos]:
                                                    pos=curpos;print pos
                                                    currentreftopick=reftableimptypes[pos];print currentreftopick;
                                    for possible in reftables:
                                            if currentreftopick in possible:reftopick=possible
                                    print reftopick
                                    refdic=reftables[reftopick]
                                    originaltable=refdic["table"]
                                    originalweights=refdic["counts"]
                                    originaltreedata=treebuilder(originaltable,originalweights)
                                    readintreecmpdata[entry][numentry][subentry][imputetablemethod]={}
                                    readintreecmpdata[entry][numentry][subentry][imputetablemethod]["unimputedoriginaltable"]=originaltreedata
                                    readintreecmpdata[entry][numentry][subentry][imputetablemethod]["unimputedoriginalimpmethod"]=reftopick
                                    reffile=open(refname,"w")
                                    reffile.write(findnewickrepresentationofatreewrapperwithlabelednodes_leavesonly(originaltreedata)+"\n")
                                    reffile.close()                                    
                                    assert(type(filterdata[entry][numentry][subentry][imputetablemethod]["tables"])==list)
                                    imptrees=open(imptreesname,"w")
                                    alltrees=[]
                                    for pos in xrange(len(filterdata[entry][numentry][subentry][imputetablemethod]["tables"])):
                                        imptreetable=filterdata[entry][numentry][subentry][imputetablemethod]["tables"][pos]
                                        impweights=weightdic[entry][numentry][subentry][imputetablemethod]["weights"][pos]
                                        imptreedata=treebuilder(imptreetable,impweights)
                                        alltrees.append(imptreedata)
                                        imptrees.write(findnewickrepresentationofatreewrapperwithlabelednodes_leavesonly(imptreedata)+"\n")
                                    imptrees.close()
                                    sysstr="java -jar TreeCmp/bin/TreeCmp.jar -r %s -d rf ms pd qt -i %s -o %s -P -N"%(refname,imptreesname,treecmpoutputname)
                                    output=os.system(sysstr)
                                    treecmpoutput=open(treecmpoutputname)
                                    extractedata=[ i.split("\t") for i in (treecmpoutput.read()).split("\n")][:-1] # remove last empty entry
                                    treecmpoutput.close()
                                    titles=extractedata[0]
                                    extractedata=[list(row) for row in zip(*extractedata)]

                                    readintreecmpdata[entry][numentry][subentry][imputetablemethod]["data"]=extractedata
                                    readintreecmpdata[entry][numentry][subentry][imputetablemethod]["trees"]=alltrees
            os.remove(refname)
    os.remove(imptreesname)# it gets overwritten in each loop so you only need to remove it once at the end
    os.remove(treecmpoutputname)
    os.rmdir(workdirname)
    print(os.getcwd(),os.listdir("."))
    return readintreecmpdata,titles

def makedicofallstatsfortrees(treedic,listofimpmetrics,titles,startpos,stoppos):
    dataseperateddic={}
    selectedmetrics=titles[startpos:stoppos]
    for pos in xrange(len(selectedmetrics)):
        dataseperateddic[selectedmetrics[pos]]={}
        for entry in treedic:
            for numentry in treedic[entry]:
                if (type(treedic[entry][numentry])==dict):
                    for subentry in treedic[entry][numentry]:
                        for impmethod in treedic[entry][numentry][subentry]:
                            assert(type(treedic[entry][numentry][subentry])==dict)
                            assert(type(treedic[entry][numentry][subentry][impmethod])==dict)
                            assert(type(treedic[entry][numentry][subentry][impmethod]["data"])==list)
                            assert(type(treedic[entry][numentry][subentry][impmethod]["unimputedoriginalimpmethod"])==str)
                            origimpmethod=treedic[entry][numentry][subentry][impmethod]["unimputedoriginalimpmethod"]
                            if (entry,numentry,subentry,impmethod,origimpmethod) not in dataseperateddic[selectedmetrics[pos]]:
                                               dataseperateddic[selectedmetrics[pos]][(entry,numentry,subentry,impmethod,origimpmethod)]=[]
                            dataseperateddic[selectedmetrics[pos]][(entry,numentry,subentry,impmethod,origimpmethod)]=map(float,treedic[entry][numentry][subentry][impmethod]["data"][pos+startpos][1:])
    return dataseperateddic

def turntreestatdicsintoplots(statisticdic,selectedmetric,dataandoriginalimptypes,colorder=None,*args):
        metricnames=[selectedmetric]
        rawentries=[entry for entry in statisticdic[selectedmetric]]
        listofinputargforimputation=[set(inputargs) for inputargs in zip(*rawentries)]
        datasetnames=list(listofinputargforimputation[0])
        blankoutstrs=[("Blanking Method: "+(entry.split("_")[2])) for entry in list(listofinputargforimputation[1])]#center text how?
        blankoutnames=list(listofinputargforimputation[1])
        imputationnames=[("Imp for Data(%s)"%(entry.split("_"))[0]) for entry in list(listofinputargforimputation[3])]
        impfororitables=[]
        for dataimpmeth in list(listofinputargforimputation[3]):
                for entry in dataandoriginalimptypes:
                        if dataimpmeth == entry[1]:impfororitables.append(entry[0])
        print impfororitables
        imputationnames=[(imputationnames[pos]+"and Original(%s)"%(impfororitables[pos].split("_")[0])) for pos in xrange(len(imputationnames))]######################################################
        imputationentrynames=list(listofinputargforimputation[3])
        if type(colorder)==list:
                imputationentrynames=colorder
                impfororitables=[]
                for dataimpmeth in imputationentrynames:
                        for entry in dataandoriginalimptypes:
                                if dataimpmeth == entry[1]:impfororitables.append(entry[0])
                print impfororitables
                imputationnames=[("Imp for Data(%s) and Original(%s)"%(imputationentrynames[entry].split("_")[0],impfororitables[entry].split("_")[0])) for entry in xrange(len(imputationentrynames))]####################
        numberofblankoutnames=list(listofinputargforimputation[2])
        numberofblankedoutentries=sorted([(type(2))(strpart) for entry in numberofblankoutnames for strpart in entry.split() if strpart.isdigit()])
        temp=[]
        for item in numberofblankedoutentries:
                for posstring in numberofblankoutnames:
                        if (" "+str(item)+" ") in posstring:temp.append(posstring)
        numberofblankoutnames=temp
        print(numberofblankoutnames)
        print(numberofblankedoutentries)

        rows=blankoutstrs
        cols=imputationnames

        pad = 10 # in points

        allfigs=[]

        #################################
        #selectedmetric=metricnames[0]

        #################################

        for selectedmetric in metricnames:

                maxy,miny=0,10000000
                maxx,minx=0,10000000

                fig, axes = plt.subplots(nrows=len(rows), ncols=len(cols))#(,figsize=(len(cols)*2,len(rows)*4))# how to make bar and whisker plot in subplots

                for ax, col in zip(axes[0], cols):
                    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                                xycoords='axes fraction', textcoords='offset points',
                                size='large', ha='center', va='baseline')#,size=fontsize)

                for ax, row in zip(axes[:,0], rows):
                    ax.annotate(row, xy=(-.22, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                                xycoords=ax.yaxis.label, textcoords='offset points',
                                size='large', ha='center', va='center',rotation=90)#,size=fontsize)

                plt.setp(axes.flat, xlabel='Number of Missing Entries', ylabel='Metric: %s'%(selectedmetric))

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                print("working on col %d and row %d"%(colpos,rowpos))
                                sorteddata=[]
                                currentimputationmethod=imputationentrynames[colpos]
                                currentcurrentblankpickmethod=blankoutnames[rowpos]
                                for nmepos in xrange(len(numberofblankedoutentries)):
                                        print("begin for this nme %d"%(numberofblankedoutentries[nmepos]))
                                        datafornmeimpandblankout=[]
                                        for accumalteddataname in statisticdic[selectedmetric]:
                                                if (imputationentrynames[colpos] in accumalteddataname) and (blankoutnames[rowpos] in accumalteddataname) and (numberofblankoutnames[nmepos] in accumalteddataname) and (impfororitables[colpos] in accumalteddataname):
                                                        for metricvalue in statisticdic[selectedmetric][accumalteddataname]: datafornmeimpandblankout.append(metricvalue)
                                        sorteddata.append(datafornmeimpandblankout)
                                        print("data for %s and %s and %d: \n"%(currentimputationmethod,currentcurrentblankpickmethod,numberofblankedoutentries[nmepos]))
                                        print(datafornmeimpandblankout)
                                boxplotdic=axes[rowpos][colpos].boxplot(sorteddata)
                                axes[rowpos][colpos].set_xticklabels(numberofblankedoutentries)
                                for cap in boxplotdic["caps"]:cap.set_linewidth(4)
                                rangey=axes[rowpos][colpos].get_ylim()
                                rangex=axes[rowpos][colpos].get_xlim()
                                if rangey[1]>maxy:maxy=rangey[1]
                                if rangey[0]<miny:miny=rangey[0]
                                if rangex[1]>maxx:maxx=rangex[1]
                                if rangex[0]<minx:minx=rangex[0]

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                axes[rowpos][colpos].set_ylim(miny,maxy)
                                axes[rowpos][colpos].set_xlim(minx,maxx)
                                
                                
                fig.tight_layout(pad=1,h_pad=.5,w_pad=.5)#  
                allfigs.append((fig, axes))

                fig.suptitle("\nSubplots for Combinations of Parameters to select which and how many Entries to Blank, and Imputation method for %s Metric"%selectedmetric, size =16)
                fig.subplots_adjust(top=.8)

                plt.show()
                print(numberofblankoutnames)
                print(sorteddata)
        return allfigs

def turntreestatdicsintoplots_strippedaxis(statisticdic,selectedmetric,dataandoriginalimptypes,colorder=None,*args):
        metricnames=[selectedmetric]
        rawentries=[entry for entry in statisticdic[selectedmetric]]
        listofinputargforimputation=[set(inputargs) for inputargs in zip(*rawentries)]
        datasetnames=list(listofinputargforimputation[0])
        blankoutstrs=[("Blanking Method: "+(entry.split("_")[2])) for entry in list(listofinputargforimputation[1])]#center text how?
        blankoutnames=list(listofinputargforimputation[1])
        imputationnames=[("Imp for Data(%s)"%(entry.split("_"))[0]) for entry in list(listofinputargforimputation[3])]
        impfororitables=[]
        for dataimpmeth in list(listofinputargforimputation[3]):
                for entry in dataandoriginalimptypes:
                        if dataimpmeth == entry[1]:impfororitables.append(entry[0])
        print impfororitables
        imputationnames=[(imputationnames[pos]+"and Original(%s)"%(impfororitables[pos].split("_")[0])) for pos in xrange(len(imputationnames))]
        imputationentrynames=list(listofinputargforimputation[3])
        if type(colorder)==list:
                imputationentrynames=colorder
                impfororitables=[]
                for dataimpmeth in imputationentrynames:
                        for entry in dataandoriginalimptypes:
                                if dataimpmeth == entry[1]:impfororitables.append(entry[0])
                print impfororitables
                imputationnames=[("Imp for Data(%s) and Original(%s)"%(imputationentrynames[entry].split("_")[0],impfororitables[entry].split("_")[0])) for entry in xrange(len(imputationentrynames))]
        numberofblankoutnames=list(listofinputargforimputation[2])
        numberofblankedoutentries=sorted([(type(2))(strpart) for entry in numberofblankoutnames for strpart in entry.split() if strpart.isdigit()])
        temp=[]
        for item in numberofblankedoutentries:
                for posstring in numberofblankoutnames:
                        if (" "+str(item)+" ") in posstring:temp.append(posstring)
        numberofblankoutnames=temp
        print(numberofblankoutnames)
        print(numberofblankedoutentries)

        rows=blankoutstrs
        cols=imputationnames

        pad = 10 # in points

        allfigs=[]

        #################################
        #selectedmetric=metricnames[0]

        #################################


        for selectedmetric in metricnames:

                maxy,miny=0,10000000
                maxx,minx=0,10000000

                fig, axes = plt.subplots(nrows=len(rows), ncols=len(cols))#(,figsize=(len(cols)*2,len(rows)*4))# how to make bar and whisker plot in subplots

                for ax, col in zip(axes[0], cols):
                    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                                xycoords='axes fraction', textcoords='offset points',
                                size='large', ha='center', va='baseline')#,size=fontsize)

                plt.setp(axes.flat, xlabel='Number of Missing Entries', ylabel='Metric: %s'%(selectedmetric))

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                print("working on col %d and row %d"%(colpos,rowpos))
                                sorteddata=[]
                                currentimputationmethod=imputationentrynames[colpos]
                                currentcurrentblankpickmethod=blankoutnames[rowpos]
                                for nmepos in xrange(len(numberofblankedoutentries)):
                                        print("begin for this nme %d"%(numberofblankedoutentries[nmepos]))
                                        datafornmeimpandblankout=[]
                                        for accumalteddataname in statisticdic[selectedmetric]:
                                                if (imputationentrynames[colpos] in accumalteddataname) and (blankoutnames[rowpos] in accumalteddataname) and (numberofblankoutnames[nmepos] in accumalteddataname) and (impfororitables[colpos] in accumalteddataname):
                                                        for metricvalue in statisticdic[selectedmetric][accumalteddataname]: datafornmeimpandblankout.append(metricvalue)
                                        sorteddata.append(datafornmeimpandblankout)
                                        print("data for %s and %s and %d: \n"%(currentimputationmethod,currentcurrentblankpickmethod,numberofblankedoutentries[nmepos]))
                                        print(datafornmeimpandblankout)
                                boxplotdic=axes[rowpos][colpos].boxplot(sorteddata)
                                axes[rowpos][colpos].set_xticklabels(numberofblankedoutentries)
                                for cap in boxplotdic["caps"]:cap.set_linewidth(4)
                                rangey=axes[rowpos][colpos].get_ylim()
                                rangex=axes[rowpos][colpos].get_xlim()
                                if rangey[1]>maxy:maxy=rangey[1]
                                if rangey[0]<miny:miny=rangey[0]
                                if rangex[1]>maxx:maxx=rangex[1]
                                if rangex[0]<minx:minx=rangex[0]

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                axes[rowpos][colpos].set_ylim(miny,maxy)
                                axes[rowpos][colpos].set_xlim(minx,maxx)
                                
                                
                fig.tight_layout(pad=1,h_pad=.5,w_pad=.5)#  
                allfigs.append((fig, axes))

                plt.show()

                print(numberofblankoutnames)
                print(sorteddata)
        return allfigs


def turntreestatdicsintoplots_strippedaxis_meanvaluesscatterplotxlog(statisticdic,selectedmetric,dataandoriginalimptypes,colorder=None,hardyrange=None,coltitle=None,*args):
        metricnames=[selectedmetric]
        rawentries=[entry for entry in statisticdic[selectedmetric]]
        listofinputargforimputation=[set(inputargs) for inputargs in zip(*rawentries)]
        datasetnames=list(listofinputargforimputation[0])
        blankoutstrs=[("Blanking Method: "+(entry.split("_")[2])) for entry in list(listofinputargforimputation[1])]#center text how?
        blankoutnames=list(listofinputargforimputation[1])
        imputationnames=[("Imp for Data(%s)"%(entry.split("_"))[0]) for entry in list(listofinputargforimputation[3])]
        impfororitables=[]
        for dataimpmeth in list(listofinputargforimputation[3]):
                for entry in dataandoriginalimptypes:
                        if dataimpmeth == entry[1]:impfororitables.append(entry[0])
        print impfororitables
        imputationnames=[(imputationnames[pos]+"and Original(%s)"%(impfororitables[pos].split("_")[0])) for pos in xrange(len(imputationnames))]
        imputationentrynames=list(listofinputargforimputation[3])
        if type(colorder)==list:
                imputationentrynames=colorder
                impfororitables=[]
                for dataimpmeth in imputationentrynames:
                        for entry in dataandoriginalimptypes:
                                if dataimpmeth == entry[1]:impfororitables.append(entry[0])
                print impfororitables
                imputationnames=[("Imp for Data(%s) and Original(%s)"%(imputationentrynames[entry].split("_")[0],impfororitables[entry].split("_")[0])) for entry in xrange(len(imputationentrynames))]
        if type(coltitle)==list:
                imputationnames=coltitle
        numberofblankoutnames=list(listofinputargforimputation[2])
        numberofblankedoutentries=sorted([(type(2))(strpart) for entry in numberofblankoutnames for strpart in entry.split() if strpart.isdigit()])
        temp=[]
        for item in numberofblankedoutentries:
                for posstring in numberofblankoutnames:
                        if (" "+str(item)+" ") in posstring:temp.append(posstring)
        numberofblankoutnames=temp
        print(numberofblankoutnames)
        print(numberofblankedoutentries)

        rows=blankoutstrs
        cols=imputationnames

        pad = 10 # in points

        allfigs=[]


        #################################
        #selectedmetric=metricnames[0]

        #################################


        for selectedmetric in metricnames:

                maxy,miny=0,10000000
                maxx,minx=0,10000000

                fig, axes = plt.subplots(nrows=len(rows), ncols=len(cols))#(,figsize=(len(cols)*2,len(rows)*4))# how to make bar and whisker plot in subplots

                for ax, col in zip(axes[0], cols):
                    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                                xycoords='axes fraction', textcoords='offset points',
                                size='large', ha='center', va='baseline')#,size=fontsize)


                plt.setp(axes.flat, xlabel='Number of Missing Entries', ylabel='Metric: %s'%(selectedmetric))

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                print("working on col %d and row %d"%(colpos,rowpos))
                                sorteddata=[]
                                currentimputationmethod=imputationentrynames[colpos]
                                currentcurrentblankpickmethod=blankoutnames[rowpos]
                                for nmepos in xrange(len(numberofblankedoutentries)):
                                        print("begin for this nme %d"%(numberofblankedoutentries[nmepos]))
                                        datafornmeimpandblankout=[]
                                        for accumalteddataname in statisticdic[selectedmetric]:
                                                if (imputationentrynames[colpos] in accumalteddataname) and (blankoutnames[rowpos] in accumalteddataname) and (numberofblankoutnames[nmepos] in accumalteddataname)  and (impfororitables[colpos] in accumalteddataname):
                                                        for metricvalue in statisticdic[selectedmetric][accumalteddataname]: datafornmeimpandblankout.append(metricvalue)
                                        sorteddata.append(datafornmeimpandblankout)
                                        print("data for %s and %s and %d: \n"%(currentimputationmethod,currentcurrentblankpickmethod,numberofblankedoutentries[nmepos]))
                                        print(datafornmeimpandblankout)
                                boxplotdic=axes[rowpos][colpos].semilogx(numberofblankedoutentries, [sum(data)/len(data) for data in sorteddata])
                                rangey=axes[rowpos][colpos].get_ylim()
                                rangex=axes[rowpos][colpos].get_xlim()
                                if rangey[1]>maxy:maxy=rangey[1]
                                if rangey[0]<miny:miny=rangey[0]
                                if rangex[1]>maxx:maxx=rangex[1]
                                if rangex[0]<minx:minx=rangex[0]

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                if type(hardyrange)==list:axes[rowpos][colpos].set_ylim(hardyrange[0],hardyrange[1])
                                else:axes[rowpos][colpos].set_ylim(miny,maxy)
                                axes[rowpos][colpos].set_xlim(minx,maxx)
                                
                                
                fig.tight_layout(pad=1,h_pad=.5,w_pad=.5)#  
                allfigs.append((fig, axes))

                plt.show()

                print (numberofblankoutnames)
                print (sorteddata)
        return allfigs


def savetreedata(currentdataset,treedic,treestatdic):
    timestr=time.asctime(time.localtime());date= timestr.split(" ");date=" ".join(date[:3]+date[4:]);
    os.chdir(currentdataset)
    workingfolder=("tree_datasets"+os.sep)
    if workingfolder[:-1] not in os.listdir("."):os.mkdir(workingfolder)
    os.chdir(workingfolder)
    readm="README.txt";readmfile=open(readm,"w");
    readmfile.write(" Generated Trees and data on %s\n"%(timestr));
    readmfile.close()
    saveddatastructfile=open("savedtreedatastructs.pkl","w");pickle.dump(treedic,saveddatastructfile);pickle.dump(treestatdic,saveddatastructfile);saveddatastructfile.close()
    os.chdir("..");os.chdir("..");


def turndatainbins(datalists,minval=0,maxval=1):
     alldata=[]
     for datalist in datalists:alldata.extend(datalist)
     sorteddata=sorted(list(set(alldata)))
     mindif=min(sorteddata[pos]-sorteddata[pos-1]  for pos in xrange(1,len(sorteddata)))
     print(mindif)
     bins=[(minval+pos*mindif) for pos in xrange((type(2))((maxval-minval)/mindif)+1)]
     print(bins) 
     bincounts=[]
     for datalist in datalists:
          bincounts+=[[datalist.count(binv) for binv in bins]]
     return bins,bincounts

def compressandmergebins(bins,bincounts,numbinstomerge=2):
     positionstomerge=[(pos,pos+numbinstomerge-1) for pos in xrange(0,len(bins)-(numbinstomerge-1),numbinstomerge)]
     if positionstomerge[-1][-1]==len(bins):extrabin=(positionstomerge[-1][-1],len(bins)-1)
     elif (positionstomerge[-1][-1]+1)<len(bins):extrabin=(positionstomerge[-1][-1]+1,len(bins)-1)# in case num bins isnt evenyl dividably by numbinstomerge
     else:extrabin=None
     newbins=[(bins[pos[0]],bins[pos[1]]) for pos in positionstomerge]
     if type(extrabin)!=type(None):# if not evenyl divided two cases to consider
          if extrabin[0]==extrabin[1]:newbins.append((bins[extrabin[0]],))
          else:newbins.append((bins[extrabin[0]],bins[-1]))
     newbincounts=[]
     for bincount in bincounts:
          newcount=[sum(bincount[pos[0]:pos[1]+1]) for pos in positionstomerge]
          if type(extrabin)!=type(None):# if not evenyl divided two cases to consider
               if extrabin[0]==extrabin[1]:newcount.append(bincount[extrabin[0]])
               else:newcount.append(sum(bincount[extrabin[0]:extrabin[1]+1]))
          newbincounts.append(newcount)
     return newbins,newbincounts



#print(
"""
To compare sets of data from three imputation methods use following set of commands, Note Fisher_Exact can use more memory than available if too many bins:
basebins=turndatainbins([one64qdata,linreg64qdata,fit64qdata])
compress2=compressandmergebins(basebins[0],basebins[1])
FisherExact.Fisher.fisher_exact(zip(*compress2[1][:2]))#compare 1st and 2nd one and lin
FisherExact.Fisher.fisher_exact(zip(*compress2[1][1:]))#compare 3rd and 2nd fit and lin
FisherExact.Fisher.fisher_exact(zip(*[compress2[1][0]]+[compress2[1][2]]))#compare 3rd and 1st fit and one
"""#)
def turntreestatdicsintoplots_strippedaxis_meanvaluesscatterplotxlogoverlap(statisticdic,selectedmetric,dataandoriginalimptypes,colors=["b","g","r"],colorder=None,hardyrange=None,coltitle=None,*args):
        import numpy as np
        metricnames=[selectedmetric]
        rawentries=[entry for entry in statisticdic[selectedmetric]]
        listofinputargforimputation=[set(inputargs) for inputargs in zip(*rawentries)]
        datasetnames=list(listofinputargforimputation[0])
        blankoutstrs=[("Blanking Method: "+(entry.split("_")[2])) for entry in list(listofinputargforimputation[1])]#center text how?
        blankoutnames=list(listofinputargforimputation[1])
        imputationnames=[("Imp for Data(%s)"%(entry.split("_"))[0]) for entry in list(listofinputargforimputation[3])]
        impfororitables=[]
        for dataimpmeth in list(listofinputargforimputation[3]):
                for entry in dataandoriginalimptypes:
                        if dataimpmeth == entry[1]:impfororitables.append(entry[0])
        print impfororitables
        imputationnames=[(imputationnames[pos]+"and Original(%s)"%(impfororitables[pos].split("_")[0])) for pos in xrange(len(imputationnames))]
        imputationentrynames=list(listofinputargforimputation[3])
        if type(colorder)==list:
                imputationentrynames=colorder
                impfororitables=[]
                for dataimpmeth in imputationentrynames:
                        for entry in dataandoriginalimptypes:
                                if dataimpmeth == entry[1]:impfororitables.append(entry[0])
                print impfororitables
                imputationnames=[("Imp for Data(%s) and Original(%s)"%(imputationentrynames[entry].split("_")[0],impfororitables[entry].split("_")[0])) for entry in xrange(len(imputationentrynames))]
        if type(coltitle)==list:
                imputationnames=coltitle
        numberofblankoutnames=list(listofinputargforimputation[2])
        numberofblankedoutentries=sorted([(type(2))(strpart) for entry in numberofblankoutnames for strpart in entry.split() if strpart.isdigit()])
        temp=[]
        for item in numberofblankedoutentries:
                for posstring in numberofblankoutnames:
                        if (" "+str(item)+" ") in posstring:temp.append(posstring)
        numberofblankoutnames=temp
        print(numberofblankoutnames)
        print(numberofblankedoutentries)

        rows=blankoutstrs
        cols=imputationnames

        pad = 10 # in points

        allfigs=[]


        #################################
        #selectedmetric=metricnames[0]

        #################################


        for selectedmetric in metricnames:

                maxy,miny=0,10000000
                maxx,minx=0,10000000

                fig, axes = plt.subplots(nrows=len(rows), ncols=1)#(,figsize=(len(cols)*2,len(rows)*4))# how to make bar and whisker plot in subplots

                #for ax, col in zip(axes[0], cols):
                #    ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                #                xycoords='axes fraction', textcoords='offset points',
                #                size='large', ha='center', va='baseline')#,size=fontsize)


                plt.setp(axes.flat, xlabel='Number of Missing Entries', ylabel='Metric: %s'%(selectedmetric))

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                print("working on col %d and row %d"%(colpos,rowpos))
                                sorteddata=[]
                                currentimputationmethod=imputationentrynames[colpos]
                                currentcurrentblankpickmethod=blankoutnames[rowpos]
                                for nmepos in xrange(len(numberofblankedoutentries)):
                                        print("begin for this nme %d"%(numberofblankedoutentries[nmepos]))
                                        datafornmeimpandblankout=[]
                                        for accumalteddataname in statisticdic[selectedmetric]:
                                                if (imputationentrynames[colpos] in accumalteddataname) and (blankoutnames[rowpos] in accumalteddataname) and (numberofblankoutnames[nmepos] in accumalteddataname)  and (impfororitables[colpos] in accumalteddataname):
                                                        for metricvalue in statisticdic[selectedmetric][accumalteddataname]: datafornmeimpandblankout.append(metricvalue)
                                        sorteddata.append(datafornmeimpandblankout)
                                        print("data for %s and %s and %d: \n"%(currentimputationmethod,currentcurrentblankpickmethod,numberofblankedoutentries[nmepos]))
                                        print(datafornmeimpandblankout)
                                boxplotdic=axes[rowpos].semilogx(numberofblankedoutentries, [sum(data)/len(data) for data in sorteddata],color=colors[colpos],label=currentimputationmethod)
                                axes[rowpos].errorbar(numberofblankedoutentries,[sum(data)/len(data) for data in sorteddata] , yerr=[np.std(data) for data in sorteddata],color=colors[colpos])                                
                                rangey=axes[rowpos].get_ylim()
                                rangex=axes[rowpos].get_xlim()
                                if rangey[1]>maxy:maxy=rangey[1]
                                if rangey[0]<miny:miny=rangey[0]
                                if rangex[1]>maxx:maxx=rangex[1]
                                if rangex[0]<minx:minx=rangex[0]

                for colpos in xrange(len(imputationentrynames)):
                        for rowpos in xrange(len(blankoutnames)):
                                if type(hardyrange)==list:axes[rowpos].set_ylim(hardyrange[0],hardyrange[1])
                                else:axes[rowpos].set_ylim(miny,maxy)
                                axes[rowpos].set_xlim(minx,maxx)
                                handles,labels=axes[rowpos].get_legend_handles_labels()
                                axes[rowpos].legend(handles,labels)
                                

                
                                
                                
                fig.tight_layout(pad=1,h_pad=.5,w_pad=.5)#  
                allfigs.append((fig, axes))

                plt.show()

                print (numberofblankoutnames)
                print (sorteddata)
        return allfigs

def turntreestatdicsintoplots_strippedaxis_meanvaluesscatterplotxlogoverlapformanydatasets(statisticdics,selectedmetric,dataandoriginalimptypes,datasetnames=["a","b","c"],colors=["b","g","r"],colorder=None,hardyrange=None,coltitle=None,*args):
        import numpy as np
        metricnames=[selectedmetric]
        initstatdic=statisticdics[0]
        statisticdic=initstatdic
        rawentries=[entry for entry in statisticdic[selectedmetric]]
        listofinputargforimputation=[set(inputargs) for inputargs in zip(*rawentries)]
        #datasetnames=list(listofinputargforimputation[0])
        blankoutstrs=[("Blanking Method: "+(entry.split("_")[2])) for entry in list(listofinputargforimputation[1])]#center text how?
        blankoutnames=list(listofinputargforimputation[1])
        imputationnames=[("%s"%(entry.split("_"))[0]) for entry in list(listofinputargforimputation[3])]
        impfororitables=[]
        for dataimpmeth in list(listofinputargforimputation[3]):
                for entry in dataandoriginalimptypes:
                        if dataimpmeth == entry[1]:impfororitables.append(entry[0])
        #print impfororitables
        imputationnames=[(imputationnames[pos]+"and Original(%s)"%(impfororitables[pos].split("_")[0])) for pos in xrange(len(imputationnames))]
        imputationentrynames=list(listofinputargforimputation[3])
        if type(colorder)==list:
                imputationentrynames=colorder
                impfororitables=[]
                for dataimpmeth in imputationentrynames:
                        for entry in dataandoriginalimptypes:
                                if dataimpmeth == entry[1]:impfororitables.append(entry[0])
                #print impfororitables
                imputationnames=[("%s"%(imputationentrynames[entry].split("_")[0])) for entry in xrange(len(imputationentrynames))]
        if type(coltitle)==list:
                imputationnames=coltitle
        numberofblankoutnames=list(listofinputargforimputation[2])
        numberofblankedoutentries=sorted([(type(2))(strpart) for entry in numberofblankoutnames for strpart in entry.split() if strpart.isdigit()])
        temp=[]
        for item in numberofblankedoutentries:
                for posstring in numberofblankoutnames:
                        if (" "+str(item)+" ") in posstring:temp.append(posstring)
        numberofblankoutnames=temp
        print(numberofblankoutnames)
        print(numberofblankedoutentries)

        rows=datasetnames
        cols=blankoutstrs
        #print datasetnames,blankoutstrs,rows,cols

        pad = 10 # in points

        allfigs=[]


        maxy,miny=0,10000000
        maxx,minx=0,10000000

        fig, axes = plt.subplots(nrows=len(statisticdics), ncols=len(blankoutstrs))#(,figsize=(len(cols)*2,len(rows)*4))# how to make bar and whisker plot in subplots
        #print axes, len(statisticdics), len(blankoutstrs),len(rows),len(cols)

        for ax, col in zip(axes[0], cols):
                ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                                xycoords='axes fraction', textcoords='offset points',
                                size='large', ha='center', va='baseline')#,size=fontsize)


        for ax, row in zip(axes[:,0], rows):
                ax.annotate(row, xy=(-.22, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                        xycoords=ax.yaxis.label, textcoords='offset points',
                        size='large', ha='center', va='center',rotation=90)#,size=fontsize)

        plt.setp(axes.flat, xlabel='Number of Missing Entries')#, ylabel='Metric: %s'%(selectedmetric))

        #for statisticdic in statisticdics:#rows

        #for selectedmetric in metricnames:#cols

        assert(len(statisticdics)==len(rows))
        assert(len(statisticdics)==len(datasetnames))
        assert(len(blankoutnames)==len(blankoutstrs))

        for colpos in xrange(len(blankoutstrs)):
                for rowpos in xrange(len(datasetnames)):#to do imputationentrynames
                        print("working on col %d and row %d"%(colpos,rowpos))
                        #currentimputationmethod=imputationentrynames[colpos]
                        currentcurrentblankpickmethod=blankoutnames[colpos]
                        statisticdic=statisticdics[rowpos]
                        selectedmetric=metricnames[0]#only one metric
                        #switch nmepos and impmeth to do each method one at a time?
                        for impmethodpos in xrange(len(imputationentrynames)):
                                sorteddata=[]
                                
                                for nmepos in xrange(len(numberofblankedoutentries)):
                                        print("begin for this nme %d and imputation method %s"%(numberofblankedoutentries[nmepos],imputationentrynames[impmethodpos]))
                                        currentimputationmethod=imputationentrynames[impmethodpos]
                                        datafornmeimpandblankout=[]
                                        for accumalteddataname in statisticdic[selectedmetric]:
                                                #print "accumalteddataname", accumalteddataname, len(accumalteddataname)
                                                #print imputationentrynames[impmethodpos],blankoutnames[colpos],numberofblankoutnames[nmepos],impfororitables[impmethodpos]
                                                if len(accumalteddataname)==4:
                                                        if "additional" in numberofblankoutnames[nmepos]:numberofblankoutnames[nmepos]=numberofblankoutnames[nmepos].replace( " additional missingentries"," missingentries")
                                                        if (imputationentrynames[impmethodpos] in accumalteddataname) and (blankoutnames[colpos] in accumalteddataname) and (numberofblankoutnames[nmepos] in accumalteddataname):
                                                                #print accumalteddataname, "hit"
                                                                for metricvalue in statisticdic[selectedmetric][accumalteddataname]: datafornmeimpandblankout.append(metricvalue)
                                                elif len(accumalteddataname)==5:
                                                        if "additional" not in numberofblankoutnames[nmepos]:numberofblankoutnames[nmepos]=numberofblankoutnames[nmepos].replace(" missingentries", " additional missingentries")
                                                        if (imputationentrynames[impmethodpos] in accumalteddataname) and (blankoutnames[colpos] in accumalteddataname) and (numberofblankoutnames[nmepos] in accumalteddataname)  and (impfororitables[impmethodpos] in accumalteddataname):
                                                                #print accumalteddataname, "hit"
                                                                for metricvalue in statisticdic[selectedmetric][accumalteddataname]: datafornmeimpandblankout.append(metricvalue)
                                                
                                                else:
                                                        assert(False)

                                        sorteddata.append(datafornmeimpandblankout)
                                        print("data for %s and %s and %d: \n"%(currentimputationmethod,currentcurrentblankpickmethod,numberofblankedoutentries[nmepos]))
                                        #print(datafornmeimpandblankout)
                                #print numberofblankedoutentries, [sum(data)/len(data) for data in sorteddata],colors[impmethodpos],currentimputationmethod, "pick color"
                                #print len(numberofblankedoutentries),len([sum(data)/len(data) for data in sorteddata])
                                #for data in sorteddata:print len(data) 
                                boxplotdic=axes[rowpos][colpos].semilogx(numberofblankedoutentries, [sum(data)/len(data) for data in sorteddata],color=colors[impmethodpos],label=currentimputationmethod)
                                axes[rowpos][colpos].errorbar(numberofblankedoutentries,[sum(data)/len(data) for data in sorteddata] , yerr=[np.std(data) for data in sorteddata],color=colors[impmethodpos])                                
                                rangey=axes[rowpos][colpos].get_ylim()
                                rangex=axes[rowpos][colpos].get_xlim()
                                if rangey[1]>maxy:maxy=rangey[1]
                                if rangey[0]<miny:miny=rangey[0]
                                if rangex[1]>maxx:maxx=rangex[1]
                                if rangex[0]<minx:minx=rangex[0]

        for colpos in xrange(len(blankoutstrs)):
                for rowpos in xrange(len(datasetnames)):#to do imputationentrynames
                        if type(hardyrange)==list:axes[rowpos][colpos].set_ylim(hardyrange[0],hardyrange[1])
                        else:axes[rowpos][colpos].set_ylim(miny,maxy)
                        axes[rowpos][colpos].set_xlim(minx,maxx)
                        handles,labels=axes[rowpos][colpos].get_legend_handles_labels()
                        axes[rowpos][colpos].legend(handles,labels,loc='upper left')
                        

        
                        
                        
        fig.tight_layout(pad=1,h_pad=.5,w_pad=.5)#  
        allfigs.append((fig, axes))

        plt.show()

        print (numberofblankoutnames)
        print (sorteddata)
        #print colors
        return allfigs
#dataandoriginalimptypes=[('LinearReg_imputedoriginal',"LinearReg_imputedtables"),('OneNN_imputedoriginal',"OneNN_imputedtables"),('FitchImp_imputedoriginal',"FitchImp_imputedtables")]
#os.chdir("Frumkin2005_treea_processedandreformateddata_Fri Feb 10 2017");os.chdir("tree_datasets");fil=open("savedtreedatastructs.pkl");a=pickle.load(fil);treeadata=pickle.load(fil);fil.close();os.chdir("..");os.chdir("..");
#os.chdir("Frumkin2005_treeb_processedandreformateddata_Fri Feb 10 2017");os.chdir("tree_datasets");fil=open("savedtreedatastructs.pkl");a=pickle.load(fil);treebdata=pickle.load(fil);fil.close();os.chdir("..");os.chdir("..");
#os.chdir("Frumkin2005_treec_processedandreformateddata_Fri Feb 10 2017");os.chdir("tree_datasets");fil=open("savedtreedatastructs.pkl");a=pickle.load(fil);treecdata=pickle.load(fil);fil.close();os.chdir("..");os.chdir("..");
#see=turntreestatdicsintoplots_strippedaxis_meanvaluesscatterplotxlogoverlapformanydatasets([treeadata,treebdata,treecdata],"R-F",dataandoriginalimptypes,datasetnames=["Tree A","Tree B","Tree C"],colors=["b","g","r"],colorder=["LinearReg_imputedtables","OneNN_imputedtables","FitchImp_imputedtables"],hardyrange=None,coltitle=None);
#seea=turntreestatdicsintoplots_strippedaxis_meanvaluesscatterplotxlogoverlapformanydatasets([treeadata,treeadata],"R-F",dataandoriginalimptypes,datasetnames=["Tree A","A"],colors=["b","g","r"],colorder=["LinearReg_imputedtables","OneNN_imputedtables","FitchImp_imputedtables"],hardyrange=None,coltitle=None);
