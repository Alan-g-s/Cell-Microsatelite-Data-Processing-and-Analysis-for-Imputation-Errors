#!/usr/bin/python2.7

import openpyxl as px
from string import ascii_uppercase

allletters=list(ascii_uppercase)

def countforlevel(level):
    """function that counts how many total letters combinations can be made from at most "level" chars at once"""
    assert((type(level).__name__)=='int')
    assert(level>-1)
    numofchars=26
    count=0
    for eachlevel in xrange(level+1):
        count+=numofchars**(eachlevel+1)
    return count

def countnumlevelsrequired(count):
    """function that determines how combinations of letters (in combination of at most "levelcounter" letters at once ) are needed to uniquely label all the cells you have"""
    levelcounter=0
    levelcount=countforlevel(levelcounter)
    while(count>levelcount):
        levelcounter+=1
        levelcount=countforlevel(levelcounter)
    return levelcounter

def makecollabelsforexcel(level):
    """function that returns a list of all column names that are used in a excel file in a recursive manner"""
    assert((type(level).__name__)=='int')
    assert(level>-1)
    if level==0:
        return list(allletters)
    else:
        newlabels=list(allletters)
        #counter=0
        priorlist=makecollabelsforexcel(level-1)
        for string in priorlist:
            for char in list(allletters):
                #counter+=1
                newlabels.append(''.join([string,char]))
        return newlabels


def turntable2excelfile(datatable,filename,sheetname):
    """function that takes in a data table as a 2d array like object and save that table in an excel file named "filename" in sheet "sheetname" """
    WW=px.Workbook()
    pp=WW.get_active_sheet()
    pp.title=sheetname

    # this functions make the column 'names' that excel uses
    levelforcols=countnumlevelsrequired(len(datatable[0]))
    excellabels=makecollabelsforexcel(levelforcols)

    # this actually saves the data to cells in an excel file
    for row in xrange(len(datatable)):
        for col in xrange(len(datatable[row])):
            pp.cell('%s%d'%(excellabels[col],row+1)).value=datatable[row][col]# necessary either way

    WW.save(filename)
