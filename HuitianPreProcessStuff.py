#!/usr/bin/env python3
# Last update: Jan 9, 2022, By Huitian Diao


import numpy as np
import pandas as pd
import scanpy as sc
import os
import csv

### Self defined functions for identifying cell antibody hashtags
def findAb(ftTsv):
    AbDict = {}
    with open(ftTsv, "r") as fin:
        rfin = csv.reader(fin, delimiter="\t")
        rowN = 0
        for row in rfin:
            rowN += 1
            if row[2] == "Antibody Capture":
                AbDict[rowN] = row[0]
    return(AbDict)

def cellAbReads(root_name, ftTsv, mtxMtx):
    outName = "%s_Cells_hashTags_count.csv"%root_name
    abDict = findAb(ftTsv)
    with open(mtxMtx, "r") as fin:
        with open(outName, "w") as fout:
            rfin = csv.reader(fin, delimiter=" ")
            wfout = csv.writer(fout, delimiter=",")
            wfout.writerow(["Cell"] + list(abDict.values()))
            
            next(rfin)
            next(rfin)
            next(rfin)
            step = 1
            outRow = [1] + [0 for x in abDict]
            for row in rfin:
                cellN = int(row[1])
                if cellN != step:
                    wfout.writerow(outRow)
                    step = cellN
                    outRow = [row[1]] + [0 for x in abDict]
                if int(row[0]) in list(abDict.keys()):
                    abIdx = list(abDict.keys()).index(int(row[0])) + 1
                    outRow[abIdx] = int(row[2])
            wfout.writerow(outRow)

def cellType(root_name, count_cutoff, amb_cutoff):
    cellTag = "%s_Cells_hashTags_count.csv"%root_name
    cellTag_norm = "%s_Cells_hashTags_count_norm.csv"%root_name
    outfile = "%s_Cells_hashTags.csv"%root_name
    
    count_df = pd.read_csv(cellTag, index_col=0)
    norm_count_df = count_df / count_df.sum() * 10000
    norm_count_df = norm_count_df.div(norm_count_df.sum(axis=1), axis=0) * 100
    norm_count_df = norm_count_df.fillna(0)
    norm_count_df.to_csv(cellTag_norm)
    
    with open(cellTag, "r") as fin:
        with open(cellTag_norm, "r") as fin_norm:
            with open(outfile, "w") as fout:
                rfin = csv.reader(fin, delimiter=",")
                rfin_norm = csv.reader(fin_norm, delimiter=",")
                all_types = next(rfin)[1:]
                next(rfin_norm)
                
                wfout = csv.writer(fout, delimiter=",")
                cell_types = []
                for row in rfin:
                    row_nu = [float(x) for x in row[1:]]
                    row_norm = next(rfin_norm)
                    row_norm_nu = [float(x) for x in row_norm[1:]]
                    rowMax = max(row_nu)
                    row_normMax = max(row_norm_nu)
                    row_type = all_types[row_norm_nu.index(row_normMax)]
                    if row_normMax/100 < amb_cutoff:
                        #print('Error in cell %s, maximum hashtag reads less than 70 percent...' %row[0])
                        #print(row)
                        row_type = "Doublet"
                    if rowMax < count_cutoff:
                        row_type = "Negative"


                    wfout.writerow([row_type])
                    cell_types.append(row_type)
    types_set = list(set(cell_types))
    
    
    sum_df = pd.DataFrame({"sample":[root_name], "total_cell_number":[len(cell_types)],
               "negative_percentage":[round(float(cell_types.count("Negative"))*100/len(cell_types), 2)],
               "doublet_percentage":[round(float(cell_types.count("Doublet"))*100/len(cell_types), 2)]})
    for i in types_set:
        if i != "Negative" and i != "Doublet":
            sum_df[i] = [round(float(cell_types.count(i))*100/len(cell_types), 2)]
    
    sum_df.to_csv("%s_hashtagSummary.csv"%root_name)
    print(sum_df)
    
    return(sum_df)
