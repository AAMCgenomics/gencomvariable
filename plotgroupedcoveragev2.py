#!/usr/bin/env python

#Copyright Alexander A. Martinez C & Gorgas Memorial Institute for Health Studies
#Written by: Alexander A. Martinez C, Genomics and proteomics research unit, Gorgas memorial #Institute For Health Studies.
#Licensed under the Apache License, Version 2.0 (the "License"); you may not use
#this work except in compliance with the License. You may obtain a copy of the
#License at:
#http://www.apache.org/licenses/LICENSE-2.0
#Unless required by applicable law or agreed to in writing, software distributed
#under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
#CONDITIONS OF ANY KIND, either express or implied. See the License for the
#specific language governing permissions and limitations under the License.

#This software uses diferent smart programs prepared and collected elsewhere each one retain its particular licence, please comply with them.



#This script takes all the fastq files within a folder and according to the pathogen chose #clean the reads, pileup them and generate several outputs. Currently the pathogen with more #application and outputs availables is SARS-CoV-2


import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import fileinput
import sys
import glob

file1=sys.argv[2]
directory=sys.argv[1]
name=sys.argv[3]
#outputcsv=sys.argv[4]
inputcsv=file1
outputcsv=directory+"/"+name+"_depth_plot.png"
outputcsv2=directory+"/"+name+"_depth_plot2.png"
#print(inputcsv)

text = open(inputcsv, "r")

sns.set_theme(style="ticks")


df = pd.read_csv(inputcsv, sep=' ')

df.columns =['Sample', 'pos', 'depth']
#print(df)
# Initialize a grid of plots with an Axes for each walk
grid = sns.FacetGrid(df, col="Sample", hue="Sample", palette="tab20c",
                    col_wrap=4, height=2)
#print(grid)
# Draw a horizontal line to show the starting point
grid.refline(y=10, linestyle=":")

grid.refline(y=100, linestyle=":")


# Draw a line plot to show the trajectory of each random walk
grid.map(plt.plot, "pos", "depth")

grid.refline(y=1000, linestyle=":")


print(df['pos'].max())
withsize=[]
if(df['pos'].max()>30000):
	withsize = df['pos'].max()/1000
elif(df['pos'].max()>100000):
	withsize = df['pos'].max()/5000
else:
	withsize = df['pos'].max()/50          


# Adjust the tick positions and labels
grid.set(xticks=np.arange(0, df['pos'].max(), step=withsize*10), yticks=[0, 10, 100, 10000],
         xlim=(1, df['pos'].max()), ylim=(3.5,df['depth'].max()), yscale="log")



print(withsize)
        

# Adjust the arrangement of the plots
grid.fig.tight_layout(w_pad=1)
grid.fig.set_figwidth(withsize)
grid.savefig(outputcsv)

grid.set(xticks=np.arange(0, df['pos'].max(), step=df['pos'].max()/10), yticks=[0, 10, 100, 10000],
         xlim=(1, df['pos'].max()), ylim=(3.5,df['depth'].max()), yscale="log")
grid.fig.set_figwidth(withsize/2)
grid.savefig(outputcsv2)