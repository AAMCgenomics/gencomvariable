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


import fileinput
import sys
import os
#print("library loaded")

tsvfile=sys.argv[1]
file1=sys.argv[2]
#print(file1)
#print(tsvfile)
#outputcsv=os.path.dirname(file1)+"/listfile.txt"

d={}
with open(tsvfile) as f:
    for line in f:
        (value, key)=line.split()
        #print(key)
        d[key]=value
        #key = int(val)
        #I have to eliminate the cities which have less than 10,000 people
        #d.append(key,value)
    #return {}
    #You want to return the created dictionary

folder=file1
if folder in d:
    print(d[folder])
else:
    print('barcodenotfound')