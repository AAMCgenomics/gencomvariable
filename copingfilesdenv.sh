#!/bin/bash
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

#This software uses several smart programs prepared and collected elsewhere each one retain its particular license, please comply with them.

dir=$1

if [ -z $dir   ]; then

echo " ############## directory argument not provided ############"

exit 1

fi


if [ ! -f $dir/*R1* ]; then
echo 
exit 1
fi


if [ ! -d $dir/denv ]; then

	mkdir $dir/denv

fi

(
for read_1 in $dir/*R1*;
do
    start=`date +%s`
    read_2=${read_1/R1/R2}
   longname=$(basename $read_1)
   c=$(echo $longname | awk -F "_S" '{print $1}')
   echo -e " Creating folder and copying sample $c"
   if [ ! -d $dir/denv/$c ]; then
	mkdir $dir/denv/$c
    fi
    cp $read_1 $dir/denv/$c/.
    cp $read_2 $dir/denv/$c/.
   
done
wait
)