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

#CONDA_ENV_NAME='analysis_env'
check_conda_and_run(){
    script_to_run=$0
    echo -e "Checking if the conda environment $1 is installed"
    conda_env=$(conda env list | grep "envs" | awk -F "*" '{print $2}' | awk -F "/" '{print $5}')
    #echo $conda_env
    if [[ $conda_env == $1 ]] ;
    then
        echo -e "  - Yay, we're in the correct conda env $1!"
        echo -e "${GREEN}Running ${script_to_run}${NC}"
        #python ${script_to_run}
    else
        #echo -e "Activating env $1"
        source activate $1 
        #echo -e "   - you need to conda activate $1 because subshells"
    fi
}

check_conda_and_run2(){
    script_to_run=$0
    #echo -e "Checking if the conda environment is $1"
    conda_env=$(conda env list | grep "envs" | awk -F "*" '{print $2}' | awk -F "/" '{print $5}')
    #echo $conda_env
    if [[ $conda_env == $1 ]] ;
    then
        echo -e "  - Dengue pipeline is installed now is starting analysis!"
        echo -e "${GREEN}Running ${script_to_run}${NC}"
        #python ${script_to_run}
    else
        
        echo -e " ################ Is required to install $1 this will take more than 10 mins ################# "
        cd /home/jovyan/databases/DENV_pipeline 
        conda env create -f environment.yml
        source activate analysis_env
        pip install /home/jovyan/databases/DENV_pipeline/.
        #echo -e "   - you need to conda activate $1 because subshells"
    fi
}

check_conda_and_run analysis_env


check_conda_and_run2 analysis_env
#source activate analysis_env

dir=$1
echo -e "working directory was set to $dir"
if [ ! -f $dir/*R1* ]; then
echo "###################################################################"
echo "#####Failed to  find illumina paired read in folder $dir###########"
echo "###################################################################"
exit 1
fi

/home/jovyan/shared/data/copingfilesdenv.sh $dir


denv_pipeline --indir $dir/denv --cores 10 --outdir $dir/denv 