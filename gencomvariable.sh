#/bin/bash

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


usage() { echo "$0 not option added to this script please read the following help or use gencomvariable.sh -h :" && grep " .)\ #" $0; exit 0; }
[ $# -eq 0 ] && usage
while getopts ":hd:a:v:m:t:c:" arg; do
  case $arg in
    a) # use -a to specify Accession reference number to use for mapping.
      #echo "accession is ${OPTARG}"
      accession=${OPTARG}
      ;;
    m) # use -m FALSE to avoid mapping to reference and and instead use a virus tailored analysis for now PTV more to come!!!.
      #echo "accession is ${OPTARG}"
      mapping=${OPTARG}
      if [ -z $mapping  ] ;
      then 
      mapping="TRUE"
      else
      mapping=${OPTARG}
      fi   
      ;;
    d) # Specify directory where R1 and R2 reads are located.
      dir=${OPTARG}
      ;;
    v) # Virus to apply especial analysis currently HIV, PTV, IMRA influenza Illumina and Minion, MPXV and a custom options are available.
      virus=${OPTARG}
      if [ virus="custom" ] ;
      then
      mapping="FALSE"
      else
      1=1
      fi
      ;;
    c) # apply a provided reference in fasta format to reads; this option must be used with options -v custom and -m FALSE.
      reference=${OPTARG}      
      ;;     
    t) # provide a tsv file for ONT Influenza sequencing first column must be sample and second barcode.
      meta=${OPTARG}      
      ;;           
    h | *) # Display help para este pipeline.
      usage
      exit 0
      ;;
  esac
done
#fleching fasta from ncbi
if [[ $mapping = "FALSE" ]] ;
then 
echo "option ### -m FALSE ### selecting skipping mapping to reference starting extra analysis if virus was indicated."
else

echo "downloading $accession as reference"
mkdir -p $dir/analysis
mkdir -p $dir/analysis/coverage
#### download esearch from genebank read this link for details https://www.ncbi.nlm.nih.gov/books/NBK179288/ ####

/home/jovyan/shared/edirect/esearch -db nuccore -query $accession | /home/jovyan/shared/edirect/efetch -format fasta > $dir/$accession.ref.fasta

REF_NAME=`cat $dir/$accession.ref.fasta | grep '>' | tr -d '>' | cut -d ' ' -f 1`

LENGTH=`tail -n +2 $dir/$accession.ref.fasta | tr -d '\n' | wc -m | xargs`
echo -e "$REF_NAME\t$LENGTH" > $dir/analysis/my.genome

#### this loop will look for fastq R1 and R2 reads localed in a $dir, currently all steps are intended to be use with paired-end reads ####


for read_1 in $dir/*R1*
do #echo $read_1
read_2=${read_1/R1/R2}
echo $read_1
echo $read_2
longname=$(basename $read_1)
echo $longname

c=$(echo $longname | awk -F "_S" '{print $1}')
    echo $c
#d=$(echo $c | awk -F "/" '{print $7}')
    #echo $d
d=$c
echo "starting $d mapping with $accession Analysis"
echo "############################################"
now="$(date +%c)"
echo -e "starting $d $accession Analysis at \t$now" >> "$dir/analysis/mensajes.log"
mkdir -p $dir/analysis/$d
/home/jovyan/shared/data/fastp/fastp -i $read_1 -I $read_2 -o $dir/analysis/$d/$d.R1.fastq.gz -O $dir/analysis/$d/$d.R2.fastq.gz -w 10 -R "fastp report for sample: $d" -h $dir/analysis/$d/$d.fastqreport.html -j $dir/analysis/$d/$d.jsonreport.json 
echo "###############finished fastq step #############################"

#### minimap2 must be installed in a dir and modified in this area. ####

/home/jovyan/shared/data/minimap2-2.24_x64-linux/minimap2 -ax sr $dir/$accession.ref.fasta $dir/analysis/$d/$d.R1.fastq.gz -2 $dir/analysis/$d/$d.R2.fastq.gz -K 2G -t 10 > $dir/analysis/$d/assembled.sam 


samtools view  -bS $dir/analysis/$d/assembled.sam  > $dir/analysis/$d/assembled1.bam

echo '######done samtools view ######'
samtools sort -@ 10 -m 10M $dir/analysis/$d/assembled1.bam -o $dir/analysis/$d/assembled.bam

##$input_dir/assembled.sam $input_dir/assembled1.bam

bamsorted=$dir/analysis/$d/assembled.bam

bedtools bamtobed -i $dir/analysis/$d/assembled.bam > $dir/analysis/$d/reads.bed
bedtools genomecov -bga -i $dir/analysis/$d/reads.bed -g $dir/analysis/my.genome | awk  '$4 < 3' > $dir/analysis/$d/zero.bed
echo '#####bamtobedfinished###'
echo '#####maskingbed###'
maskFastaFromBed -fi $dir/$accession.ref.fasta -bed $dir/analysis/$d/zero.bed -fo $dir/analysis/$d/masked.fasta

echo '####masked reference ready###'
bcftools mpileup -Ou  -f $dir/analysis/$d/masked.fasta  $dir/analysis/$d/assembled.bam | bcftools call --ploidy 1 -mv -Oz -o $dir/analysis/$d/test.vcf.gz
echo "###Mpileup ready###"
bcftools index $dir/analysis/$d/test.vcf.gz
cat $dir/analysis/$d/masked.fasta | bcftools consensus $dir/analysis/$d/test.vcf.gz > $dir/analysis/$d/new_consensus.fasta

sequencename=$(echo $d'_'$accession)


    #echo ">$d" > $c/$d.fas
echo ">$sequencename" > $dir/analysis/$d/$d.fas
echo ">$d" > $dir/analysis/$d/$d.code.fas

tail -n +2 $dir/analysis/$d/new_consensus.fasta >> $dir/analysis/$d/$d.fas
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $dir/analysis/$d/$d.fas | tr "\t" "\n" | sed -r '/^>/! s/N+$|^N+//g' |  fold -w 80 >  $dir/analysis/$d/$d.Ntrimmed.fas

tail -n +2 $dir/analysis/$d/new_consensus.fasta >> $dir/analysis/$d/$d.code.fas
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $dir/analysis/$d/$d.code.fas | tr "\t" "\n" | sed -r '/^>/! s/N+$|^N+//g' |  fold -w 80 >  $dir/analysis/$d/$d.Ncodetrimmed.fas

echo "############################################"
echo "finishing $d consensus generation with bcftools mpileup"
echo "############################################"

samtools depth -aa $bamsorted | awk -v sample=$d '{$1=sample ; print;}' >  $dir/analysis/$d/$d.tsv
samtools coverage $bamsorted >  $dir/analysis/$d/$d.coverage
samtools stats $bamsorted | grep ^SN | cut -f 2- | grep 'average length:' >  $dir/analysis/$d/$d.readlength
awk -v sample=$d '{$1=sample ; print;}' $dir/analysis/$d/$d.tsv > $dir/analysis/coverage/coverage_file_$d.tsv
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/$d/ $dir/analysis/$d/$d.tsv $d
    
done
cat $dir/analysis/coverage/coverage_file_*.tsv > $dir/analysis/coverage/depth_file.tsv 
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/coverage $dir/analysis/coverage/depth_file.tsv $accession
fi 

if [ -z $virus ] ;
then
echo "Not additional analysis specified, if your objective is to analyze HIV illumina reads for genotyping please add #### -v HIV #### option to the script"


exit 1
fi
#####################################################
if [ $virus = "HIV" ] ;
then
###### for HIV currently we just are using the sierrapy pipe for anotating fasta generated in https://hydra.canada.ca/pages/home #####    
for fasta in $dir/analysis/*/*Ncodetrimmed.fas
do
c=$(echo $fasta | awk -F ".Ncodetrimmed.fas" '{print $1}')
    #echo $c
name=$(echo $c | awk -F "/" '{print $8}')
    #echo $d
#echo $fasta
sierrapy fasta $fasta -q /home/jovyan/shared/data/queryfileV1-1.gpl -o $dir/analysis/$name/$name.json
python /home/jovyan/shared/data/jsonconverter_finalv4.py $dir/analysis/$name/$name.0.json $dir/analysis/$name/$name.xml
done

#####################################################
elif [ $virus = "SARS2" ] ;
then

mkdir -p $dir/analysis
mkdir -p $dir/analysis/coverage
reference=/home/jovyan/shared/data/SARS2.fas

REF_NAME=`cat $reference | grep '>' | tr -d '>' | cut -d ' ' -f 1`

LENGTH=`tail -n +2 $reference | tr -d '\n' | wc -m | xargs`
echo -e "$REF_NAME\t$LENGTH" > $dir/analysis/my.genome
#cat $dir/analysis/my.genome
for read_1 in $dir/*R1*
do #echo $read_1
read_2=${read_1/R1/R2}
echo $read_1
echo $read_2
longname=$(basename $read_1)
echo $longname

c=$(echo $longname | awk -F "_S" '{print $1}')
    echo $c
#d=$(echo $c | awk -F "/" '{print $7}')
    echo $d
d=$c
echo "starting $d mapping with $REF_NAME fasta reference provided"
echo "############################################"
now="$(date +%c)"

echo -e "starting $d $REF_NAME Analysis at \t$now" >> "$dir/analysis/mensajes.log"
mkdir -p $dir/analysis/$d
/home/jovyan/shared/data/fastp/fastp -i $read_1 -I $read_2 -o $dir/analysis/$d/$d.R1.fastq.gz -O $dir/analysis/$d/$d.R2.fastq.gz -w 10 -R "fastp report for sample: $d" -h $dir/analysis/$d/$d.fastqreport.html -j $dir/analysis/$d/$d.jsonreport.json
echo "###############finished fastq step for $d reads 1 and 2 #############################"

/home/jovyan/shared/data/minimap2-2.24_x64-linux/minimap2 -ax sr $reference $dir/analysis/$d/$d.R1.fastq.gz -2 $dir/analysis/$d/$d.R2.fastq.gz > $dir/analysis/$d/assembled.sam 


samtools view  -bS $dir/analysis/$d/assembled.sam  > $dir/analysis/$d/assembled1.bam

samtools sort $dir/analysis/$d/assembled1.bam -o $dir/analysis/$d/assembled.bam

#$input_dir/assembled.sam $input_dir/assembled1.bam

bamsorted=$dir/analysis/$d/assembled.bam

bedtools bamtobed -i $dir/analysis/$d/assembled.bam > $dir/analysis/$d/reads.bed
bedtools genomecov -bga -i $dir/analysis/$d/reads.bed -g $dir/analysis/my.genome | awk  '$4 < 3' > $dir/analysis/$d/zero.bed
echo '#####bamtobedfinished###'
echo '#####maskingbed###'

maskFastaFromBed -fi $reference -bed $dir/analysis/$d/zero.bed -fo $dir/analysis/$d/masked.fasta

echo '####masked reference ready###'
bcftools mpileup -Ou  -f $dir/analysis/$d/masked.fasta  $dir/analysis/$d/assembled.bam | bcftools call --ploidy 1 -mv -Oz -o $dir/analysis/$d/test.vcf.gz
echo "###Mpileup ready###"
bcftools index $dir/analysis/$d/test.vcf.gz
cat $dir/analysis/$d/masked.fasta | bcftools consensus $dir/analysis/$d/test.vcf.gz > $dir/analysis/$d/new_consensus.fasta

sequencename=$(echo $d'_'$REF_NAME)

echo ">$sequencename" > $dir/analysis/$d/$d.fas
echo ">$d" > $dir/analysis/$d/$d.code.fas

tail -n +2 $dir/analysis/$d/new_consensus.fasta >> $dir/analysis/$d/$d.fas
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $dir/analysis/$d/$d.fas | tr "\t" "\n" | sed -r '/^>/! s/N+$|^N+//g' |  fold -w 80 >  $dir/analysis/$d/$d.Ntrimmed.fas

tail -n +2 $dir/analysis/$d/new_consensus.fasta >> $dir/analysis/$d/$d.code.fas
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $dir/analysis/$d/$d.code.fas | tr "\t" "\n" | sed -r '/^>/! s/N+$|^N+//g' |  fold -w 80 >  $dir/analysis/$d/$d.Ncodetrimmed.fas

echo "############################################"
echo "finishing $d consensus generation with bcftools mpileup"
echo "############################################"

samtools depth -aa $bamsorted | awk -v sample=$d '{$1=sample ; print;}' >  $dir/analysis/$d/$d.tsv
samtools coverage $bamsorted >  $dir/analysis/$d/$d.coverage
samtools stats $bamsorted | grep ^SN | cut -f 2- | grep 'average length:' >  $dir/analysis/$d/$d.readlength
awk -v sample=$d '{$1=sample ; print;}' $dir/analysis/$d/$d.tsv > $dir/analysis/coverage/coverage_file_$d.tsv
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/$d/ $dir/analysis/$d/$d.tsv $d
 


done


#####################################################
elif [ $virus = "custom" ] ;
then
      if [ -z $reference  ] ;
      then 
      echo "Using a custom mapping option however a reference was not provided using -c argument please provide it an try again"
      exit 1
      else
      echo "Starting a custom mapping using provided  $reference against reads localted in $dir"
      fi  
mkdir -p $dir/analysis
mkdir -p $dir/analysis/coverage

REF_NAME=`cat $reference | grep '>' | tr -d '>' | cut -d ' ' -f 1`

LENGTH=`tail -n +2 $reference | tr -d '\n' | wc -m | xargs`
echo -e "$REF_NAME\t$LENGTH" > $dir/analysis/my.genome
#cat $dir/analysis/my.genome
for read_1 in $dir/*R1*
do #echo $read_1
read_2=${read_1/R1/R2}
echo $read_1
echo $read_2
longname=$(basename $read_1)
echo $longname

c=$(echo $longname | awk -F "_S" '{print $1}')
    echo $c
#d=$(echo $c | awk -F "/" '{print $7}')
    echo $d
d=$c
echo "starting $d mapping with $REF_NAME fasta reference provided"
echo "############################################"
now="$(date +%c)"

echo -e "starting $d $REF_NAME Analysis at \t$now" >> "$dir/analysis/mensajes.log"
mkdir -p $dir/analysis/$d
/home/jovyan/shared/data/fastp/fastp -i $read_1 -I $read_2 -o $dir/analysis/$d/$d.R1.fastq.gz -O $dir/analysis/$d/$d.R2.fastq.gz -w 10 -R "fastp report for sample: $d" -h $dir/analysis/$d/$d.fastqreport.html -j $dir/analysis/$d/$d.jsonreport.json
echo "###############finished fastq step for $d reads 1 and 2 #############################"

/home/jovyan/shared/data/minimap2-2.24_x64-linux/minimap2 -ax sr $reference $dir/analysis/$d/$d.R1.fastq.gz -2 $dir/analysis/$d/$d.R2.fastq.gz > $dir/analysis/$d/assembled.sam 

echo "Minimap finished"

samtools view  -bS $dir/analysis/$d/assembled.sam  > $dir/analysis/$d/assembled1.bam

echo "view finished"

samtools sort  $dir/analysis/$d/assembled1.bam -o $dir/analysis/$d/assembled.bam


echo "sort finished"
#$input_dir/assembled.sam $input_dir/assembled1.bam

bamsorted=$dir/analysis/$d/assembled.bam

bedtools bamtobed -i $dir/analysis/$d/assembled.bam > $dir/analysis/$d/reads.bed
bedtools genomecov -bga -i $dir/analysis/$d/reads.bed -g $dir/analysis/my.genome | awk  '$4 < 3' > $dir/analysis/$d/zero.bed
echo '#####bamtobedfinished###'
echo '#####maskingbed###'

maskFastaFromBed -fi $reference -bed $dir/analysis/$d/zero.bed -fo $dir/analysis/$d/masked.fasta

echo '####masked reference ready###'
bcftools mpileup -Ou  -f $dir/analysis/$d/masked.fasta  $dir/analysis/$d/assembled.bam | bcftools call --ploidy 1 -mv -Oz -o $dir/analysis/$d/test.vcf.gz
echo "###Mpileup ready###"
bcftools index $dir/analysis/$d/test.vcf.gz
cat $dir/analysis/$d/masked.fasta | bcftools consensus $dir/analysis/$d/test.vcf.gz > $dir/analysis/$d/new_consensus.fasta

sequencename=$(echo $d'_'$REF_NAME)

echo ">$sequencename" > $dir/analysis/$d/$d.fas
echo ">$d" > $dir/analysis/$d/$d.code.fas

tail -n +2 $dir/analysis/$d/new_consensus.fasta >> $dir/analysis/$d/$d.fas
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $dir/analysis/$d/$d.fas | tr "\t" "\n" | sed -r '/^>/! s/N+$|^N+//g' |  fold -w 80 >  $dir/analysis/$d/$d.Ntrimmed.fas

tail -n +2 $dir/analysis/$d/new_consensus.fasta >> $dir/analysis/$d/$d.code.fas
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $dir/analysis/$d/$d.code.fas | tr "\t" "\n" | sed -r '/^>/! s/N+$|^N+//g' |  fold -w 80 >  $dir/analysis/$d/$d.Ncodetrimmed.fas

echo "############################################"
echo "finishing $d consensus generation with bcftools mpileup"
echo "############################################"

samtools depth -aa $bamsorted | awk -v sample=$d '{$1=sample ; print;}' >  $dir/analysis/$d/$d.tsv
samtools coverage $bamsorted >  $dir/analysis/$d/$d.coverage
samtools stats $bamsorted | grep ^SN | cut -f 2- | grep 'average length:' >  $dir/analysis/$d/$d.readlength
awk -v sample=$d '{$1=sample ; print;}' $dir/analysis/$d/$d.tsv > $dir/analysis/coverage/coverage_file_$d.tsv
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/$d/ $dir/analysis/$d/$d.tsv $d
 



done


##############################################################
elif  [ $virus = "PTV" ] ; 
then

mkdir -p $dir/analysis
mkdir -p $dir/analysis/coverage

#/home/jovyan/edirect/esearch -db nuccore -query "MK896481.1" | /home/jovyan/edirect/efetch -format fasta > $dir/PTV_gene_S.ref.fasta
cp /home/jovyan/shared/data/PTgeneSreference.fas $dir/PTV_gene_S.ref.fasta

REF_NAME1=`cat $dir/PTV_gene_S.ref.fasta | grep '>' | tr -d '>' | cut -d ' ' -f 1`

LENGTH1=`tail -n +2 $dir/PTV_gene_S.ref.fasta | tr -d '\n' | wc -m | xargs`
echo -e "$REF_NAME1\t$LENGTH1" > $dir/analysis/PTV_gene_S.genome




#/home/jovyan/edirect/esearch -db nuccore -query "MK896482.1" | /home/jovyan/edirect/efetch -format fasta > $dir/PTV_gene_M.ref.fasta
cp /home/jovyan/shared/data/PTgeneMreference.fas $dir/PTV_gene_M.ref.fasta

REF_NAME2=`cat $dir/PTV_gene_M.ref.fasta | grep '>' | tr -d '>' | cut -d ' ' -f 1`

LENGTH2=`tail -n +2 $dir/PTV_gene_M.ref.fasta | tr -d '\n' | wc -m | xargs`
echo -e "$REF_NAME2\t$LENGTH2" > $dir/analysis/PTV_gene_M.genome

#/home/jovyan/edirect/esearch -db nuccore -query "MK896483.1" | /home/jovyan/edirect/efetch -format fasta > $dir/PTV_gene_L.ref.fasta
cp /home/jovyan/shared/data/PTgeneLreference.fas $dir/PTV_gene_L.ref.fasta

REF_NAME3=`cat $dir/PTV_gene_L.ref.fasta | grep '>' | tr -d '>' | cut -d ' ' -f 1`

LENGTH3=`tail -n +2 $dir/PTV_gene_L.ref.fasta | tr -d '\n' | wc -m | xargs`
echo -e "$REF_NAME3\t$LENGTH3" > $dir/analysis/PTV_gene_L.genome

for read_1 in $dir/*R1*
do #echo $read_1
read_2=${read_1/R1/R2}
echo $read_1
echo $read_2
longname=$(basename $read_1)
echo $longname


c=$(echo $longname | awk -F "_S" '{print $1}')
    #echo $c
d=$c
#d=$(echo $c | awk -F "/" '{print $7}')
    #echo $d
echo "starting $d Punta toro Virus mapping Analysis"
echo "############################################"
now="$(date +%c)"
echo -e "starting $d SARS-Cov-2019 Analysis at \t$now" >> "$dir/analysis/mensajes.log"
mkdir -p $dir/analysis/$d

for ref in $dir/*fasta
do 
e=$(echo $ref | awk -F ".ref.fasta" '{print $1}')
    #echo $c
gene=$(echo $e | awk -F "/" '{print $7}')
    #echo $d
/home/jovyan/shared/data/minimap2-2.24_x64-linux/minimap2 -ax sr $ref $read_1 -2 $read_2 > $dir/analysis/$d/$gene.assembled.sam 


samtools view  -bS $dir/analysis/$d/$gene.assembled.sam  > $dir/analysis/$d/$gene.assembled1.bam

samtools sort $dir/analysis/$d/$gene.assembled1.bam -o $dir/analysis/$d/$gene.assembled.bam

echo "after samtools sort ##### "

bamsorted=$dir/analysis/$d/$gene.assembled.bam

bedtools bamtobed -i $dir/analysis/$d/$gene.assembled.bam > $dir/analysis/$d/$gene.reads.bed

echo "#####after bamtobed #####"


bedtools genomecov -bga -i $dir/analysis/$d/$gene.reads.bed -g $dir/analysis/$gene.genome | awk  '$4 < 3' > $dir/analysis/$d/$gene.zero.bed

echo '#####maskingbed###'

maskFastaFromBed -fi $ref -bed $dir/analysis/$d/$gene.zero.bed -fo $dir/analysis/$d/$gene.masked.fasta

echo "####masked $gene reference ready###"

bcftools mpileup -Ou  -f $dir/analysis/$d/$gene.masked.fasta  $dir/analysis/$d/$gene.assembled.bam | bcftools call --ploidy 1 -mv -Oz -o $dir/analysis/$d/test.vcf.gz
echo "###Mpileup ready###"
bcftools index $dir/analysis/$d/test.vcf.gz
cat $dir/analysis/$d/$gene.masked.fasta | bcftools consensus $dir/analysis/$d/test.vcf.gz > $dir/analysis/$d/$gene.consensus.fasta

sequencename=$(echo $d'_'$gene)
echo ">$sequencename" > $dir/analysis/$d/$sequencename.fas
tail -n +2 $dir/analysis/$d/$gene.consensus.fasta >> $dir/analysis/$d/$sequencename.fas

awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $dir/analysis/$d/$sequencename.fas | tr "\t" "\n" | sed -r '/^>/! s/N+$|^N+//g' |  fold -w 80 >  $dir/analysis/$d/$sequencename.Ntrimmed.fas

echo "############################################"
echo "finishing $d consensus with $gene generation with bcftools mpileup"
echo "############################################"

samtools depth -aa $bamsorted | awk -v sample=$sequencename '{$1=sample ; print;}' >  $dir/analysis/$d/$sequencename.tsv
samtools coverage $bamsorted >  $dir/analysis/$d/$sequencename.coverage
samtools stats $bamsorted | grep ^SN | cut -f 2- | grep 'average length:' >  $dir/analysis/$d/$sequencename.readlength
awk -v sample=$sequencename '{$1=sample ; print;}' $dir/analysis/$d/$sequencename.tsv > $dir/analysis/coverage/coverage_file_$sequencename.tsv
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/$d/ $dir/analysis/$d/$sequencename.tsv $sequencename



echo $ref
done
done
cat $dir/analysis/coverage/coverage_file_*gene_L.tsv > $dir/analysis/coverage/depth_file_gene_L.tsv 
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/coverage $dir/analysis/coverage/depth_file_gene_L.tsv  gene_L

cat $dir/analysis/coverage/coverage_file_*gene_M.tsv > $dir/analysis/coverage/depth_file_gene_M.tsv 
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/coverage $dir/analysis/coverage/depth_file_gene_M.tsv  gene_M

cat $dir/analysis/coverage/coverage_file_*gene_S.tsv > $dir/analysis/coverage/depth_file_gene_S.tsv 
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/coverage $dir/analysis/coverage/depth_file_gene_S.tsv  gene_S



echo "###### finished PTV pipeline at $date #########"

elif [ $virus = "adipo" ] ;
then

mkdir -p $dir/analysis
mkdir -p $dir/analysis/coverage

reference="/home/jovyan/shared/data/adinopectin.reference2.fasta"

REF_NAME=`cat $reference | grep '>' | tr -d '>' | cut -d ' ' -f 1`

LENGTH=`tail -n +2 $reference | tr -d '\n' | wc -m | xargs`
echo -e "$REF_NAME\t$LENGTH" > $dir/analysis/my.genome
#cat $dir/analysis/my.genome
for read_1 in $dir/*R1*
do #echo $read_1
read_2=${read_1/R1/R2}
echo $read_1
echo $read_2
longname=$(basename $read_1)
echo $longname

c=$(echo $longname | awk -F "_S" '{print $1}')
    #echo $c
d=$c
#d=$(echo $c | awk -F "/" '{print $7}')
    echo $d
echo "starting $d mapping with $REF_NAME fasta reference provided"
echo "############################################"
now="$(date +%c)"

echo -e "starting $d $REF_NAME Analysis at \t$now" >> "$dir/analysis/mensajes.log"
mkdir -p $dir/analysis/$d
/home/jovyan/shared/data/fastp/fastp -i $read_1 -I $read_2 -o $dir/analysis/$d/$d.R1.fastq.gz -O $dir/analysis/$d/$d.R2.fastq.gz -w 10 -R "fastp report for sample: $d" -h $dir/analysis/$d/$d.fastqreport.html -j $dir/analysis/$d/$d.jsonreport.json -f 25 -F 25
echo "###############finished fastq step for $d reads 1 and 2 #############################"

#gunzip $dir/analysis/$d/$d.R1.fastq.gz $dir/analysis/$d/$d.R2.fastq.gz
echo "###############starting trim primers steps for $d reads 1 and 2 #############################"

#/home/jovyan/shared/data/minimap2-2.24_x64-linux/minimap2 -ax sr /home/jovyan/shared/data/adinopectin.reference2.fasta $dir/analysis/$d/$d.R1.fastq.gz -2 $dir/analysis/$d/$d.R2.fastq.gz  | python /home/jovyan/shared/data/Alt_nCov2019_primers/tools/trim_primers/trim_primer_parts.py /home/jovyan/shared/data/Alt_nCov2019_primers/Primers/adipo_v1/adinopectin.primer.bed $dir/analysis/$d/$d.R1.cut.fastq $dir/analysis/$d/$d.R2.cut.fastq --gzip

#pigz -p 10 $dir/analysis/$d/$d.R1.fastq $dir/analysis/$d/$d.R2.fastq 

echo "Starting minimap mapping"

/home/jovyan/shared/data/minimap2-2.24_x64-linux/minimap2 -ax sr /home/jovyan/shared/data/adinopectin.reference2.fasta $dir/analysis/$d/$d.R1.fastq.gz -2 $dir/analysis/$d/$d.R2.fastq.gz  > $dir/analysis/$d/assembled.sam 


samtools view  -bS $dir/analysis/$d/assembled.sam  > $dir/analysis/$d/assembled1.bam

samtools sort $dir/analysis/$d/assembled1.bam -o $dir/analysis/$d/assembled.bam

#$input_dir/assembled.sam $input_dir/assembled1.bam

bamsorted=$dir/analysis/$d/assembled.bam

bedtools bamtobed -i $dir/analysis/$d/assembled.bam > $dir/analysis/$d/reads.bed
bedtools genomecov -bga -i $dir/analysis/$d/reads.bed -g $dir/analysis/my.genome | awk  '$4 < 3' > $dir/analysis/$d/zero.bed
echo '#####bamtobedfinished###'
echo '#####maskingbed###'

maskFastaFromBed -fi /home/jovyan/shared/data/adinopectin.reference2.fasta -bed $dir/analysis/$d/zero.bed -fo $dir/analysis/$d/masked.fasta

echo '####masked reference ready###'
bcftools mpileup -Ou  -f $dir/analysis/$d/masked.fasta  $dir/analysis/$d/assembled.bam | bcftools call --ploidy 1 -mv -Oz -o $dir/analysis/$d/test.vcf.gz
echo "###Mpileup ready###"
bcftools index $dir/analysis/$d/test.vcf.gz
cat $dir/analysis/$d/masked.fasta | bcftools consensus $dir/analysis/$d/test.vcf.gz > $dir/analysis/$d/new_consensus.fasta

sequencename=$(echo $d'_'$REF_NAME)

echo ">$sequencename" > $dir/analysis/$d/$d.fas
echo ">$d" > $dir/analysis/$d/$d.code.fas

tail -n +2 $dir/analysis/$d/new_consensus.fasta >> $dir/analysis/$d/$d.fas
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $dir/analysis/$d/$d.fas | tr "\t" "\n" | sed -r '/^>/! s/N+$|^N+//g' |  fold -w 80 >  $dir/analysis/$d/$d.Ntrimmed.fas

tail -n +2 $dir/analysis/$d/new_consensus.fasta >> $dir/analysis/$d/$d.code.fas
awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $dir/analysis/$d/$d.code.fas | tr "\t" "\n" | sed -r '/^>/! s/N+$|^N+//g' |  fold -w 80 >  $dir/analysis/$d/$d.Ncodetrimmed.fas

echo "############################################"
echo "finishing $d consensus generation with bcftools mpileup"
echo "############################################"

samtools depth -aa $bamsorted | awk -v sample=$d '{$1=sample ; print;}' >  $dir/analysis/$d/$d.tsv
samtools coverage $bamsorted >  $dir/analysis/$d/$d.coverage
samtools stats $bamsorted | grep ^SN | cut -f 2- | grep 'average length:' >  $dir/analysis/$d/$d.readlength
awk -v sample=$d '{$1=sample ; print;}' $dir/analysis/$d/$d.tsv > $dir/analysis/coverage/coverage_file_$d.tsv
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/$d/ $dir/analysis/$d/$d.tsv $d
 

done

cat $dir/analysis/coverage/coverage_file_*.tsv > $dir/analysis/coverage/depth_file.tsv 
python /home/jovyan/shared/data/plotgroupedcoveragev2.py $dir/analysis/coverage $dir/analysis/coverage/depth_file.tsv $REF_NAME

echo "###### finished adipo pipeline at $date #########"



elif [ $virus = "mpxv" ] ;
then

echo "###### starting mpxv pipeline at $date #########"
mkdir -p $dir/analysis
mkdir -p $dir/analysis/coverage

reference="/home/jovyan/shared/data/NC_063383.1.ref.fasta"

REF_NAME=`cat $reference | grep '>' | tr -d '>' | cut -d ' ' -f 1`

LENGTH=`tail -n +2 $reference | tr -d '\n' | wc -m | xargs`
echo -e "$REF_NAME\t$LENGTH" > $dir/analysis/my.genome
#cat $dir/analysis/my.genome
for read_1 in $dir/*R1*
do #echo $read_1
read_2=${read_1/R1/R2}
echo $read_1
echo $read_2
longname=$(basename $read_1)
echo $longname

c=$(echo $longname | awk -F "_S" '{print $1}')
    #echo $c
d=$c
#d=$(echo $c | awk -F "/" '{print $7}')
    echo $d
echo "starting $d mapping with $REF_NAME fasta reference provided"
echo "############################################"
now="$(date +%c)"

echo -e "starting $d $REF_NAME Analysis at \t$now" >> "$dir/analysis/mensajes.log"
mkdir -p $dir/analysis/$d
/home/jovyan/shared/data/fastp/fastp -i $read_1 -I $read_2 -o $dir/analysis/$d/$d.R1.fastq.gz -O $dir/analysis/$d/$d.R2.fastq.gz -w 10 -R "fastp report for sample: $d" -h $dir/analysis/$d/$d.fastqreport.html -j $dir/analysis/$d/$d.jsonreport.json -f 25 -F 25
echo "###############finished fastq step for $d reads 1 and 2 #############################"


echo "Please update the MPXV pipeline for adding the reference and the loop for analyzing illumina generated reads."

done


elif [ $virus = "discovery" ] ;
then
start=`date +%m`

now="$(date +%c)"
echo "###### starting pathogen discovery pipeline at $now #########"
#mkdir -p $dir/analysis
#mkdir -p $dir/analysis/coverage
echo "Will create a folder $dir/analysis"
N=4
pathogen_dis(){
   /home/jovyan/shared/data/pathogen_discovery_adapted.sh \
   --sample_num $1 \
   --read_type "illumina-paired" \
   --input_dir $dir \
   --kraken_db "/home/jovyan/databases/kraken2/kraken2_human_db" \
   --read_1 $2 \
   --output_dir $dir/pathogen_discovery_results
}

(
for read_1 in $dir/*R1*;
do
    start=`date +%s`
   ((i=i%N)); ((i++==0)) && wait
   longname=$(basename $read_1)
   c=$(echo $longname | awk -F "_S" '{print $1}')
   
   
   pathogen_dis $c $longname & 
done
wait
)
end="$(date +%c)"
end2=`date +%s`

echo "###### pathogen discovery pipeline finished at $now  it took `expr $end2 - $start`  minutes #########"

elif [ $virus = "discoveryviral" ] ;
then
start=`date +%m`

now="$(date +%c)"
echo "###### starting pathogen discovery pipeline at $now #########"
#mkdir -p $dir/analysis
#mkdir -p $dir/analysis/coverage
echo "Will create a folder $dir/analysis"
N=4
pathogen_dis(){
   /home/jovyan/shared/data/pathogen_discovery_viral.sh \
   --sample_num $1 \
   --read_type "illumina-paired" \
   --input_dir $dir \
   --kraken_db "/home/jovyan/databases/kraken2/kraken2_human_db" \
   --read_1 $2 \
   --output_dir $dir/pathogen_viral_results
}

(
for read_1 in $dir/*R1*;
do
    start=`date +%s`
   ((i=i%N)); ((i++==0)) && wait
   longname=$(basename $read_1)
   c=$(echo $longname | awk -F "_S" '{print $1}')
   
   
   pathogen_dis $c $longname & 
done
wait
)
end="$(date +%c)"
end2=`date +%s`

echo "###### pathogen discovery pipeline finished at $now  it took `expr $end2 - $start`  minutes #########"



elif [ $virus = "discoverynano" ] ;
then
start=`date +%m`

now="$(date +%c)"
echo "###### starting pathogen discovery pipeline at $now #########"
#mkdir -p $dir/analysis
#mkdir -p $dir/analysis/coverage
echo "Will create a folder $dir/analysis"
N=4
pathogen_dis(){
   /home/jovyan/shared/data/pathogen_discovery_adapted.sh \
   --sample_num $1 \
   --read_type "nanopore" \
   --input_dir $dir \
   --kraken_db "/home/jovyan/databases/kraken2/kraken2_human_db" \
   --read_1 $2 \
   --output_dir $dir/pathogen_discovery_results
}

(
for read_1 in $dir/*gz;
do
    start=`date +%s`
   ((i=i%N)); ((i++==0)) && wait
   longname=$(basename $read_1)
   c=$(echo $longname | awk -F ".fastq" '{print $1}')
   
   
   pathogen_dis $c $longname & 
done
wait
)
end="$(date +%c)"
end2=`date +%s`

echo "###### pathogen discovery pipeline finished at $now  it took `expr $end2 - $start`  seconds #########"


elif [ $virus = "INF_illumina" ] ;
then

mkdir -p $dir/analysis
mkdir -p $dir/analysis/fastq
mkdir -p $dir/analysis/coverage

for read_1 in $dir/*R1*
do #echo $read_1
read_2=${read_1/R1/R2}
echo $read_1
echo $read_2
longname=$(basename $read_1)
echo $longname

c=$(echo $longname | awk -F "_S" '{print $1}')
    echo $c
#d=$(echo $c | awk -F "/" '{print $7}')
    echo $d
d=$c
echo "starting $d mapping with $REF_NAME fasta reference provided"
echo "############################################"
now="$(date +%c)"

echo -e "starting $d $REF_NAME Analysis at \t$now" >> "$dir/analysis/mensajes.log"
#mkdir -p $dir/analysis/$d
/home/jovyan/shared/data/fastp/fastp -i $read_1 -I $read_2 -o $dir/analysis/fastq/$d.R1.fastq.gz -O $dir/analysis/fastq/$d.R2.fastq.gz -w 10 -R "fastp report for sample: $d" -h $dir/analysis/fastq/$d.fastqreport.html -j $dir/analysis/fastq/$d.jsonreport.json
echo "###############finished fastq step for $d reads 1 and 2 #############################"

/home/jovyan/shared/data/flu-amd/IRMA FLU $dir/analysis/fastq/$d.R1.fastq.gz $dir/analysis/fastq/$d.R2.fastq.gz $d $dir/analysis

echo "###############finished IRMA pipeline for $d reads 1 and 2 #############################"

done

elif [ $virus = "INF_ONT" ] ;
then

mkdir -p $dir/allresults
mkdir -p $dir/allresults/fastas
mkdir -p $dir/allresults/coverages


for folder in $dir/barcode* 
do #echo $read_1
#read_2=${read_1/R1/R2}
#echo $read_1
#echo $read_2
#echo $folder
longname=$(basename $folder)
#echo  $longname

python /home/jovyan/shared/data/diclist.py $meta $longname > $dir/output
sample=$(cat $dir/output)
if [ $sample = "barcodenotfound" ]; then
echo "#### $longname not found in metadata"
else 
echo "using $longname for sample $sample in the analysis"
gunzip $dir/$longname/*gz
mkdir $dir/$longname/merged
cat $dir/$longname/*fastq > $dir/$longname/merged/$sample.$longname.fastq
pigz $dir/$longname/merged/$sample.$longname.fastq
/home/jovyan/shared/data/flu-amd/IRMA FLU-minion  $dir/$longname/merged/$sample.$longname.fastq.gz $dir/$longname/$sample
Rscript /home/jovyan/shared/data/coverage4allgenes.R $dir/$longname/$sample/tables $sample  $dir/allresults/coverages
cat $dir/$longname/$sample/amended_consensus/*_1.fa >> $dir/allresults/fastas/PB2_all.fasta
cat $dir/$longname/$sample/amended_consensus/*_2.fa >> $dir/allresults/fastas/PB1_all.fasta
cat $dir/$longname/$sample/amended_consensus/*_3.fa >> $dir/allresults/fastas/PA_all.fasta
cat $dir/$longname/$sample/amended_consensus/*_4.fa >> $dir/allresults/fastas/HA_all.fasta
cat $dir/$longname/$sample/amended_consensus/*_5.fa >> $dir/allresults/fastas/NP_all.fasta
cat $dir/$longname/$sample/amended_consensus/*_6.fa >> $dir/allresults/fastas/NA_all.fasta
cat $dir/$longname/$sample/amended_consensus/*_7.fa >> $dir/allresults/fastas/M_all.fasta
cat $dir/$longname/$sample/amended_consensus/*_8.fa >> $dir/allresults/fastas/NS_all.fasta
mkdir $dir/$longname/rawfiles
mv $dir/$longname/*fastq $dir/$longname/rawfiles/.
pigz $dir/$longname/rawfiles/*fastq
echo "###############finished IRMA pipeline for $longname of sample $sample #############################"
fi 


rm $dir/output
done

elif [ $virus = "DENV" ] ;
then
start=`date +%m`

now="$(date +%c)"
echo "###### starting Dengue pipeline at $now #########"
#mkdir -p $dir/analysis
#mkdir -p $dir/analysis/coverage
echo "Will use folder $dir for read analysis"
N=4
DENVpipe(){
   /home/jovyan/shared/data/denvpiperun.sh $1
}

DENVpipe $dir

end="$(date +%c)"
end2=`date +%s`

echo "###### Dengue pipeline finished at $now  it took `expr $end2 - $start`  minutes #########"


echo "############### Finished gencomvariable.sh pipeline for organism $accession $virus in directory $dir #############################"
fi




