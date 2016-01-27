#!/bin/sh
#Created by: Chan Chee Seng
#Version: 1.1
#Revision history: 10 April 2013
#                  16 August 2013 added parsing of the config_casava-1.8.2.txt to extract the lanes to map
#                  5 June 2014 : fork pipeline for Single Cell RNAseq pipeline
#Description: This script prepares the directory for STAR mapping pipeline with the run directory.  Run this script after the CASAVA Unalign step and must run within the UnalignedDir directory
#Usage: Usage: $0 runidWithoutFlowcelID UnalignedDir jobID

source /home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.config

#dollar 1 is run number without flowcell ID
#dollar 2 is absolute path to SampleSheet.csv
#dollar 3 is relative path to Unaligned directory
#dollar 4 is job ID
#-r run directory name including flowcell id
#-s is absolute path to SampleSheet.csv
#-f is relative path to Unaligned directory
#-j is job ID
#-p is whether consumer is Production or not. Optional, when ad-hoc analysis do not use this parameter
#-c is elm logger status code: 5,6,7,10,11,12,13
#-C is user comments for annotating reruns
#-P if Pre-processing the fastq reads

#main
#definition of all options variables
runDirName="undef"
sampleSheetPath="undef"
UnalignedDir="undef"
elmLoggerStatusStart="undef"
jobId="undef"
consumer="undef"
comments=""
userCustomConfigFile="undef"
userCustomSTARConfigFile="undef"
userQueueName="undef"
#Preprocessing="undef"
#if ( ! getopts "c:f:j:p:r:s:u:U:q:C:P:" FLAG )
if ( ! getopts "c:f:j:p:r:s:u:U:q:C:" FLAG )
then
  echo [`date`] [`hostname`] [$0] No options provided.
  echo [`date`] [`hostname`] [$0] Usage: $0 generateStarExpressionconfiguration.sh -r run directory name -s sampleSheet.csv -f UnalignedDir -u userCustomConfigFile -U userCustomSTARConfigFile -j jobID -p "Production|Community" -q "All users others than userrig must provide queue name excluding rnaseqpipeline.q" -c elmLoggerStatusCode
#  echo [`date`] [`hostname`] [$0] Usage: $0 generateStarExpressionconfiguration.sh -r run directory name -s sampleSheet.csv -f UnalignedDir -u userCustomConfigFile -U userCustomSTARConfigFile -j jobID -p "Production|Community" -q "All users others than userrig must provide queue name excluding rnaseqpipeline.q" -c elmLoggerStatusCode -P Pre-processing fastq reads
  exit 1
fi
#while getopts "c:f:j:p:r:s:u:U:q:C:P:" FLAG
while getopts "c:f:j:p:r:s:u:U:q:C:" FLAG
do
  case $FLAG in
  c) elmLoggerStatusStart=${OPTARG}
    elmLoggerStatusComplete=$((elmLoggerStatusStart+1))
    elmLoggerStatusCode=${elmLoggerStatusStart}
    ;;
  f) UnalignedDir=${OPTARG}
    ;;
  j) jobId=${OPTARG}
    ;;
  p) consumer=${OPTARG}
    ;;
  r) runDirName=${OPTARG}
    ;;
  s) sampleSheetPath=${OPTARG}
    sampleSheetPathDirOnly=$(dirname ${sampleSheetPath})
    ;;
  u) userCustomConfigFile=${OPTARG}
    ;;
  U) userCustomSTARConfigFile=${OPTARG}
    ;;
  q) userQueueName=${OPTARG}
    ;;
  C) comments=${OPTARG}
    ;;
# P) Preprocessing=${OPTARG}
#   ;; 
  *) echo [`date`] [`hostname`] [$0] Error!! Invalid option >&2
    exit 1
    ;;
  esac
done
#if any of the essential options are not provided, flag error and exit with failure
if [ ${elmLoggerStatusStart} == "undef" ] || [ ${UnalignedDir} == "undef" ] || [ ${consumer} == "undef" ] || [ ${runDirName} == "undef" ] || [ ${sampleSheetPath} == "undef" ] || [ ${consumer} == "undef" ]
then
  echo [`date`] [`hostname`] [$0] Not all options provided.
  echo [`date`] [`hostname`] [$0] Usage: $0 -r run directory name -s sampleSheet.csv -f UnalignedDir -u userCustomConfigFile -U userCustomSTARConfigFile -j jobID -p "Production|Community" -q userQueueName -c elmLoggerStatusCode
#  echo [`date`] [`hostname`] [$0] Usage: $0 -r run directory name -s sampleSheet.csv -f UnalignedDir -u userCustomConfigFile -U userCustomSTARConfigFile -j jobID -p "Production|Community" -q userQueueName -c elmLoggerStatusCode -P Preprocessing
  exit 1
fi

echo [`date`] [`hostname`] [$0] parameters received: elmLoggerStatusStart=${elmLoggerStatusStart} UnalignedDir=${UnalignedDir} userCustomConfigFile=${userCustomConfigFile} userCustomSTARConfigFile=${userCustomSTARConfigFile} consumer=${consumer} runDirName=${runDirName} sampleSheetPath=${sampleSheetPath} userQueueName=${userQueueName}
echo [`date`] [`hostname`] [$0] All parameters=$@
echo [`date`] [`hostname`] [$0] Additional parameters set up to this point: elmLoggerStatusComplete=${elmLoggerStatusComplete} elmLoggerStatusCode=${elmLoggerStatusCode}


#If customConfigFile is defined and not null value, source user provided parameters
if [ ${userCustomConfigFile} == "undef" ]
then
  echo [`date`] [`hostname`] ${userCustomConfigFile}="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionDefaultParameters.config"
  echo [`date`] [`hostname`] ${userCustomConfigFile}=${userCustomConfigFile}
  ${userCustomConfigFile}="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionDefaultParameters.config"
  if [ ${userCustomSTARConfigFile} == "undef" ]
  then 
  echo [`date`] [`hostname`] ${userCustomSTARConfigFile}="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/Parameters1.in"
  echo [`date`] [`hostname`] ${userCustomSTARConfigFile}=${userCustomSTARConfigFile}
  ${userCustomSTARConfigFile}="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/Parameters1.in"
  else
       echo [`date`] [`hostname`] ${userCustomSTARConfigFile}
  fi
else
  userCustomConfigFile=${userCustomConfigFile}
  echo [`date`] [`hostname`] ${userCustomConfigFile}
  if [ ${userCustomSTARConfigFile} == "undef" ]
  then 
  echo [`date`] [`hostname`] ${userCustomSTARConfigFile}="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/Parameters1.in"
  echo [`date`] [`hostname`] ${userCustomSTARConfigFile}=${userCustomSTARConfigFile}
  ${userCustomSTARConfigFile}="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/Parameters1.in"
  else
       echo [`date`] [`hostname`] ${userCustomSTARConfigFile}
  fi
fi

####if [ $# -lt 4 ]
####then
####   echo [`date`] [`hostname`] [$0] Usage: $0 runidWithoutFlowcelID SampleSheet.csv UnalignedDir jobID
####   exit 1
####fi

####runId=${1}
####sampleSheetPath=${2}
runId=$(echo ${runDirName}|cut -f1 -d"_")
sampleSheetName=$(basename ${sampleSheetPath})
indexType=$(echo ${sampleSheetName}|cut -f2 -d"_")
machId=$(echo ${runId}|cut -f1 -d"-")
####UnalignedDir=${3}

###=====================================================================
#test if variable is set to null string or set to some value
if [ -z "${consumer:+zzz}" ]
then
  #variable is not set
  echo [`date`] [`hostname`] [$0] "consumer is not set; consumer set to Community"
  pipelineConsumer="Community"
elif [ ${consumer} == "Production" ]
then
  #variable is set to Production
  echo [`date`] [`hostname`] [$0] consumer is Production
  pipelineConsumer="Production"
elif [ ${consumer} == "Community" ]
then
  #variable is set to Community
  echo [`date`] [`hostname`] [$0] consumer is Community
  pipelineConsumer="Community"
else
  #variable is set to invalid value
  echo [`date`] [`hostname`] [$0] invalid value for consumer. consumer is set to consumer=${consumer}. Program exiting...
  exit 1
fi
###=====================================================================

todayDate=`date +%Y%m%dT%H%M%S`
#Get the full path to the current UnalignedDir directory
runDirPath=$(pwd -L)
#return to run directory
cd ..

touch runParameters.txt

echo [`date`] [`hostname`] [$0] generateStarExpressionconfiguration.sh -r ${runDirName} -s ${sampleSheetPath} -f ${UnalignedDir} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -j ${jobId} -p ${pipelineConsumer} -c ${elmLoggerStatusCode} ${userQueueName:+-q "$userQueueName"} >>runParameters.txt

#Check if SampleSheet.csv exists here
if [ ! -f ${sampleSheetPath} ]
then
  echo [`date`] [`hostname`] [$0] ${sampleSheetPath} not found in this directory: $(pwd). Program exiting...
  exit 1
fi

#Check if the UnalignedDir exists here
if [ -d ${UnalignedDir} ]
then
  exptDir=${UnalignedDir}
  outputDir=$(echo ${exptDir}|sed "s/\/.*$//")
else
  echo [`date`] [`hostname`] [$0] ${UnalignedDir} not found in this directory
  exit 1
fi

if [ -f ${sampleSheetPathDirOnly}/config_casava-1.8.2.txt ]
then
##  laneToMap=$(grep -w "ANALYSIS star_expression" ${sampleSheetPathDirOnly}/config_casava-1.8.2.txt|grep -v "ANALYSIS none"|grep -v "ANALYSIS star_expression_stranded"|cut -f1 -d":")

  laneToMap=$(grep -w "ANALYSIS star_expression" ${sampleSheetPathDirOnly}/config_casava-1.8.2.txt|grep -v "ANALYSIS none"|cut -f1 -d":")
  laneToMapstrand=$(grep -w "ANALYSIS star_expression_stranded" ${sampleSheetPathDirOnly}/config_casava-1.8.2.txt|grep -v "ANALYSIS none"|cut -f1 -d":")
  laneToMapVariantcalling=$(grep -w "ANALYSIS star_variant" ${sampleSheetPathDirOnly}/config_casava-1.8.2.txt|grep -v "ANALYSIS none"|cut -f1 -d":")
  laneToMapVariantcallingstrand=$(grep -w "ANALYSIS star_variant_stranded" ${sampleSheetPathDirOnly}/config_casava-1.8.2.txt|grep -v "ANALYSIS none"|cut -f1 -d":")
else
  echo [`date`] [`hostname`] [$0] ${sampleSheetPathDirOnly}/config_casava-1.8.2.txt is not found in this run directory
  exit 1
fi

stranded="-S Strandspecific"
unstranded="-S NonStrandspecific"

# for lanes map with star_expression
set -o pipefail
if [ ! -z "${laneToMap:+yyy}" ]
then
  mkdir output_${runId}_StarExpression_${todayDate}_${indexType}
  touch runParameters.txt
  cd output_${runId}_StarExpression_${todayDate}_${indexType}
  touch Jobstatus.txt
  echo [`date`] [`hostname`] [$0] 'Mapping these lane(s):' ${laneToMap}
  echo [`date`] [`hostname`] [$0] ${machId}
  if [ ${machId} == "HS006" ] || [ ${machId} == "HS007" ] || [ ${machId} == "NS001" ]
  then
  echo [`date`] [`hostname`] [$0] ${laneToMap}
  /home/userrig/pipelines/whatApps/whatsinthisrunReverseI2.sh -r ${runId}|grep -E "${runId}	[${laneToMap}]"|sort>config_star_expression.txt.tmp
  else
  /home/userrig/pipelines/whatApps/whatsinthisrun.pl ${runId}|grep -E "${runId}	[${laneToMap}]"|sort>config_star_expression.txt.tmp
  fi

  for eachLibrary in $(sed '0,/^Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description$/d' ../${sampleSheetName})
  do
    laneId=$(echo ${eachLibrary}|cut -f1 -d",")
    libraryId=$(echo ${eachLibrary}|cut -f2 -d","|cut -f2 -d"_")
#    barcode=$(echo ${eachLibrary}|cut -f3 -d","|cut -f2 -d"_")
    barcode=$(echo ${eachLibrary}|cut -f3 -d","|awk ' BEGIN {FS = "-"}; {if ($3 == "") print $2; else print $2"-"$3}')
    if [ ! -z ${barcode} ]
    then
    if [ ${barcode} == "NoIndex" ]
    then
      awk -v laneId=${laneId} -v libraryId=${libraryId} '{if ($2==laneId && $3==libraryId) print $0}' config_star_expression.txt.tmp >>config_star_expression.txt
    else
      awk -v laneId=${laneId} -v libraryId=${libraryId} -v barcode=${barcode} '{if ($2==laneId && $5==barcode && $7==libraryId) print $0}' config_star_expression.txt.tmp >>config_star_expression.txt
    fi
    else
    echo "skip barcode"
    fi
  done
  if [ ! -s config_star_expression.txt ]
  then
    echo [`date`] [`hostname`] [$0] config_star_expression.txt is EMPTY! Skipping star mapping...
    printf "[`date`] [`hostname`] [$0] [`pwd`] changing back to previous directory:"
    cd -
    ####exit 1
  else

    echo [`date`] [`hostname`] [$0] Please proceed with the next command "cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${userQueueName:+-q "$userQueueName"} ${unstranded}"
    cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${comments:+-C "$comments"} ${userQueueName:+-q "$userQueueName"} ${unstranded} >RNAseqSTARExpressionPipeline_${todayDate}.log
    echo [`date`] [`hostname`] [$0] "cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${userQueueName:+-q "$userQueueName"} ${unstranded}" >>../runParameters.txt
    bwaMappingExitStatus=$?
    if [ ${bwaMappingExitStatus} -eq 0 ]
    then
      JobIdWaitStarCompleted=$(qsub -q bwapipeline.q -cwd -terse ${DebugSwitch} ${pipelineLoc}/waitOnStarExpressionPipelineCompletedEvent.sh)
      echo [`date`] [`hostname`] [$0] [`pwd`] JobIdWaitStarCompleted=${JobIdWaitStarCompleted}
    fi
    #Now we can call the next pipeline
count=`cat config_star_expression.txt |wc -l`
echo ${count} Number of Libraries in this Job >Jobstatus.txt

    printf "[`date`] [`hostname`] [$0] [`pwd`] changing back to previous directory:"
    cd -
  fi
else
  echo [`date`] [`hostname`] [$0] No star_expression analysis detected. Proceeding to process star_varinat...
fi

#################################################### for lanes to map with star_expression_stranded
if [ ! -z "${laneToMapstrand:+yyy}" ]
then
  mkdir output_${runId}_StarExpressionStranded_${todayDate}_${indexType}
  touch runParameters.txt
  cd output_${runId}_StarExpressionStranded_${todayDate}_${indexType}
  touch Jobstatus.txt
  echo [`date`] [`hostname`] [$0] 'Mapping these lane(s):' ${laneToMapstrand}
  echo [`date`] [`hostname`] [$0] ${machId}
  if [ ${machId} == "HS006" ] || [ ${machId} == "HS007" ] || [ ${machId} == "NS001" ]
  then
  echo [`date`] [`hostname`] [$0] ${laneToMapstrand}
  /home/userrig/pipelines/whatApps/whatsinthisrunReverseI2.sh -r ${runId}|grep -E "${runId}	[${laneToMapstrand}]"|sort>config_star_expression.txt.tmp
  else
  /home/userrig/pipelines/whatApps/whatsinthisrun.pl ${runId}|grep -E "${runId}	[${laneToMapstrand}]"|sort>config_star_expression.txt.tmp
  fi

  for eachLibrary in $(sed '0,/^Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description$/d' ../${sampleSheetName})
  do
    laneId=$(echo ${eachLibrary}|cut -f1 -d",")
    libraryId=$(echo ${eachLibrary}|cut -f2 -d","|cut -f2 -d"_")
#    barcode=$(echo ${eachLibrary}|cut -f3 -d","|cut -f2 -d"_")
    barcode=$(echo ${eachLibrary}|cut -f3 -d","|awk ' BEGIN {FS = "-"}; {if ($3 == "") print $2; else print $2"-"$3}')
    if [ ! -z ${barcode} ]
    then
    if [ ${barcode} == "NoIndex" ]
    then
      awk -v laneId=${laneId} -v libraryId=${libraryId} '{if ($2==laneId && $3==libraryId) print $0}' config_star_expression.txt.tmp >>config_star_expression.txt
    else
      awk -v laneId=${laneId} -v libraryId=${libraryId} -v barcode=${barcode} '{if ($2==laneId && $5==barcode && $7==libraryId) print $0}' config_star_expression.txt.tmp >>config_star_expression.txt
    fi
    else
    echo "skip barcode"
    fi
  done
  if [ ! -s config_star_expression.txt ]
  then
    echo [`date`] [`hostname`] [$0] config_star_expression.txt is EMPTY! Skipping star mapping...
    printf "[`date`] [`hostname`] [$0] [`pwd`] changing back to previous directory:"
    cd -
    ####exit 1
  else

    echo [`date`] [`hostname`] [$0] Please proceed with the next command "cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${userQueueName:+-q "$userQueueName"} ${stranded}"
    cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${comments:+-C "$comments"} ${userQueueName:+-q "$userQueueName"} ${stranded} >RNAseqSTARExpressionPipeline_${todayDate}.log
   echo [`date`] [`hostname`] [$0] "cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${userQueueName:+-q "$userQueueName"} ${stranded}" >>../runParameters.txt

    bwaMappingExitStatus=$?
    if [ ${bwaMappingExitStatus} -eq 0 ]
    then
      JobIdWaitStarCompleted=$(qsub -q bwapipeline.q -cwd -terse ${DebugSwitch} ${pipelineLoc}/waitOnStarExpressionPipelineCompletedEvent.sh)
      echo [`date`] [`hostname`] [$0] [`pwd`] JobIdWaitStarCompleted=${JobIdWaitStarCompleted}
    fi
    #Now we can call the next pipeline
count=`cat config_star_expression.txt |wc -l`
echo ${count} Number of Libraries in this Job >Jobstatus.txt

    printf "[`date`] [`hostname`] [$0] [`pwd`] changing back to previous directory:"
    cd -
  fi
else
  echo [`date`] [`hostname`] [$0] No star_expression_stranded analysis detected. Proceeding to process star_varinat...
fi

##############################################################for lanes to map with star_variant
if [ ! -z "${laneToMapVariantcalling:+yyy}" ]
then
  mkdir output_${runId}_StarVariant_${todayDate}_${indexType}
  touch runParameters.txt
  cd output_${runId}__StarVariant_${todayDate}_${indexType}
  touch Jobstatus.txt
  echo [`date`] [`hostname`] [$0] 'Mapping these lane(s):' ${laneToMapVariantcalling}
  echo [`date`] [`hostname`] [$0] ${machId}
  if [ ${machId} == "HS006" ] || [ ${machId} == "HS007" ] || [ ${machId} == "NS001" ]
  then
  echo [`date`] [`hostname`] [$0] ${laneToMapVariantcalling}
  /home/userrig/pipelines/whatApps/whatsinthisrunReverseI2.sh -r ${runId}|grep -E "${runId}	[${laneToMapVariantcalling}]"|sort>config_star_expression.txt.tmp
  else
  /home/userrig/pipelines/whatApps/whatsinthisrun.pl ${runId}|grep -E "${runId}	[${laneToMapVariantcalling}]"|sort>config_star_expression.txt.tmp
  fi

  for eachLibrary in $(sed '0,/^Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description$/d' ../${sampleSheetName})
  do
    laneId=$(echo ${eachLibrary}|cut -f1 -d",")
    libraryId=$(echo ${eachLibrary}|cut -f2 -d","|cut -f2 -d"_")
#    barcode=$(echo ${eachLibrary}|cut -f3 -d","|cut -f2 -d"_")
    barcode=$(echo ${eachLibrary}|cut -f3 -d","|awk ' BEGIN {FS = "-"}; {if ($3 == "") print $2; else print $2"-"$3}')
    if [ ! -z ${barcode} ]
    then
    if [ ${barcode} == "NoIndex" ]
    then
      awk -v laneId=${laneId} -v libraryId=${libraryId} '{if ($2==laneId && $3==libraryId) print $0}' config_star_expression.txt.tmp >>config_star_expression.txt
    else
      awk -v laneId=${laneId} -v libraryId=${libraryId} -v barcode=${barcode} '{if ($2==laneId && $5==barcode && $7==libraryId) print $0}' config_star_expression.txt.tmp >>config_star_expression.txt
    fi
    else
    echo "skip barcode"
    fi
  done
  if [ ! -s config_star_expression.txt ]
  then
    echo [`date`] [`hostname`] [$0] config_star_expression.txt is EMPTY! Skipping star mapping...
    ####exit 1
  else
  echo [`date`] [`hostname`] [$0] Please proceed with the next command "cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${userQueueName:+-q "$userQueueName"} ${unstranded}"
  cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${comments:+-C "$comments"} ${userQueueName:+-q "$userQueueName"} ${unstranded} >RNAseqStarVariantPipeline_${todayDate}.log
  echo [`date +%s.%N`] [`hostname`] [$0] "cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart}" ${userQueueName:+-q "$userQueueName"} ${unstranded}" >>../runParameters.txt
    bwaMappingExitStatus=$?
    if [ ${bwaMappingExitStatus} -eq 0 ]
    then
      JobIdWaitStarCompleted=$(qsub -q bwapipeline.q -cwd -terse ${DebugSwitch} ${pipelineLoc}/waitOnStarExpressionPipelineCompletedEvent.sh)
      echo [`date`] [`hostname`] [$0] JobIdWaitStarCompleted=${JobIdWaitStarCompleted}

      echo [`date`] [`hostname`] [$0] SNV analysis is requested after WTS. Submitting SNV pipeline
      echo [`date`] [`hostname`] [$0] [`pwd`] Call SNV pipeline from WTS. commandLine="qsub -q bwapipeline.q -cwd -terse -v PATH ${DebugSwitch} -N generateSNVconfiguration -hold_jid ${JobIdWaitStarCompleted} /mnt/AnalysisPool/libraries/SNV_pipeline/bin/generateSNVconfiguration.pl --config_file config_star_expression.txt --BWA_bam_output_folder $(pwd) --SNV_output_folder $(pwd) --PCR_status yes --type dna --elmLoggerStatus 5 --consumer production"

      JobIdGenerateSnvConfig=$(qsub -q bwapipeline.q -cwd -terse -v PATH ${DebugSwitch} -N generateSNVconfiguration -hold_jid ${JobIdWaitStarCompleted} /mnt/AnalysisPool/libraries/SNV_pipeline/bin/generateSNVconfiguration.pl --config_file config_star_expression.txt --BWA_bam_output_folder $(pwd) --SNV_output_folder $(pwd) --PCR_status yes --type dna --elmLoggerStatus 5 --consumer production)
      echo [`date`] [`hostname`] [$0] [`pwd`] generateSNVconfiguration.pl submitted with jobId JobIdGenerateSnvConfig=${JobIdGenerateSnvConfig}
    fi
    #Now we can call the next pipeline
count=`cat config_star_expression.txt |wc -l`
echo ${count} Number of Libraries in this Job >Jobstatus.txt

    printf "[`date`] [`hostname`] [$0] [`pwd`] changing back to previous directory:"
    cd -
  fi
else
  echo [`date`] [`hostname`] [$0] No star_variant analysis detected. Program exiting...
fi

##############################################################for lanes to map with star_variant
  if [ ! -z "${laneToMapVariantcallingstrand:+yyy}" ]
  then
  mkdir output_${runId}_StarVariantStranded_${todayDate}_${indexType}
  touch runParameters.txt
  cd output_${runId}__StarVariantStranded_${todayDate}_${indexType}
  touch Jobstatus.txt
  echo [`date`] [`hostname`] [$0] 'Mapping these lane(s):' ${laneToMapVariantcallingstrand}
  echo [`date`] [`hostname`] [$0] ${machId}
  if [ ${machId} == "HS006" ] || [ ${machId} == "HS007" ] || [ ${machId} == "NS001" ]
  then
  echo [`date`] [`hostname`] [$0] ${laneToMapVariantcallingstrand}
  /home/userrig/pipelines/whatApps/whatsinthisrunReverseI2.sh -r ${runId}|grep -E "${runId}	[${laneToMapVariantcallingstrand}]"|sort>config_star_expression.txt.tmp
  else
  /home/userrig/pipelines/whatApps/whatsinthisrun.pl ${runId}|grep -E "${runId}	[${laneToMapVariantcallingstrand}]"|sort>config_star_expression.txt.tmp
  fi

  for eachLibrary in $(sed '0,/^Lane,Sample_ID,Sample_Name,Sample_Plate,Sample_Well,I7_Index_ID,index,I5_Index_ID,index2,Sample_Project,Description$/d' ../${sampleSheetName})
  do
    laneId=$(echo ${eachLibrary}|cut -f1 -d",")
    libraryId=$(echo ${eachLibrary}|cut -f2 -d","|cut -f2 -d"_")
#    barcode=$(echo ${eachLibrary}|cut -f3 -d","|cut -f2 -d"_")
    barcode=$(echo ${eachLibrary}|cut -f3 -d","|awk ' BEGIN {FS = "-"}; {if ($3 == "") print $2; else print $2"-"$3}')
    if [ ! -z ${barcode} ]
    then
    if [ ${barcode} == "NoIndex" ]
    then
      awk -v laneId=${laneId} -v libraryId=${libraryId} '{if ($2==laneId && $3==libraryId) print $0}' config_star_expression.txt.tmp >>config_star_expression.txt
    else
      awk -v laneId=${laneId} -v libraryId=${libraryId} -v barcode=${barcode} '{if ($2==laneId && $5==barcode && $7==libraryId) print $0}' config_star_expression.txt.tmp >>config_star_expression.txt
    fi
    else
    echo "skip barcode"
    fi
  done
  if [ ! -s config_star_expression.txt ]
  then
    echo [`date`] [`hostname`] [$0] config_star_expression.txt is EMPTY! Skipping star mapping...
    ####exit 1
  else
  echo [`date`] [`hostname`] [$0] Please proceed with the next command "cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${userQueueName:+-q "$userQueueName"} ${stranded}"
  cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${comments:+-C "$comments"} ${userQueueName:+-q "$userQueueName"} ${stranded} >RNAseqStarVariantPipeline_${todayDate}.log
  echo [`date +%s.%N`] [`hostname`] [$0] "cat config_star_expression.txt|/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${runDirPath} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -p ${pipelineConsumer} -c ${elmLoggerStatusStart}" ${userQueueName:+-q "$userQueueName"} ${stranded}" >>../runParameters.txt
    bwaMappingExitStatus=$?
    if [ ${bwaMappingExitStatus} -eq 0 ]
    then
      JobIdWaitStarCompleted=$(qsub -q bwapipeline.q -cwd -terse ${DebugSwitch} ${pipelineLoc}/waitOnStarExpressionPipelineCompletedEvent.sh)
      echo [`date`] [`hostname`] [$0] JobIdWaitStarCompleted=${JobIdWaitStarCompleted}

      echo [`date`] [`hostname`] [$0] SNV analysis is requested after WTS. Submitting SNV pipeline
      echo [`date`] [`hostname`] [$0] [`pwd`] Call SNV pipeline from WTS. commandLine="qsub -q bwapipeline.q -cwd -terse -v PATH ${DebugSwitch} -N generateSNVconfiguration -hold_jid ${JobIdWaitStarCompleted} /mnt/AnalysisPool/libraries/SNV_pipeline/bin/generateSNVconfiguration.pl --config_file config_star_expression.txt --BWA_bam_output_folder $(pwd) --SNV_output_folder $(pwd) --PCR_status yes --type dna --elmLoggerStatus 5 --consumer production"

      JobIdGenerateSnvConfig=$(qsub -q bwapipeline.q -cwd -terse -v PATH ${DebugSwitch} -N generateSNVconfiguration -hold_jid ${JobIdWaitStarCompleted} /mnt/AnalysisPool/libraries/SNV_pipeline/bin/generateSNVconfiguration.pl --config_file config_star_expression.txt --BWA_bam_output_folder $(pwd) --SNV_output_folder $(pwd) --PCR_status yes --type dna --elmLoggerStatus 5 --consumer production)
      echo [`date`] [`hostname`] [$0] [`pwd`] generateSNVconfiguration.pl submitted with jobId JobIdGenerateSnvConfig=${JobIdGenerateSnvConfig}
    fi
    #Now we can call the next pipeline
count=`cat config_star_expression.txt |wc -l`
echo ${count} Number of Libraries in this Job >Jobstatus.txt

  fi
else
  echo [`date`] [`hostname`] [$0] No star_variant_stranded analysis detected. Program exiting...
  exit 0
fi
