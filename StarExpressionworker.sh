#!/bin/sh
#Copyright © 2014 Genome Institute of Singapore. All rights reserved.
#Created by: Chan Chee Seng
#Version: 3.0
#Revision history: 14 July 2013
#                  03 September 2013: Change the libraryFileSizeStr to byte format
#                  21 October 2013: Change variable name from "cufflinks_ReferenceSequence" to "cufflinks_FragmentBiasCorrection"
#                  27 November 2013: Updated gis logging messages for tophat completed and pipeLineStage repeated stage names
#                  5 June 2014: fork pipeline for Single Cell RNAseq pipeline
#                  15 July 2014: Introduce Read group info and mask file (rRNA,Mt_tRNA,tRNA ucsc) in cufflinks
#                  1 January 2015: fix bug for processing of cuffdiff output
#                  13 January 2015: Introduced bam index
#Description: tophat and cufflinks
#Usage: StarExpressionworker.sh <full path to fastq directory> <assembly> <identifier>
#Usage: StarExpressionworker.sh <pathToCustomConfigFile>
#-f full path to fastq directory
#-a assembly
#-i identifier
#-p is whether consumer is Production or not. Optional, when ad-hoc analysis do not use this parameter
#-c is elm logger status code: 5,6,7,10,11,12,13
#-u pathToUserCustomConfigFile
#-C is user comments for annotating reruns
#-S is stranded information

# -- our name ---
#$ -N RNAseqSTARExpression
#$ -S /bin/sh

#Start of function definition
usage() {
  echo [`date`] [`hostname`] [$0] Usage: StarExpressionWorker.sh -f full path to fastq.gz file -a assembly -i identifier -p "Production|Community" -c elmLoggerStatusCode -U STARParameterFile -u RNAseq Config File -S Strandspecific NonStrandspecific
}
send_elm_logger_status() {
  if [ ${DEBUG} -eq 0 ]
  then
    # actual parameter list for the logger
    # elm_logger_status_rnaseq_pipeline.pl runId laneNum runType libraryId elm_logger_status elmComments pipeLine_Stage pipeLine_Stage_status tool_version command_used unique_identifier cpu_Count library_file_size server_node

    echo [`date`] [`hostname`] [$0] ${pipelineLoggerCommand} --runId ${runId} --laneNum ${laneNum} --runType ${runType} --libraryId ${libraryId} --elmLoggerStatus ${elmLoggerStatusCode} --elmComments "${elmComments}" --pipeLineName "${pipeLineName}" --pipeLineVersion "${pipeLineVersion}" --pipeLineStage "${pipeLineStage}" --pipeLineStageStatus "${pipeLineStageStatus}" --toolVersion "${toolVersion}" --commandLine "${commandLine}" --identifier ${identifier} --cpuCount ${cpuCount} --libraryFileSize ${libraryFileSizeStr} --consumer Production
                                    ${pipelineLoggerCommand} --runId ${runId} --laneNum ${laneNum} --runType ${runType} --libraryId ${libraryId} --elmLoggerStatus ${elmLoggerStatusCode} --elmComments "${elmComments}" --pipeLineName "${pipeLineName}" --pipeLineVersion "${pipeLineVersion}" --pipeLineStage "${pipeLineStage}" --pipeLineStageStatus "${pipeLineStageStatus}" --toolVersion "${toolVersion}" --commandLine "${commandLine}" --identifier ${identifier} --cpuCount ${cpuCount} --libraryFileSize ${libraryFileSizeStr} --consumer Production
  fi
}

#main
#definition of all options variables
fastqFileDir="undef"
assembly="undef"
identifier="undef"
consumer="undef"
elmLoggerStatusStart="undef"
userCustomConfigFile="undef"
userCustomSTARConfigFile="undef"
comments=""
stranded="undef"

if ( ! getopts "a:c:f:i:p:u:U:C:S:" FLAG )
then
  echo [`date`] [`hostname`] [$0] No options provided.
  usage
  exit 1
fi
while getopts "a:c:f:i:p:u:U:C:S:" FLAG
do
  case $FLAG in
  a) assembly=${OPTARG}
    ;;
  c) elmLoggerStatusStart=${OPTARG}
    elmLoggerStatusComplete=$((elmLoggerStatusStart+1))
    elmLoggerStatusTroubleshooting=7
    ;;
  f) fastqFileDir=${OPTARG}
    ;;
  i) identifier=${OPTARG}
    ;;
  p) consumer=${OPTARG}
    ;;
  u) userCustomConfigFile=${OPTARG}
    customConfigFile=${userCustomConfigFile}
    pipelineConsumer="Community"
    ;;
  U) userCustomSTARConfigFile=${OPTARG}
    customSTARConfigFile=${userCustomSTARConfigFile}
    pipelineConsumer="Community"
    ;;
  C) comments=${OPTARG}
    ;;
  S) stranded=${OPTARG}
    ;;
  *) echo [`date`] [`hostname`] [$0] Error!! Invalid option >&2
    usage
    exit 1
    ;;
  esac
done

#if [ ${userCustomConfigFile} != "undef" ]
echo [`date`] [`hostname`] [$0] userCustomConfigFile=${userCustomConfigFile}
echo [`date`] [`hostname`] [$0] userCustomSTARConfigFile=${userCustomSTARConfigFile}
if [ ${userCustomConfigFile} != "undef" ] && [ ${userCustomSTARConfigFile} != "undef" ]
then
  echo [`date`] [`hostname`] [$0] Community version of script activated
  consumer="Community"
  elmLoggerStatusStart=5
  elmLoggerStatusComplete=$((elmLoggerStatusStart+1))
  elmLoggerStatusTroubleshooting=7
#  fastqFileDir="userDefined"
#  identifier="userDefined"
#  assembly="userDefined"
  if [ ! -f ${customConfigFile} ] && [ ! -f ${customSTARConfigFile} ]
  then
    echo [`date`] [`hostname`] [$0] configuration file does not exist. Program exiting...
    exit 1
  fi
fi

#if any of the essential options are not provided, flag error and exit with failure
if [ ${elmLoggerStatusStart} == "undef" ] || [ ${fastqFileDir} == "undef" ] || [ ${identifier} == "undef" ] || [ ${consumer} == "undef" ] || [ ${assembly} == "undef" ] || [ ${customConfigFile} == "undef" ] || [ ${customSTARConfigFile} == "undef" ] || [ ${stranded} == "undef" ]
then
  echo [`date`] [`hostname`] [$0] Not all options provided.
  usage
  exit 1
fi
echo [`date`] [`hostname`] [$0] parameters received: elmLoggerStatusStart=${elmLoggerStatusStart} fastqFileDir=${fastqFileDir} customConfigFile=${customConfigFile} customSTARConfigFile=${customSTARConfigFile} identifier=${identifier} consumer=${consumer} assembly=${assembly} comments=${comments} stranded=${stranded}
echo [`date`] [`hostname`] [$0] All parameters=$@
echo [`date`] [`hostname`] [$0] Additional parameters set up to this point: elmLoggerStatusComplete=${elmLoggerStatusComplete} elmLoggerStatusCode=${elmLoggerStatusCode}
echo [`date`] [`hostname`] [$0] Identifier in the folder name -- "${identifier}"

source /home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.config
source ${customConfigFile}

if [ ${elmLoggerEvent} == "Production" ]
  then
        elmlog="EVENTLOG"
  else
        elmlog="EVENTLOG-DEV"
fi
username=`printenv USER`
identifier=`date +%s.%N`
echo [`date`] [`hostname`] [$0] Person submit the job -- "${username}"
echo [`date`] [`hostname`] [$0] Identifier of this library send to elm -- "${identifier}"
libraryFileSizeStr=$(du -sbcD ${fastqFileDir}/*.fastq*|tail -1|cut -f1)
echo [`date`] [`hostname`] [$0] libraryFileSize ${libraryFileSizeStr}

echo [`date`] [`hostname`] [$0] ${customConfigFile}
#To detect the number of parameters passed into this program
####if [ $# -eq 3 ]
####then
####  echo [`date`] [`hostname`] [$0] Default pipeline activated
####  fastqFileDir=${1}
####  assembly=${2}
####  identifier=${3}
####else
####  echo [`date`] [`hostname`] [$0] Invalid number of parameters in command line
####  exit 1
####fi

if [ ${consumer} == "Production" ]
then
  pipelineConsumer="Production"
elif [ ${consumer} == "Community" ]
then
  pipelineConsumer="Community"
else
  echo [`date`] [`hostname`] [$0] 4 parameter is invalid. Program exiting...
  exit 1
fi

#initialise elmLoggerStatusCode of pipeline as 5 "In Analysis".  If there is error, elmLoggerStatusCode is 7 "Troubleshooting"
elmLoggerStatusCode=${elmLoggerStatusStart}
elmCommentsInit=${comments:+$comments:}"STAR RSEM CUFFDIFF worker stage"
elmComments=${elmCommentsInit}:${fastqFileDir}
# new parameters
pipeLineStageHeader="WTS Expression analysis stage"
cpuCount=20

####runDirName=$(echo ${fastqFileDir}|awk -F/ '{print $(NF-4)}')
runDirName=$(echo ${fastqFileDir}|sed -e "s/\/fastqLinks$//" -e "s/^.*\///"|cut -f1 -d".")
runType=$(echo ${runDirName}|cut -d"-" -f2)
runId=$(echo ${runDirName}|cut -f1 -d"_")
assembly=$(echo ${fastqFileDir}|sed -e "s/\/fastqLinks$//" -e "s/^.*\///"|cut -f2 -d".")
libraryId=$(echo ${fastqFileDir}|sed -e "s/\/fastqLinks$//" -e "s/^.*\///"|cut -f3 -d".")
####fastqFile=$(echo ${fastqFileDir}|sed -e "s/.*\///" -e "s/.fastq.gz//")
####readNumber=$(echo ${fastqFile}|cut -d"_" -f4)
####libraryId=$(echo ${fastqFile}|cut -f1 -d"_")
####laneNum=$(echo ${fastqFile}|cut -f3 -d"_"|sed "s/L//"|sed "s/0//g")
echo [`date`] [`hostname`] [$0] runType=${runType} fastqFileDir=${fastqFileDir} fastqFile=${fastqFile} runDirName=${runDirName} assembly=${assembly} readNumber=${readNumber}
laneID=$(ls ${fastqFileDir}|head -1|cut -f3 -d"_"|sed "s/^L0*//")
echo [`date`] [`hostname`] [$0] send to elm LaneID:${laneID}

#The following can be overridden in the user supplied configuration file
#tophat_GtfFile="--GTF ${gtfDir}/${assembly}/gtf/${assembly}_annotation.gtf"
tophat_GtfFilePath="${gtfDir}/${assembly}/gtf/${assembly}_annotation.gtf"
tophat_GenomeBowtieIndex="${genomeDir}/${assembly}/bowtie2_path/base/${assembly}"
tophat_TranscriptomeIndex="--transcriptome-index ${gtfDir}/${assembly}/gtf/${assembly}_annotation"
GenomeFastaFile="${genomeDir}/${assembly}/${assembly}.fa"

rsem_TranscriptomeIndex="${genomeDir}/${assembly}/RSEM_index/RSEMtr_${assembly}"
HTSeq_GtfFile="${gtfDir}/${assembly}/gtf/${assembly}_annotation.gtf"
rnaseqc_GtfFile="${gtfDir}/${assembly}/gtf/${assembly}_RNASeqQCannotation.gtf"

starGenomeDir="--genomeDir ${genomeDir}/${assembly}/star/"
#star_GTFFileParameter="--sjdbGTFfile"
star_GtfFilePath="${gtfDir}/${assembly}/gtf/${assembly}_annotation.gtf"
star_ParameterList="--parametersFiles"
starOutput="--outFileNamePrefix"
starInput="--readFilesIn"
star_ReadGroup="--outSAMattrRGline ID:STAR_2.4.2a   SM:${libraryId}       LB:${libraryId}       PU:LaneId      CN:Genome_Institute_of_Singapore      PI:200  DT:$(date +%Y-%m-%dT%H:%M:%S)  PL:ILLUMINA"

softlinkFileName=`readlink ${star_GtfFilePath}`
echo -e "@CO\tANNOTATIONFILE:${softlinkFileName}" >commentsheader.txt
#softlinkFileName="readlink"
#star_AnnotationSourcePath="$(echo ${softlinkFileName} -f ${star_GtfFilePath}"
#star_AnnotationSourceFile="$(echo ${star_AnnotationSourcePath}|cut -d"/" -f8)"
#echo -e "@CO\tANNOTATIONFILE:${star_AnnotationSourceFile}" >commentsheader.txt
cp -pi commentsheader.txt ${fastqFileDir}
star_Samcomment="--outSAMheaderCommentFile ${fastqFileDir}/commentsheader.txt"

cufflinks_GtfFilePath="${gtfDir}/${assembly}/gtf/${assembly}_annotation.gtf"

#<<COMMENT1
#If customConfigFile is defined and not null value, source user provided parameters
#if [ -z "${customConfigFile:+zzz}" ]
#then
#  echo [`date`] [`hostname`] [$0] Using default parameters for tophat and cufflinks
#  echo [`date`] [`hostname`] source /home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionDefaultParameters.config
#  echo [`date`] [`hostname`] customConfigFile=${customConfigFile}
#  source /home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionDefaultParameters.config
#else
#  echo [`date`] [`hostname`] [$0] Using user supplied StarExpression parameters
#  echo [`date`] [`hostname`] source ${customConfigFile}
#  source ${customConfigFile}
#  assembly=$(echo ${GenomeBowtieIndex}|sed "s/^.*\///")
#  libraryId=${fastqFileDir}
#  identifier=`date +%s.%N`
#fi
#COMMENT1
#######################identifier=`date +%s.%N`

#test if variable is set or not, don't care if the variable is set to null string or set to some value
if [ -z "${debugMode+zzz}" ]
then
  #variable is not set
  DEBUG=0
  echo DEBUG is turned off
else
  #variable is set
  DEBUG=1
  echo DEBUG is turned on
fi

touch ../Jobstatus.txt
echo [`date`] [`hostname`] [$0] ${libraryId} >>../Jobstatus.txt
#Check if assembly file is available
##if [ ! -f ${GenomeFastaFile} ]
if [ ! -f ${GenomeFastaFile} ] || [ ! -f ${tophat_GtfFilePath} ]
then
##  echo [`date`] [`hostname`] [$0] ${assembly}.fa: ${GenomeFastaFile} is not found. Exit program.
echo [`date`] [`hostname`] [$0] ${assembly.fa}: ${GenomeFastaFile} or ${tophat_GtfFile} is not found. Exit program.
echo ${libraryId} Failed >>../Jobstatus.txt
  elmLoggerStatusCode=${elmLoggerStatusTroubleshooting}
  elmComments=${elmCommentsInit}:"No reference found"
  pipeLineStage="WTS Pipeline initializing" # new parameter
  pipeLineStageStatus="Troubleshooting"
####  pipeLineStatus="Troubleshooting" # new parameter
  send_elm_logger_status
  exit 100
fi

#create output directories for tophat and cufflinks
mkdir ./STAR_output
mkdir ./RNASeqQC_output
mkdir ./HTSeq_output
mkdir ./RSEM_output

#Generate a list of fastq.gz files
oldCurrentWorkDir=$(pwd)
cd ${fastqFileDir}
echo Now in Sample directory $(pwd)
libraryFileSizeStr=$(du -sbcD *.fastq*|tail -1|cut -f1)
for eachFastqGzFile in *_R?_*.fastq*
do
  readNumber=$(echo ${eachFastqGzFile}|cut -d"_" -f4)
  if [ ${readNumber} == "R1" ]
  then
    fastqGzR1List=${fastqGzR1List:+$fastqGzR1List,}${eachFastqGzFile}
  elif [ ${readNumber} == "R2" ]
  then
    fastqGzR2List=${fastqGzR2List:+$fastqGzR2List,}${eachFastqGzFile}
  else
    echo Invalid read number. Program exiting...
    elmLoggerStatusCode=${elmLoggerStatusTroubleshooting}
    elmComments=${elmCommentsInit}:"Invalid read number"
  pipeLineStage="WTS pipeline initializing" # new parameter
  pipeLineStageStatus="Troubleshooting"
####    pipeLineStatus="Troubleshooting" # new parameter
    send_elm_logger_status
echo ${libraryId} Failed >>../../Jobstatus.txt
    exit 1
  fi
done
echo [`date`] [`hostname`] [$0] fastqGzR1List=${fastqGzR1List}
echo [`date`] [`hostname`] [$0] fastqGzR2List=${fastqGzR2List}

#Get the lane numnber of the first fastq.gz file
laneNum=$(echo ${fastqGzR1List}|cut -f1 -d","|cut -f3 -d"_"|sed "s/^L0*//")
echo [`date`] [`hostname`] [$0] lane number=${laneNum}

if [ ! -z ${cufflinks_MaskFile} ]
then
cufflinks_MaskFilePath="${cufflinks_MaskFile} ${gtfDir}/${assembly}/gtf/${assembly}_annotation_maskfile.gtf"
else
cufflinks_MaskFilePath=""
  echo [`date`] [`hostname`] [$0] ${cufflinks_MaskFilePath}
fi

if [ ! -z ${cufflinks_FragmentBiasCorrection} ]
then
cufflinks_FragmentBiasCorrectionPath="${cufflinks_FragmentBiasCorrection} ${genomeDir}/${assembly}/bowtie2_path/base/${assembly}.fa"
else
cufflinks_FragmentBiasCorrectionPath=""
  echo [`date`] [`hostname`] [$0] ${cufflinks_FragmentBiasCorrectionPath}
fi

#cd ${oldCurrentWorkDir}
echo [`date`] [`hostname`] [$0] Star started

commandLine="${STARCommand} ${starGenomeDir} ${star_ReadGroup} ${star_Samcomment} ${star_ParameterList} ${customSTARConfigFile} ${starOutput} ${oldCurrentWorkDir}/STAR_output/${libraryId}_${assembly}_ ${starInput} ${fastqGzR1List} ${fastqGzR2List}"
echo [`date`] [`hostname`] [$0] ${commandLine}
pipeLineStage=${pipeLineStageHeader}:" Running Star"
pipeLineStageStatus="Started"
elmComments=${elmCommentsInit}:"running Star"
toolVersion=${starToolVersion}
send_elm_logger_status
             
	${STARCommand} ${starGenomeDir} ${star_ReadGroup} ${star_Samcomment} ${star_ParameterList} ${customSTARConfigFile} ${starOutput} ${oldCurrentWorkDir}/STAR_output/${libraryId}_${assembly}_ ${starInput} ${fastqGzR1List} ${fastqGzR2List}
echo [`date`] [`hostname`] [$0] star completed with exit code=$?
pipeLineStageStatus="Completed"
send_elm_logger_status

StarOutputFileName="${oldCurrentWorkDir}/STAR_output/${libraryId}_${assembly}_Aligned.sortedByCoord.out.bam"
StarTranscriptomeOutputFileName="${oldCurrentWorkDir}/STAR_output/${libraryId}_${assembly}_Aligned.toTranscriptome.out.bam"

####STAR Mapping Pass Or Failed
#if [ ! -f ${oldCurrentWorkDir}/${libraryId}_${assembly}Aligned.sortedByCoord.out.bam ]; then
if [ ! -s ${StarOutputFileName} ] || [ ! -s ${StarTranscriptomeOutputFileName} ]; then ##file exists and has a size greater then zero
##Mail Alert
##echo [`date`] [`hostname`] [$0] RNAseq STAR FAILED for ${libraryId} ${oldCurrentWorkDir} |mail -s "The Star mapping of ${libraryId} Failed under [`pwd`] " chancs@gis.a-star.edu.sg,gnanakkan@gis.a-star.edu.sg
${pipelineLoc}/sendBclReport.pl "echo [`date`] [`hostname`] [$0] Failed: WTS Expression Analysis of ${libraryId} in ${oldCurrentWorkDir} at STAR mapping" "Failed: WTS Expressiom Analysis for ${libraryId} at STAR mapping"
	echo "Star output bam not found!!!"
	echo "Resubmit the pipeline"
#	rm -f STAR_output
echo ${libraryId} Failed >>../../Jobstatus.txt
        exit 1
fi

cd ${oldCurrentWorkDir}/STAR_output

#*****Creating bam index
echo [`date`] [`hostname`] [$0] creating Bam index
commandLine="${samToolsCommand} index ${StarOutputFileName}"
echo [`date`] [`hostname`] [$0] ${commandLine}
   ${samToolsCommand} index ${StarOutputFileName}

ln -s ${StarOutputFileName} ${oldCurrentWorkDir}/STAR_output/${libraryId}_${assembly}_SNV.final.bam
ln -s ${StarOutputFileName}.bai ${oldCurrentWorkDir}/STAR_output/${libraryId}_${assembly}_SNV.final.bam.bai

cd ${oldCurrentWorkDir}

#*****Creating fastqc
echo [`date`] [`hostname`] [$0] creating fastqc Report
commandLine="${FastQCCommand} --format bam_mapped --threads 16 ${StarOutputFileName}"
echo [`date`] [`hostname`] [$0] commandLine="${commandLine}"
pipeLineStage=${pipeLineStageHeader}:" Running FASTQC"
pipeLineStageStatus="Started"
elmComments=${elmCommentsInit}:"Tophat completed, running FASTQC"
toolVersion="${fastqcToolVersion}"
cpuCount=16
send_elm_logger_status

      ${FastQCCommand} --format bam_mapped --threads 16 --outdir ${oldCurrentWorkDir}/  ${StarOutputFileName}
echo [`date`] [`hostname`] [$0] FASTQC completed with exit code=$?
elmComments=${elmCommentsInit}:"FASTQC completed"
pipeLineStageStatus="Completed"
send_elm_logger_status
#echo [`date`] [`hostname`] [$0] fastqc command completed with exit code=$?
#exit 0


#*****RNASeqQC Report
echo [`date`] [`hostname`] [$0] creating RNASeqQC Report
commandLine="java -Xmx3g -XX:ParallelGCThreads=18 -jar ${rnaseqc_dir} -n 1000 -s ${libraryId}|${StarOutputFileName}|RNASeqQC -t ${rnaseqc_GtfFile} -r ${GenomeFastaFile}  -noDoC -o ${oldCurrentWorkDir}/RNASeqQC_output/rnaseqcReport/ "
echo [`date`] [`hostname`] [$0] ${commandLine}
pipeLineStage=${pipeLineStageHeader}:" Running RNASeqQC"
pipeLineStageStatus="Started"
elmComments=${elmCommentsInit}:"FASTQC completed, running RNASeqQC"
toolVersion="${RNASeQCToolVersion}"
cpuCount=18
send_elm_logger_status
cmd_rnaseqc="java -Xmx3g -XX:ParallelGCThreads=18 -jar ${rnaseqc_dir} \
		-n 1000 -s \"${libraryId}|${StarOutputFileName}|RNASeqQC\" -t \
		${rnaseqc_GtfFile} -r ${GenomeFastaFile}  -noDoC -o ${oldCurrentWorkDir}/RNASeqQC_output/rnaseqcReport/ "
     ${cmd_rnaseqc}
echo [`date`] [`hostname`] [$0] RNASeqQC completed with exit code=$?
elmComments=${elmCommentsInit}:"RNASeqQC completed"
pipeLineStageStatus="Completed"
send_elm_logger_status


#*******Creating BAM Statistics
${samToolsCommand} stats ${StarOutputFileName} >${libraryId}_samtools_stats.txt
echo [`date`] [`hostname`] [$0] bam statistics completed with exit code=$?

grep ^SN ${libraryId}_samtools_stats.txt | cut -f 2- >${libraryId}_samtools_SummaryNumbers.txt
${plotBamStatsCommand} -p Samtools_stats/ ${libraryId}_samtools_stats.txt
echo [`date`] [`hostname`] [$0] bam plot completed with exit code=$?
mv -i *samtools* Samtools_stats/

#**********Raw Read Count using HTSeq-count
echo [`date`] [`hostname`] [$0] HTSeq-count started
commandLine="${HTSeqCommand} ${HTSeqParametersList} ${StarOutputFileName} ${HTSeq_GtfFile} _GREATER_SIGN_ HTSeq_output/${libraryId}_${assembly}_htseqdefault.counts.txt"
echo [`date`] [`hostname`] [$0] ${commandLine}
pipeLineStage=${pipeLineStageHeader}:" Running HTSeq-count"
pipeLineStageStatus="Started"
elmComments=${elmCommentsInit}:"RNASeqQC completed, running HTSeq-count"
toolVersion="${HTSeqcountToolVersion}"
cpuCount=1
send_elm_logger_status

             ${HTSeqCommand} ${HTSeqParametersList} ${StarOutputFileName} ${HTSeq_GtfFile} > HTSeq_output/${libraryId}_${assembly}_htseqdefault.counts.txt
echo [`date`] [`hostname`] [$0] HTSeq-count completed with exit code=$?
elmComments=${elmCommentsInit}:"HTSeq-count completed"
pipeLineStageStatus="Completed"
send_elm_logger_status

HTSeqOutput="HTSeq_output/${libraryId}_${assembly}_htseqdefault.counts.txt"

####HTSeq Pass Or Failed
if [ ! -s ${HTSeqOutput} ]; then ##file exists and has a size greater then zero
##Mail Alert
${pipelineLoc}/sendBclReport.pl "echo [`date`] [`hostname`] [$0] Failed: WTS Expression of ${libraryId} in ${oldCurrentWorkDir} at HTSeq count" "Failed: WTS Expression for ${libraryId} at HTSeq count"
        echo "HTSeq output not found!!!"
        echo "Resubmit the pipeline"
echo ${libraryId} Failed HTSeq count >>../Jobstatus.txt
######  exit 1
fi

echo [`date`] [`hostname`] [$0] RSEM Analysis started
if [ -z ${fastqGzR2List} ]; then

commandLine="${RSEMcalculateexpressionCommand} ${rsemParametersListSE} ${StarTranscriptomeOutputFileName} ${rsem_TranscriptomeIndex} RSEM_output/RSEM_${libraryId}_${assembly}"
echo [`date`] [`hostname`] [$0] ${commandLine}
pipeLineStage=${pipeLineStageHeader}:" Running RSEM Calculate Expression"
pipeLineStageStatus="Started"
elmComments=${elmCommentsInit}:"HTSeq-count completed, running RSEM Calculate Expression"
toolVersion="${RSEMToolVersion}"
cpuCount=1
send_elm_logger_status
    ${RSEMcalculateexpressionCommand} ${rsemParametersListSE} ${StarTranscriptomeOutputFileName} ${rsem_TranscriptomeIndex} RSEM_output/RSEM_${libraryId}_${assembly} >RSEM_output/Log.rsem
    ${RSEMbam2wig} RSEM_output/RSEM_${libraryId}_${assembly}.genome.sorted.bam RSEM_output/RSEM_${libraryId}_${assembly}.sorted.wig RSEM_output/RSEM_${libraryId}_${assembly} #To use this command the previous cmd should generate the bam files for that use --output-genome
    ${RSEMPlot} RSEM_output/RSEM_${libraryId}_${assembly} RSEM_output/RSEM_${libraryId}_${assembly}.pdf
echo [`date`] [`hostname`] [$0] RSEM Calculate Expression completed with exit code=$?
elmComments=${elmCommentsInit}:"RSEM Calculate Expression completed"
pipeLineStageStatus="Completed"
send_elm_logger_status

else

commandLine="${RSEMcalculateexpressionCommand} ${rsemParametersListPE} ${StarTranscriptomeOutputFileName} ${rsem_TranscriptomeIndex} RSEM_output/RSEM_${libraryId}_${assembly}"
echo [`date`] [`hostname`] [$0] ${commandLine}
pipeLineStage=${pipeLineStageHeader}:" Running RSEM Calculate Expression"
pipeLineStageStatus="Started"
elmComments=${elmCommentsInit}:"HTSeq-count completed, running RSEM Calculate Expression"
toolVersion="${RSEMToolVersion}"
cpuCount=1
send_elm_logger_status

    ${RSEMcalculateexpressionCommand} ${rsemParametersListPE} ${StarTranscriptomeOutputFileName} ${rsem_TranscriptomeIndex} RSEM_output/RSEM_${libraryId}_${assembly} > RSEM_output/Log.rsem
    ${RSEMbam2wig} RSEM_output/RSEM_${libraryId}_${assembly}.genome.sorted.bam RSEM_output/RSEM_${libraryId}_${assembly}.sorted.wig RSEM_output/RSEM_${libraryId}_${assembly} #To use this command the previous cmd should generate the bam files to this use --output-genome
    ${RSEMPlot} RSEM_output/RSEM_${libraryId}_${assembly} RSEM_output/RSEM_${libraryId}_${assembly}.pdf
echo [`date`] [`hostname`] [$0] RSEM Calculate Expression completed with exit code=$?
elmComments=${elmCommentsInit}:"RSEM Calculate Expression completed"
pipeLineStageStatus="Completed"
send_elm_logger_status
fi

RSEMOutput="RSEM_output/RSEM_${libraryId}_${assembly}.genes.results"
####RSEM Pass Or Failed
if [ ! -s ${RSEMOutput} ]; then ##file exists and has a size greater then zero
##Mail Alert
${pipelineLoc}/sendBclReport.pl "echo [`date`] [`hostname`] [$0] Failed: WTS Expression Analysis for ${libraryId} in ${oldCurrentWorkDir} failed at RSEM" "Failed: WTS Expression Analysis for ${libraryId} at RSEM"
        echo "RSEM output not found!!!"
        echo "Resubmit the pipeline"
echo ${libraryId} Failed RSEM >>../Jobstatus.txt
        exit 1
fi

echo [`date`] [`hostname`] [$0] Cuffdiff started
commandLine="${cuffdiffCommand} ${cufflinksParametersList} ${cufflinks_GtfFilePath} ${cufflinks_MaskFilePath} ${cufflinks_FragmentBiasCorrectionPath} ${StarOutputFileName} ${StarOutputFileName} --output-dir ./Cuffdiff_output"
echo [`date`] [`hostname`] [$0] ${commandLine}

pipeLineStage=${pipeLineStageHeader}:" Running Cuffdiff"
pipeLineStageStatus="Started"
elmComments=${elmCommentsInit}:"star completed, running Cuffdiff"
toolVersion=${cuffdiffToolVersion}
cpuCount=1
send_elm_logger_status

if [ ${stranded} == "-S NonStrandspecific" ]; then
     echo [`date`] [`hostname`] [$0] ${cuffdiffCommand} ${cufflinksParametersList} ${cufflinks_GtfFilePath} ${cufflinks_MaskFilePath} ${cufflinks_FragmentBiasCorrectionPath} ${StarOutputFileName} ${StarOutputFileName} --output-dir ./Cuffdiff_output
     ${cuffdiffCommand} ${cufflinksParametersList} ${cufflinks_GtfFilePath} ${cufflinks_MaskFilePath} ${cufflinks_FragmentBiasCorrectionPath} ${StarOutputFileName} ${StarOutputFileName} --output-dir ./Cuffdiff_output
else
     echo [`date`] [`hostname`] [$0] ${cuffdiffCommand} ${cufflinksParametersListStranded} ${cufflinks_GtfFilePath} ${cufflinks_MaskFilePath} ${cufflinks_FragmentBiasCorrectionPath} ${StarOutputFileName} ${StarOutputFileName} --output-dir ./Cuffdiff_output
     ${cuffdiffCommand} ${cufflinksParametersListStranded} ${cufflinks_GtfFilePath} ${cufflinks_MaskFilePath} ${cufflinks_FragmentBiasCorrectionPath} ${StarOutputFileName} ${StarOutputFileName} --output-dir ./Cuffdiff_output
fi

#      ${cuffdiffCommand} ${cufflinksParametersList} ${cufflinks_GtfFilePath} --mask-file ${gtfDir}/${assembly}/gtf/${assembly}_annotation_maskfile.gtf --frag-bias-correct ${genomeDir}/${assembly}/bowtie2_path/base/${assembly}.fa ${StarOutputFileName} ${StarOutputFileName} --output-dir ./Cuffdiff_output
echo [`date`] [`hostname`] [$0] Cuffdiff completed with exit code=$?
elmComments=${elmCommentsInit}:"Cuffdiff completed"

####Cuffdiff Pass OR Fail Mail alert
if [ ! -s ${oldCurrentWorkDir}/Cuffdiff_output/gene_exp.diff ]; then
${pipelineLoc}/sendBclReport.pl "echo [`date`] [`hostname`] [$0] Failed: WTS Expression Analysis for ${libraryId} in ${oldCurrentWorkDir} failed at cuffdiff" "Failed: WTS Expression Analysis for ${libraryId} at Cuffdiff Analysis"
#echo [`date`] [`hostname`] [$0] STAR RNAseq Cuffdiff failed for ${libraryId} |mail -s "The Cuffdiff analysis Failed ${libraryId} under [`pwd`] " chancs@gis.a-star.edu.sg,gnanakkan@gis.a-star.edu.sg
       echo "Cuffdiff output zero size!"
       echo "resubmit pipeline"
echo ${libraryId} Failed Cuffdiff >>../Jobstatus.txt
       exit 1
else
echo ${libraryId} Passed >>../Jobstatus.txt
fi

cut -f1,5,7,10 Cuffdiff_output/genes.fpkm_tracking |grep -v tracking_id |sort >Cuffdiff_output/genes.fpkm_tracking.col15710

##grep -E -w -v "tracking_id|q2" Cuffdiff_output/genes.read_group_tracking |cut -f1,4 |sort | cut -f1-4,6 | paste - Cuffdiff_output/genes.fpkm_tracking.col15710 | awk '{print $1"\t"$4"\t"$5"\t"$2"\t"$6}' >Cuffdiff_output/${libraryId}_genes_FPKM_Rawreadcount_GIS.txt
grep -E -w -v "tracking_id|q2" Cuffdiff_output/genes.read_group_tracking |cut -f1,4 |sort | cut -f1-4,6 | paste - Cuffdiff_output/genes.fpkm_tracking.col15710 | awk 'BEGIN { FS = "\t" } ;{print $1"\t"$4"\t"$5"\t"$2"\t"$6}' >Cuffdiff_output/${libraryId}_genes_FPKM_Rawreadcount_GIS.txt
sed -i '1s/^/Gene_Id\tGene_Name\tLocus\tReadCount\tFPKM\n/' Cuffdiff_output/${libraryId}_genes_FPKM_Rawreadcount_GIS.txt

#####################Create Report HTML
####commandLine="NA"
elmLoggerStatusCode=${elmLoggerStatusComplete}
pipeLineStageStatus="Completed"
send_elm_logger_status
echo [`date`] [`hostname`] [$0] StarExpressionworkerCommunity.sh completed
exit 0
