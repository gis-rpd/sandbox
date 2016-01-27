#!/bin/sh
#Created by: Chan Chee Seng
#Version: 2.1
#Revision history: 12 July 2013
#                  21 July 2013 : Pipeline now able to analyse for any sample that comes from multiple projects
#                  27 November 2013: Renamed the ending pipeLineStageStatus from Completed to Ended
#                  4 February 2014: Added condition to check if there are no rnaseq libraries to analyse, which in turn will not create sample output directory
#                  5 June 2014: fork pipeline for Single Cell RNAseq pipeline
#Description: RNAseq pipeline that works for CASAVA1.8.2 directory only.  Uses tophat2 and cufflinks
#Usage: cat "whatsinthisrun.pl output" | RNAseqTophatCufflinksPipeline.sh [path to UnalignedDir directory] [Production]
#-f path to UnalignedDir directory
#-p is whether consumer is Production or not. Optional, when ad-hoc analysis do not use this parameter
#-c is elm logger status code: 5,6,7,10,11,12,13
#-C is user comments for annotating reruns
#-S is stranded information

#Start of function definition
send_elm_logger_status() {
  if [ ${DEBUG} -eq 0 ]
  then
  # actual parameter list for the logger
  # elm_logger_status_rnaseq_pipeline.pl runId laneNum runType libraryId elm_logger_status elmComments pipeLine_Stage pipeLine_Stage_status tool_version command_used unique_identifier cpu_Count library_file_size server_node

    echo [`date`] [`hostname`] [$0] ${pipelineLoggerCommand} --runId ${runId} --laneNum ${laneNum} --runType ${runType} --libraryId ${libraryId} --elmLoggerStatus ${elmLoggerStatusCode} --elmComments "${elmComments}" --pipeLineName "${pipeLineName}" --pipeLineVersion "${pipeLineVersion}" --pipeLineStage "${pipeLineStage}" --pipeLineStageStatus "${pipeLineStageStatus}" --toolVersion "${toolVersion}" --commandLine "${commandLine}" --identifier ${identifier} --cpuCount ${cpuCount} --libraryFileSize ${libraryFileSizeStr} --consumer Production
                                    ${pipelineLoggerCommand} --runId ${runId} --laneNum ${laneNum} --runType ${runType} --libraryId ${libraryId} --elmLoggerStatus ${elmLoggerStatusCode} --elmComments "${elmComments}" --pipeLineName "${pipeLineName}" --pipeLineVersion "${pipeLineVersion}" --pipeLineStage "${pipeLineStage}" --pipeLineStageStatus "${pipeLineStageStatus}" --toolVersion "${toolVersion}" --commandLine "${commandLine}" --identifier ${identifier} --cpuCount ${cpuCount} --libraryFileSize ${libraryFileSizeStr} --consumer Production
  fi
}
#End of function definition

#main
#definition of all options variables
runDirName="undef"
unalignedDirName="undef" #${1}
consumer="undef" #${2}
elmLoggerStatusStart="undef" #${3}
comments="" #${4}
stranded="undef"

if ( ! getopts "c:f:p:u:U:q:C:S:r:" FLAG )
then
  echo [`date`] [`hostname`] [$0] No options provided.
  echo [`date`] [`hostname`] [$0] Usage: RNAseqStarExpressionPipeline.sh -r runDirName -f /path/to/UnalignedDirFastqFiles -p "Production|Community" -u userCustomConfigFile -U userCustomSTARConfigFile -q "All users others than userrig must provide queue name excluding rnaseqpipeline.q" -c elmLoggerStatusCode -S "Strandspecific|NonStrandspecific"
  exit 1
fi

while getopts "c:f:p:u:U:q:C:S:r:" FLAG
do
  case $FLAG in
    c) elmLoggerStatusStart=${OPTARG}
      elmLoggerStatusComplete=$((elmLoggerStatusStart+1))
      elmLoggerStatusTroubleshooting=7
      ;;
    f) unalignedDirName=${OPTARG}
      ;;
    p) consumer=${OPTARG}
      ;;
    r) runDirName=${OPTARG}
      ;;
    u) userCustomConfigFile=${OPTARG}
      ;;
    U) userCustomSTARConfigFile=${OPTARG}
      ;;
    q) userQueueName=${OPTARG}
      ;;
    C) comments=${OPTARG}
      ;;
    S) stranded=${OPTARG}
      ;;
    *) echo [`date`] [`hostname`] [$0] Error!! Invalid option >&2
      exit 1
      ;;
  esac
done
#If customConfigFile is defined and not null value, source user provided parameters
if [ -z "${userCustomConfigFile:+zzz}" ]
then
  echo [`date`] [`hostname`] userCustomConfigFile="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionDefaultParameters.config"
  userCustomConfigFile="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionDefaultParameters.config"
  echo ${userCustomConfigFile}
        if [ -z "${userCustomSTARConfigFile:+zzz}" ]
        then
        echo [`date`] [`hostname`] userCustomSTARConfigFile="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/Parameters1.in"
        userCustomSTARConfigFile="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/Parameters1.in"
        echo [`date`] [`hostname`] ${userCustomSTARConfigFile}
        else
               echo [`date`] [`hostname`] ${userCustomSTARConfigFile}
        fi
else
  userCustomConfigFile=${userCustomConfigFile}
  echo [`date`] [`hostname`] ${userCustomConfigFile}
        if [ -z "${userCustomSTARConfigFile:+zzz}" ]
        then
        echo [`date`] [`hostname`] userCustomSTARConfigFile="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/Parameters1.in"
        userCustomSTARConfigFile="/home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/Parameters1.in"
        echo [`date`] [`hostname`] ${userCustomSTARConfigFile}
        else
                echo [`date`] [`hostname`] ${userCustomSTARConfigFile}
        fi
fi

touch Jobstatus.txt
touch runParameters.txt
echo [`date +%s.%N`] [`hostname`] [$0] "RNAseqStarExpressionPipeline.sh -r ${runDirName} -f ${unalignedDirName} -p ${consumer} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -q ${userQueueName} -c ${elmLoggerStatusCode} -S ${stranded}" >>../runParameters.txt

#if any of the essential options are not provided, flag error and exit with failure
if [ ${elmLoggerStatusStart} == "undef" ] || [ ${unalignedDirName} == "undef" ] || [ ${runDirName} == "undef" ] || [ ${consumer} == "undef" ] || [ ${stranded} == "undef" ]
then
  echo [`date`] [`hostname`] [$0] Not all options provided.
  echo [`date`] [`hostname`] [$0] Usage: RNAseqStarExpressionPipeline.sh -r ${runDirName} -f /path/to/UnalignedDir -p "Production|Community" ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -q ${userQueueName} -c elmLoggerStatusCode -S stranded
  exit 1
fi
echo [`date`] [`hostname`] [$0] parameters received: elmLoggerStatusStart=${elmLoggerStatusStart} unalignedDirName=${unalignedDirName} userCustomConfigFile=${userCustomConfigFile} userCustomSTARConfigFile=${userCustomSTARConfigFile} consumer=${consumer} stranded=${stranded}
echo [`date`] [`hostname`] [$0] All parameters=$@
echo [`date`] [`hostname`] [$0] Additional parameters set up to this point: elmLoggerStatusComplete=${elmLoggerStatusComplete} elmLoggerStatusCode=${elmLoggerStatusCode}

#Generating a unique identifider using timestamp in seconds since 1970-01-01 00:00:00 UTC and nanoseconds
identifier=`date +%s.%N`
echo [`date`] [`hostname`] identifier: ${identifier}

source /home/userrig/pipelines/WholeTranscriptomesequenceAnalysis/RNAseqStarExpressionPipeline.config
#initialise elmLoggerStatusCode of pipeline as 5 "In Analysis".  If there is error, elmLoggerStatusCode is 7 "Troubleshooting"
elmLoggerStatusCode=${elmLoggerStatusStart}
elmCommentsInit=${comments:+$comments:}"Initialising WTSExpression pipeline"
pipeLineStage="WTSExpression pipeline main program"
pipeLineStageStatus="Initiated"
toolVersion="WTSExpression pipeline main program"
cpuCount=1

####consumer=${2}
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

#test if variable is set or not, don't care if the variable is set to null string or set to some value
if [ -z "${debugMode+zzz}" ]
then
  #variable is not set
  DEBUG=0
  DebugSwitch=""
  echo DEBUG is turned off
else
  #variable is set
  DEBUG=1
  DebugSwitch="-v debugMode"
  echo DEBUG is turned on
fi

#Check for correct number of parameters
if [ $# -lt 1 ]
then
  echo [`date`] [`hostname`] [$0] No argument detected.  Program exiting.
  elmLoggerStatusCode=7
  exit 1
else
#sequenceLocation is now up to the Unaligned directory and NOT to the run directory as previously was the case
####  sequenceLocation=$(dirname ${1})
#Trim last trailing forward slash
  sequenceLocation=$(echo ${unalignedDirName}|sed "s/\/$//")
  if [ ${runDirName} == "undef" ]
  then
    runDirName=$(dirname ${sequenceLocation}|awk -F/ '{print $(NF-1)}')
  fi
  echo [`date`] [`hostname`] [$0] sequenceLocation: ${sequenceLocation}
fi

#Check if sequenceLocation contains the fastq output directory
if [ ! -d ${sequenceLocation} ]
then
  echo [`date`] [`hostname`] [$0] ${sequenceLocation} does not exist
  elmLoggerStatusCode=7
  exit 1
else
  echo [`date`] [`hostname`] [$0] ${sequenceLocation} is found here
fi

#Check for absolute path to libraryArchive
if [ -d /mnt/AnalysisPool/libraryArchive ]
then
  libraryArchive="/mnt/AnalysisPool/libraryArchive"
else
  libraryArchive="/mnt/libraryArchive"
fi

#Check if current environemnt supports SGE or LSF
if [ ! -z "${SGE_ROOT+zzz}" ] && [ ${SGE_ROOT} == ${sge_root} ]
then
  echo [`date`] [`hostname`] [$0] UGE detected, using UGE job submit command string
  if [ ! -z ${userQueueName} ]
  then
    queueName="-q ${userQueueName}"
  else
    queueName="-q rnaseqpipeline.q"
  fi
#  jobSubmitString="qsub -pe OpenMP 20 ${queueName} -cwd -v PATH -terse ${DebugSwitch}"
  jobSubmitString="qsub -pe OpenMP 20 -l h_rt=10:00:00,mem_free=40G,hostname=n073 ${queueName} -cwd -v PATH -terse ${DebugSwitch}"

  echo [`date`] [`hostname`] [$0] jobSubmitString=${jobSubmitString}
elif [ ! -z "${LSF_SERVERDIR+zzz}" ] && [ ${LSF_SERVERDIR} == ${lsf_serverdir} ]
then
  echo [`date`] [`hostname`] [$0] LSF detected, using LSF job submit command string
  jobSubmitString="bsub -q long -cwd"
  echo [`date`] [`hostname`] [$0] jobSubmitString=${jobSubmitString}
else
  echo [`date`] [`hostname`] [$0] No scheduler detected!
  jobSubmitString=""
  echo [`date`] [`hostname`] [$0] jobSubmitString=${jobSubmitString}
fi

while read line
do
  libraryId=`echo $line|cut -f3 -d" "`
  projectId=${libraryId}
  runId=`echo $line|cut -f1 -d" "`
  laneNum=`echo $line|cut -f2 -d" "`
  assembly=`echo $line|cut -f4 -d" "`
  machname=`echo ${runId}|cut -f1 -d"-"`
  runType=$(echo ${runId}|cut -d"-" -f2)
  barCode="NoIndex"
  echo [`date`] [`hostname`] [$0] ${projectId} ${libraryId} ${runId} ${laneNum} $machname
  if [ -d `ls -d ${sequenceLocation}` ]
  then
    if [[ ${libraryId} =~ MUX ]]
    then
      assembly=`echo $line|cut -f6 -d" "`
      libraryId=`echo $line|cut -f7 -d" "`
      barCode=`echo $line|cut -f5 -d" "`
      echo [`date`] [`hostname`] [$0] MUX library: ${libraryId} ${barCode}
    fi
    #Check if assembly file is available
    if [ ! -f ${genomeDir}/$assembly/bwa_path/nucleotide/$assembly.fa ]
    then
      echo [`date`] [`hostname`] [$0] $assembly.fa is not found. Exit program.
echo No assembly ${assembly} Failed >>Jobstatus.txt
      exit 1
    fi
    echo [`date`] [`hostname`] [$0] ${projectId} ${libraryId} ${runId} ${laneNum}
    fastqFileDir="${sequenceLocation}/Project_${projectId}/Sample_${libraryId}"
    #Check if sample directory exist here
    if [ -d ${fastqFileDir} ]
    then
      echo Sample directory exists
    else
      echo [`date`] [`hostname`] [$0] ${fastqFileDir} does not exist. Program exiting.
echo fastqFileDir Failed >>Jobstatus.txt
      exit 1
    fi
    if [ -d ${runDirName}.${assembly}.${libraryId}.${identifier}/fastqLinks ]
    then
      cd ${runDirName}.${assembly}.${libraryId}.${identifier}/fastqLinks
      echo [`date`] [`hostname`] [$0] I am currently in this directory: $(pwd)
    else
      mkdir ${runDirName}.${assembly}.${libraryId}.${identifier}
      mkdir ${runDirName}.${assembly}.${libraryId}.${identifier}/fastqLinks
      #Added new library into ${runDirName}.${assembly}.${libraryId} only if library directory does not already exist
      echo [`date`] [`hostname`] [$0] runDirNameAssemblyLibraryId: ${runDirNameAssemblyLibraryId}
      echo [`date`] [`hostname`] [$0] currently processing: ${runDirName}.${assembly}.${libraryId}
      runDirNameAssemblyLibraryId="${runDirNameAssemblyLibraryId:+$runDirNameAssemblyLibraryId }${runDirName}.${assembly}.${libraryId}"
      cd ${runDirName}.${assembly}.${libraryId}.${identifier}/fastqLinks
      echo [`date`] [`hostname`] [$0] Created directories. I am currently in this directory: $(pwd)
    fi
    #Check if fastq files are gzipped
    if [ $(ls -d ${fastqFileDir}/${libraryId}[_-]${barCode}*L00${laneNum}_R*|grep "fastq.gz"|wc -l) -gt 0 ]
    then
      fastqFormat="fastq.gz"
      echo [`date`] [`hostname`] [$0] fastq files are compressed. fastqFormat=${fastqFormat}
    elif [ $(ls -d ${fastqFileDir}/${libraryId}[_-]${barCode}*L00${laneNum}_R*|grep fastq|wc -l) -gt 0 ]
    then
      fastqFormat="fastq"
      echo [`date`] [`hostname`] [$0] fastq files are NOT compressed. fastqFormat=${fastqFormat}
    fi
    echo [`date`] [`hostname`] [$0] ${fastqFileDir} ${runDirName}.${assembly}.${libraryId}
    for eachFastqFile in ${fastqFileDir}/${libraryId}[_-]${barCode}*L00${laneNum}_R*.${fastqFormat}
    do
      ln -s ${eachFastqFile}
    done
    cd ../..
  else
    echo [`date`] [`hostname`] [$0] Library cannot be found!
    elmLoggerStatusCode=7
    pipeLineStageStatus="Troubleshooting"
    elmComments=${elmCommentsInit}:"Library cannot be found"
    send_elm_logger_status
echo No ${libraryId} Failed >>Jobstatus.txt
    exit 1
  fi
done

echo [`date`] [`hostname`] [$0] I am currently in this directory: $(pwd)

if [ $(ls -d $(pwd)/${runDirName}.*|wc -l) -eq 0 ]
then
  echo [`date`] [`hostname`] [$0] No rnaseq sample output directory found here. Program exiting...
echo rnaseq sample output directory Failed >>Jobstatus.txt
  exit 1
fi
#######################################################################
touch Jobstatus.txt
count=`cat config_star_expression.txt |cut -f3|sort|uniq|wc -l`
echo ${count} Number of MUX/NonMUX in this Job >>Jobstatus.txt
count=`cat config_star_expression.txt |cut -f7|sort|uniq|wc -l`
echo ${count} Number of Libraries in this Job >>Jobstatus.txt
#######################################################################

for eachLibraryDir in $(pwd)/${runDirName}.*
do
  cd ${eachLibraryDir}
  libraryId=$(echo ${eachLibraryDir}|sed "s/^.*\///"|cut -f3 -d".")
  assembly=$(echo ${eachLibraryDir}|sed "s/^.*\///"|cut -f2 -d".")
  laneNum=$(ls ${eachLibraryDir}/fastqLinks|head -1|cut -f3 -d"_"|sed "s/^L0*//")
 echo [`date`] [`hostname`] [$0] ${jobSubmitString} ${pipelineLoc}/StarExpressionworker.sh -f ${eachLibraryDir}/fastqLinks -a ${assembly} -i ${identifier} -p ${pipelineConsumer} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -c ${elmLoggerStatusStart} ${comments:+-C \"$comments\"} -S ${stranded} >>../runParameters.txt

 commandLine="${jobSubmitString} ${pipelineLoc}/StarExpressionworker.sh -f ${eachLibraryDir}/fastqLinks -a ${assembly} -i ${identifier} -p ${pipelineConsumer} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -c ${elmLoggerStatusStart} ${comments:+-C $comments} -S ${stranded}"

 JobIdWorker=$(${jobSubmitString} ${pipelineLoc}/StarExpressionworker.sh -f ${eachLibraryDir}/fastqLinks -a ${assembly} -i ${identifier} -p ${pipelineConsumer} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -c ${elmLoggerStatusStart} ${comments:+-C "$comments"} -S ${stranded})
  sleep 10 
  libraryFileSizeStr=$(du -sbcD ${eachLibraryDir}/fastqLinks/*.gz|tail -1|cut -f1)
  elmComments=${elmCommentsInit}:"Submitting RNAseq jobs to scheduler"
  send_elm_logger_status
  cd ..
  JobIdWorkerAll=${JobIdWorkerAll:+$JobIdWorkerAll,}${JobIdWorker}
done

#wait for all StarExpressionworker jobs to complete
echo [`date`] [`hostname`] [$0] commandLine=${jobSubmitString} -m e -M "gnanakkan@gis.a-star.edu.sg","lees@gis.a-star.edu.sg" -hold_jid ${JobIdWorkerAll} ${pipelineLoc}/postEndedStatus.sh -f $(pwd) -a ${assembly} -i ${identifier} -r ${runDirName} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${comments:+-C "$comments"}
#                     JobIdpostEndedStatus=$(${jobSubmitString} -m e -M "gnanakkan@gis.a-star.edu.sg","lees@gis.a-star.edu.sg" -hold_jid ${JobIdWorkerAll} ${pipelineLoc}/postEndedStatus.sh -f $(pwd) -a ${assembly} -i ${identifier} -r ${runDirName} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${comments:+-C "$comments"})

jobSubmitString1="qsub -pe OpenMP 1 -q bwapipeline.q -cwd -v PATH -terse ${DebugSwitch}"
                    JobIdpostEndedStatus=$(${jobSubmitString1} -m e -M "gnanakkan@gis.a-star.edu.sg","lees@gis.a-star.edu.sg" -hold_jid ${JobIdWorkerAll} ${pipelineLoc}/postEndedStatus.sh -f $(pwd) -a ${assembly} -i ${identifier} -r ${runDirName} -p ${pipelineConsumer} -c ${elmLoggerStatusStart} ${comments:+-C "$comments"})
echo [`date`] [`hostname`] [$0] JobIdpostEndedStatus=${JobIdpostEndedStatus}
############################################### Resubmit Failed jobs ###############################################
#for eachLibraryDir in $(pwd)/${runDirName}.*
#do
#  cd ${eachLibraryDir}
#  echo ${eachLibraryDir}
#  libraryId=$(echo ${eachLibraryDir}|sed "s/^.*\///"|cut -f3 -d".")
#  assembly=$(echo ${eachLibraryDir}|sed "s/^.*\///"|cut -f2 -d".")
#  laneNum=$(ls ${eachLibraryDir}/fastqLinks|tail -1|cut -f3 -d"_"|sed "s/^L0*//")
#  echo ${laneNum}
#  failedJobId=$(echo ${eachLibraryDir}/RNAseqSTARExpression.o* |cut -f6 -d"."|sed -e "s/o//")
#  resubmitJobId=$(qacct -j ${failedJobId} |grep "exit_status" |cut -f3 -d" ")
#  OR use JobIdWorker (contain qsub jobid)
#  if [ "${resubmitJobId}" -ne 0 ]
#  then
#       rm -rf RNAseqSTARExpression.o* STAR_output RNASeqQC_output HTSeq_output Cuffdiff_output RSEM_output Samtools_stats
#       echo [`date`] [`hostname`] [$0] ${jobSubmitString} ${pipelineLoc}/StarExpressionworker.sh -f ${eachLibraryDir}/fastqLinks -a ${assembly} -i ${identifier} -p ${pipelineConsumer} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -c ${elmLoggerStatusStart} ${comments:+-C \"$comments\"} >>../runParameters.txt
#       commandLine="${jobSubmitString} ${pipelineLoc}/StarExpressionworker.sh -f ${eachLibraryDir}/fastqLinks -a ${assembly} -i ${identifier} -p ${pipelineConsumer} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -c ${elmLoggerStatusStart} ${comments:+-C $comments}"
#       JobIdWorker=$(${jobSubmitString} ${pipelineLoc}/StarExpressionworker.sh -f ${eachLibraryDir}/fastqLinks -a ${assembly} -i ${identifier} -p ${pipelineConsumer} ${userCustomConfigFile:+-u "$userCustomConfigFile"} ${userCustomSTARConfigFile:+-U "$userCustomSTARConfigFile"} -c ${elmLoggerStatusStart} ${comments:+-C "$comments"})
#       libraryFileSizeStr=$(du -sbcD ${eachLibraryDir}/fastqLinks/*.gz|tail -1|cut -f1)
#       elmComments=${elmCommentsInit}:"Submitting RNAseq jobs to scheduler"
#       send_elm_logger_status
#   fi
#   cd ..
#done

commandLine="NA"
# send logger to say pipeline submission completed
elmComments=${comments:+$comments:}"WTSExpression submission is completed"
pipeLineStageStatus="Completed"
send_elm_logger_status

exit 0
