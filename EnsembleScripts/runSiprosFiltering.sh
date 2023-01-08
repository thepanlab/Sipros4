#!/bin/bash 

#-------------------------------------------------------------------#

# Get path of run script
if [ -L $0 ] ; then
    exePath=$(dirname $(readlink -f $0)) ;
else
    exePath=$(dirname $0) ;
fi;

SipFolder=""
ConfigureFile=""
OutputFolder=""
TabFile=""

# Read Parameter
while [[ $# > 0 ]]
do
key="$1"
case $key in
    -h|--help)              # help output
    echo -e "Usage:\n"
    echo -e "   runSiprosPostprocessing.sh [OPTION]...<PARAM>...\n\n"
    echo -e "<PARAMS>\n"
    echo -e "   -in\t SIP file directory (directory containing search results from Sipros Ensemble).\n"
    echo -e "   -c\t configuration file.\n"
    echo -e "   -o\t output folder (If not exist, will create one).\n"
    echo -e "<OPTIONS>\n"
    echo -e "   -h\t help.\n"
    exit 1
    ;;
    -o)			# output directory
    OutputFolder="$2"
    shift # past argument
    ;;
    -in)					# Forward paired end read file -- single file
    SipFolder="$2"
    shift # past argument
    ;;
    -t)
    TabFile="$2"
    shift
    ;;
    -c)			# Output file prefix
    ConfigureFile="$2"
    shift # past argument
    ;;
    *)
    echo "ERROR: Unidentified user variable $key"
    exit 1        				# unknown option
    ;;
esac
shift # past argument or value
done

if [ -z "$SipFolder" ] && [ -z "$ConfigureFile" ]  && [ -z "$OutputFolder" ] ; then
   echo "Input incomplete. Not all required parameters specified. Run with -h for help. Exiting..."
   exit 1
fi

if [ ! -d ${OutputFolder} ]; then
	mkdir -p ${OutputFolder}
fi


# Generate PSM table
if [ "$TabFile" == "" ]; then
TabFile=$(python2 ${exePath}/sipros_psm_tabulating.py -i ${SipFolder}/ -c ${ConfigureFile} -o ${OutputFolder}/)
fi

# PSM Filtering
# echo ${TabFile}
python2 ${exePath}/sipros_ensemble_filtering.py -i ${TabFile} -c ${ConfigureFile} -o ${OutputFolder}/

# Protein Assembly
python2 ${exePath}/sipros_peptides_assembling.py -w ${OutputFolder}/ -c ${ConfigureFile}

# Protein SIP mode clustering
python2 ${exePath}/ClusterSip.py -w ${OutputFolder}/ -c ${ConfigureFile}
