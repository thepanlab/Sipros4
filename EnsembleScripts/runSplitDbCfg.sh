#!/bin/bash 

# Get path of run script
if [ -L $0 ] ; then
    exePath=$(dirname $(readlink -f $0)) ;
else
    exePath=$(dirname $0) ;
fi;

DbPath=""
CfgPath=""
DbOutputFolder=""
CfgOutputFolder=""
numSplit=""


# Read Parameter
while [[ $# > 0 ]]
do
key="$1"
case $key in
    -h|--help)              # help output
    echo -e "Usage:\n"
    echo -e "   runSplitDbCfg.sh [OPTION]...<PARAM>...\n\n"
    echo -e "<PARAMS>\n"
    echo -e "   -d\t protein database filename.\n"
    echo -e "   -c\t configuration filename.\n"
    echo -e "   -od\t output folder of of split database (If not exist, will create one).\n"
    echo -e "   -oc\t output folder of of split configuration file (If not exist, will create one).\n"
    echo -e "   -n\t number of designed split.\n"
    echo -e "<OPTIONS>\n"
    echo -e "   -h\t help.\n"
    exit 1
    ;;
    -od)			# output directory database
    DbOutputFolder="$2"
    shift # past argument
    ;;
    -oc)			# output directory configuration
    CfgOutputFolder="$2"
    shift # past argument
    ;;
    -d)					# database
    DbPath="$2"
    shift # past argument
    ;;
    -c)					# configuration
    CfgPath="$2"
    shift # past argument
    ;;
    -n)
    numSplit="$2"
    shift
    ;;
    *)
    echo "ERROR: Unidentified user variable $key"
    exit 1        				# unknown option
    ;;
esac
shift # past argument or value
done

if [ -z "$DbPath" ] && [ -z "$CfgPath" ]  && [ -z "$DbOutputFolder" ] && [ -z "$CfgOutputFolder" ] && [ -z "$numSplit" ]; then
   echo "Input incomplete. Not all required parameters specified. Run with -h for help. Exiting..."
   exit 1
fi

if [ ! -d ${DbOutputFolder} ]; then
	mkdir -p ${DbOutputFolder}
fi

if [ ! -d ${CfgOutputFolder} ]; then
	mkdir -p ${CfgOutputFolder}
fi

#Check required binaries are in the bin directory
PYSPLITDB="${exePath}/split_db.py"
if [ -f $PYSPLITDB ] ; then
   echo "$PYSPLITDB exists."
else
   echo "$PYSPLITDB does not exist in script directory."
fi

PYSPLITCFG="${exePath}/configure_subdb.py"
if [ -f $PYSPLITCFG ] ; then
   echo "$PYSPLITCFG exists."
else
   echo "$PYSPLITCFG does not exist in script directory."
fi

# split database
FileName=$(python ${PYSPLITDB} ${DbPath} ${DbOutputFolder} ${numSplit})

# split configuration file
python ${PYSPLITCFG} ${CfgPath} ${CfgOutputFolder} ${numSplit} ${FileName}

echo "Done."