#!/bin/bash
# gets 2 integer as inputs for a random seed, use batch cluster and proc id
seed=`expr $1 \* \( $2 + 1 \) % 10000`
seed=`expr \( $1 \* 10 + $2 \) % 1000000`
source /cvmfs/sft.cern.ch/lcg/views/LCG_92/x86_64-slc6-gcc62-opt/setup.sh
export PYTHIA8=/afs/cern.ch/user/j/jmontejo/work/LLPtriggertiming/Pythia
export PYTHIA8DATA=/afs/cern.ch/work/j/jmontejo/LLPtriggertiming/Pythia/share/Pythia8/xmldoc
folder="/afs/cern.ch/user/j/jmontejo/work/LLPtriggertiming/Delphes/"
mypufile="generatePU${seed}.cmnd"
sed s/"Random:setSeed = 10"/"Random:setSeed = $seed"/g $folder/examples/Pythia8/generatePileUp.cmnd > $mypufile
sed -i 's/numberOfEvents = 100 /numberOfEvents = 6000 /g' $mypufile #6M events is the limit for one pileup file due to size
cat $mypufile

$folder/DelphesPythia8 $folder/cards/converter_card.tcl $mypufile Minbias.${1}.${2}.root
