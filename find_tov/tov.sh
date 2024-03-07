#!/bin/bash

echo "ifile (.dat)  ofile (.root)"

datafile=$1
tm=$2
bins=$3
#outfile=$2
prefix="/data/verney/bedo201511/"

#for f in "$prefix/bedo201511_avec_disque_IN2P3_2015-12-03-04:38:44_RUN150.dat" "$prefix/bedo201511_avec_disque_IN2P3_2015-12-03-05:51:23_RUN151.dat"
#do
#datafile= $f
#echo "Processing $datafile"

filelist=$1


while read -r line
do
name="$line"
datafile="$prefix$name"
echo "name of next file $datafile"
name1="$(cut -d'_' -f6<<<$name)"
runn="$(cut -d'.' -f1<<<$name1)"
outfile="$runn""_tov.root"
echo "outfile $outfile"

#exit
maketov()
{
echo "maketov $particle slot" $nslot

if test -f input.txt
then rm input.txt
fi

echo $outfile >>input.txt
# echo 1.5 150 >>input.txt #A=99
echo $tm $bins >>input.txt #A=103
echo $datafile>>input.txt
echo 2 >>input.txt
echo 3 >>input.txt
echo 4 >>input.txt
echo "4 $nslot" >>input.txt
echo $particle>>input.txt
root -l -q -b tov_BE.C++
}

# datafile=/data/data/bedo201511/tov/bedo201511_avec_disque_IN2P3_2015-11-28-22_17_28_RUN58.dat
# outfile=run_58_test.root

if test -f $outfile 
then rm $outfile 
fi

for i in {0..3}
do
   case $i in
   0) 
      particle="_beta_"
      nslot=$i
      maketov
   ;;
   1) 
      particle="_gamma_"
      nslot=$i
      maketov
   ;;
   2) 
      nslot=$i      
      particle="_n_nim_"
      maketov
   ;;
   3) 
      particle="_n_num_"
      nslot=$i
      maketov
   esac
done
done < "$filelist"
