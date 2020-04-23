#!/bin/bash

# base name of nra2
# TODO pass from cmdline and print usage otherwise!
nra2="test"

# collect all graphs:
modes=$(ls ${nra2}*_seq_*pfm | sed -e "s/${nra2}\(.*\)_seq_.*/\1/" | sort | uniq)

ref="${nra2}ref_fb00.pfm"
base=$(dirname $0)
diff=${base}/pfmdiff
for mode in $modes
do
> ${mode}.txt
for i in ${nra2}${mode}_seq_*.pfm
do
  echo $i
  rmse=$(${diff} ${ref} ${i} /dev/null | cut -d: -f2)
  usertime=$(cat ${i%.pfm}.txt | grep elapsed|awk '{print $10}'|tr -d 's')
  spp=$(cat ${i%.pfm}.txt | grep "samples per pixel"|awk '{print $6}')
  echo $spp $usertime $rmse >> ${mode}.txt
done
done

# plot equal sample/time plots:
for i in 1 2
do
file=$(mktemp ./conv.XXXXXX)
cat > $file << EOF
set term pdf
set output "conv${i}.pdf"
set ylabel "log rmse"
set xlabel "$([ $i == 1 ] && echo "samples")$([ $i == 2 ] && echo "time [s]")"
set size square
set logscale xy
plot \\
EOF

for mode in $modes
do
  echo "'"${mode}.txt"'" u ${i}:3 title "'"${mode}"'" w l lw 3 ,\\ >> $file
done
  echo '' >> $file
gnuplot $file
rm $file
done
