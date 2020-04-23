#!/bin/bash

allpass=1
N=0
while read output; do
  N=$((N+1))
  ebsdf[$N]=$(echo $output | cut -d" " -f 2)
  bsdf[$N]=$(echo $output | cut -d" " -f 3)
  epdf[$N]=$(echo $output | cut -d" " -f 4)
  pdf[$N]=$(echo $output | cut -d" " -f 5)
  diffpdf[$N]=$(echo ${pdf[$N]} - ${epdf[$N]} | bc -l)
  diffbsdf[$N]=$(echo ${bsdf[$N]} - ${ebsdf[$N]} | bc -l)
  pdfpass[$N]=$(echo "${diffpdf[$N]}*${diffpdf[$N]} < 0.00001" | bc -l)
  bsdfpass[$N]=$(echo "${diffbsdf[$N]}*${diffbsdf[$N]} < 0.00001" | bc -l)
  if [[ ${pdfpass[$N]} == "0" ]]
  then
    allpass=0
  fi
  if [[ ${bsdfpass[$N]} == "0" ]]
  then
    allpass=0
  fi
done < output.log

mogrify -format png *.pgm
for i in $(seq -f%02g 1 1 $count)
do
  composite ${i}_ebsdf.png ${i}_bsdf.png -compose difference ${i}_bsdf_diff.png
  composite ${i}_epdf.png ${i}_pdf.png -compose difference ${i}_pdf_diff.png
  mogrify -auto-level ${i}_bsdf_diff.png
  mogrify -auto-level ${i}_pdf_diff.png
done
cd ..

# only pass if all the subimages pass.
pass="pass_$allpass"
cat > $test.res << EOF
<h1 class="$pass" onclick="toggle('$test');">
$header
</h1>
<div id="$test" class="content">
<p>bsdf battle-test running sample(), estimating the bsdf the pdf using a histogram.
also testing eval() and pdf() calls directly.</p>
EOF

for i in $(seq -f%02g 1 1 $count)
do
  title="angle $i"
  img="$test/$i"
  echo "<div><a class='dia' rel='lightbox[viewer]' title='$title' href='${img}_ebsdf.png'><span></span><img src='${img}_ebsdf.png' alt='$title' class='img'/></a><h3>estimated bsdf $title</h3><h3>sum ${ebsdf[$i]}</h3>&nbsp;</div>" >> ${test}.res
  echo "<div><a class='dia' rel='lightbox[viewer]' title='$title' href='${img}_bsdf.png'><span></span><img src='${img}_bsdf.png' alt='$title' class='img'/></a><h3>bsdf eval $title</h3><h3>sum ${bsdf[$i]}</h3>&nbsp;</div>" >> ${test}.res
  echo "<div class='diffimgpass_${bsdfpass[$i]}'><a class='dia' rel='lightbox[viewer]' title='$title' href='${img}_bsdf_diff.png'><span></span><img src='${img}_bsdf_diff.png' alt='$title' class='img'/></a><h3>bsdf diff $title</h3><h3>diff ${diffbsdf[$i]}</h3>&nbsp;</div>" >> ${test}.res
  echo "<div><a class='dia' rel='lightbox[viewer]' title='$title' href='${img}_epdf.png'><span></span><img src='${img}_epdf.png' alt='$title' class='img'/></a><h3>estimated pdf $title</h3><h3>sum ${epdf[$i]}</h3>&nbsp;</div>" >> ${test}.res
  echo "<div><a class='dia' rel='lightbox[viewer]' title='$title' href='${img}_pdf.png'><span></span><img src='${img}_pdf.png' alt='$title' class='img'/></a><h3>pdf eval $title</h3><h3>sum ${pdf[$i]}</h3>&nbsp;</div>" >> ${test}.res
  echo "<div class='diffimgpass_${pdfpass[$i]}'><a class='dia' rel='lightbox[viewer]' title='$title' href='${img}_pdf_diff.png'><span></span><img src='${img}_pdf_diff.png' alt='$title' class='img'/></a><h3>pdf diff $title</h3><h3>diff ${diffpdf[$i]}</h3>&nbsp;</div>" >> ${test}.res
done
echo "&nbsp;</div>" >> $test.res
