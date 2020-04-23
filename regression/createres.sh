#!/bin/bash

function build_config
{
  cp ../config.mk ../config.mk.regressionbk
  cp $1/$2 ../config.mk
  # be sure we use the appropriate path lengths etc. nuke it all:
  make -C .. clean
  make -j16 -C .. corona modules $(cat $1/extratargets 2>/dev/null) >& $1/buildlog
  cp ../config.mk.regressionbk ../config.mk
}

function collect_results
{
  returncode=$2
  # compare to reference
  maxerror="0.11"
  if test -f $1/maxerror
  then
    maxerror=$(cat $1/maxerror)
  fi
  rmse=$(../tools/img/pfmdiff $1/testrender_fb00.pfm $1/reference.pfm /dev/null | cut -d':' -f2)
  pass=pass_$(echo "$(echo $rmse| sed -e 's/[eE]/*10^/') < $maxerror" | bc -l)
  if [ $returncode != 0 ]
  then
    pass=pass_crash
  fi

  txt=$(cat $1/testrender_fb00.txt)

  # create timing graph, if applicable:
  if cat $1/testrender_fb00.txt | grep elapsed >& /dev/null
  then
    wallclock=$(cat $1/testrender_fb00.txt | grep elapsed|awk '{print $4}')
    usertime=$(cat $1/testrender_fb00.txt | grep elapsed|awk '{print $10}')
    # put time in db along with git hash and git log timestamp, make unique and sort
    echo "$(git log -n 1 --pretty=format:"%h %ad" --date=raw) $wallclock" >> $1/wallclock
    echo "$(git log -n 1 --pretty=format:"%h %ad" --date=raw) $usertime " >> $1/usertime
    # sort by field 1, delim ' '
    # first sort by time, so the second unique pass will eliminate the earlier ones, only keep the lastest.
    sort -r -t" " -k2,2 $1/wallclock > $1/tmp
    # make unique by git hash
    sort -u -t" " -k1,1 $1/tmp > $1/tmp2
    # sort by time again
    sort -t" " -k2,2 $1/tmp2 > $1/wallclock
    rm -f $1/tmp $1/tmp2
    sort -r -t" " -k2,2 $1/usertime > $1/tmp
    sort -u -t" " -k1,1 $1/tmp > $1/tmp2
    sort -t" " -k2,2 $1/tmp2 > $1/usertime
    rm -f $1/tmp $1/tmp2

  gnuplot << EOF
set term svg
set output '$1/timings.svg'
set ylabel "time [s]"
set style fill solid 1.0 border -1
set style data histograms
set size square
set key outside
set xtics rotate by -90
plot '$1/wallclock' u 4:xticlabels(1) title 'wallclock' lt rgb "#444444", \
     '$1/usertime'  u 4:xticlabels(1) title 'usertime' lt rgb "#777777"
EOF
  fi

  # output result
  header=$1
  if test -f $1/title
  then
    header=$(cat $1/title)
  fi
  issue=$1
  if test -f $1/issue
  then
    issue=$(cat $1/issue)
  fi

  cat > $1.res << EOF
<h1 class="$pass" onclick="toggle('$1');">
$header
</h1>
<div id="$1" class="content">
<p>$issue</p>
<div><a class="dia" rel="lightbox[viewer]" title="render" href="$1/testrender.png"><span></span><img src="$1/testrender.png" alt="render" class="img"/></a><h3>render</h3>&nbsp;</div>
<div><a class="dia" rel="lightbox[viewer]" title="reference" href="$1/reference.png"><span></span><img src="$1/reference.png" alt="reference" class="img"/></a><h3>reference</h3>&nbsp;</div>
<br/>
<h3>rmse=$rmse, max allowed=$maxerror</h3>
<h3><a href="$1/log">log file</a> -- <a href="$1/buildlog">build log file</a></h3>
<pre>
$txt
</pre>
EOF
  if cat $1/testrender_fb00.txt | grep elapsed >& /dev/null
  then
    echo "<img src='$1/timings.svg' alt='timings' />" >> $1.res
  fi
  echo "</div>" >> $1.res
}

if [ "$2" == "-u" ]
then
  # only re-reate .res file from what we find in the directory,
  # don't actually re-run the test.
  collect_results $1 0
  exit
fi

# first check for completely custom script
if test -x $1/main.sh
then
  cd $1
  test=$1 ./main.sh
  exit
fi

# if there is none, run the usual render invocation.

# run test. customize arguments to corona and
# allow linking the scene file instead of manually copying it.
args="-x -s 16 -w 512 -h 288 -b 0"
referencepass="0"
if test -f $1/args
then
  args=$(cat $1/args)
fi
if test -f $1/scene
then
  scene=$(cat $1/scene)
  cp ${scene}.nra2 $1/test.nra2
  cp ${scene}01.cam $1/test01.cam
  # if the scene is the same, so is the reference image:
  cp ${scene%/*}/reference* $1
fi
if test ! -f $1/reference.png || test ! -f $1/reference.pfm
then
  echo "[$1] could not find reference image, re-creating it!"
  referencepass="1"
  # get custom reference build and args:
  args="$args --batch 128 -s 1024"
  if test -f $1/ref_args
  then
    args=$(cat $1/ref_args)
    echo "[$1] using reference args $args"
  fi
  if test -f $1/ref_config.mk
  then
    echo "[$1] using reference config file $1/ref_config.mk"
    build_config $1 ref_config.mk
  else
    build_config $1 config.mk
  fi
else
  # build regular non-reference version
  build_config $1 config.mk
fi
echo "../corona $1/test.nra2 $args" > $1/log
../corona $1/test.nra2 $args >> $1/log 2>&1 
returncode=$?

# convert to srgb if pfm is in xyz
if [ $returncode == 0 ]
then
if grep 'COL_camera=xyz' $1/config.mk >& /dev/null
then
  convert -flip $1/testrender_fb00.pfm -alpha Off -gamma 1.0 -set colorspace XYZ -colorspace sRGB -depth 8  $1/testrender.png
else
  convert -flip $1/testrender_fb00.pfm -alpha Off -gamma 1.0 -set colorspace RGB -colorspace sRGB -depth 8  $1/testrender.png
fi
# recreate reference if missing
if test ! -f $1/reference.png || test ! -f $1/reference.pfm
then
  cp $1/testrender.png $1/reference.png
  cp $1/testrender_fb00.pfm $1/reference.pfm
  cp $1/testrender_fb00.txt $1/reference.txt
fi
else
  cp style/crash.png $1/testrender.png
fi

if [ $referencepass == "0" ]
then
  collect_results $1 $returncode
else
  cat > $1.res << EOF
<h1 class="$fail" onclick="toggle('$1');">
$1
</h1>
<div id="$1" class="content">
<p>pending re-render after new reference has been created.</p>
EOF
fi

