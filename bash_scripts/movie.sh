#!/bin/bash

if [ $# -lt 3 ]; then
    echo "usage: movie output_movie FPS files"
    kill -INT $$
fi
tmpdir=`mktemp -d`
output=$1
fps=$2
shift 2
files="$@"
count=0
suffix=`echo $1 | cut -f2 -d.`
width=`identify -format "%w" $1`
height=`identify -format "%h" $1`
# Make sure that the width is even (yuv420p format for QT)
width=`expr $width \/ 2 \* 2`
height=`expr $height \/ 2 \* 2`
cwd=`pwd`
for f in $files; do
    fn=`echo $count $suffix | gawk '{printf "img%06d.%s", $1, $2}'`
    ln -s $cwd/$f ${tmpdir}/${fn}
    count=`expr $count + 1`
done
ffmpeg -r $fps -f image2 -i ${tmpdir}/img%06d.${suffix} -vcodec libx264 -pix_fmt yuv420p \
       -s ${width}x${height} -crf 25 -y $output
rm -r $tmpdir
