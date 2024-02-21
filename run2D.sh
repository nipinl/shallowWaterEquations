
postProcess () {

gridSize2d=`tail -n 1 *.r|cut -f 1`
gridSize=`tail -n 1 *.r|cut -f 3`
tSteps=`tail -n 1 *.r|cut -f 4`
gridSizep=$gridSize"p"
gridSized=$gridSize"d"
gridSize2dp=$gridSize2d"p"
gridSize2dd=$gridSize2d"d"
echo "tSteps : " $tSteps
echo "gridSize : " $gridSize
echo "gridSize : " $gridSize2d
file=`ls *.ux`
cp $file $file.bkp
for i in $(seq 1 $tSteps)
do sed -n 1,$gridSizep <$file > file_$i
sed -i 1,$gridSized $file
t=`tail -n 1 file_$i|cut -d "," -f 1`
echo midY-$t
cat<<EOF|gnuplot
#set term postscript eps lw 3 solid 20 colour
set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 600, 400
set datafile separator ","
set size 1.0,0.5
set output 'h_midY$i.png'
set xlabel "x axis" offset 0
set ylabel "height" offset 0
set yrange [0.5:2.5]
set pointintervalbox 0.2
set grid
plot "file_$i" u 2:6 t 'U at MidYplane, t=$t' with linespoints
EOF
done
file=`ls *.uy`
cp $file $file.bkp
for i in $(seq 1 $tSteps)
do sed -n 1,$gridSizep <$file > file_$i
sed -i 1,$gridSized $file
t=`tail -n 1 file_$i|cut -d "," -f 1`
echo midX-$t
cat<<EOF|gnuplot
#set term postscript eps lw 3 solid 20 colour
set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 600, 400
set datafile separator ","
set size 1.0,0.5
set output 'h_midX$i.png'
set xlabel "y axis" offset 0
set ylabel "height" offset 0
set yrange [0.5:2.5]
set pointintervalbox 0.2
set grid
plot "file_$i" u 2:6 t 'U at MidXplane, t=$t' with linespoints
EOF
done
#2D surface plot

file=`ls *.u2`
cp $file $file.bkp
for i in $(seq 1 $tSteps)
do sed -n 1,$gridSize2dp <$file > file_$i
sed -i 1,$gridSize2dd $file
t=`tail -n 1 file_$i|cut -d "," -f 1`
echo 2D-$t
cat<<EOF|gnuplot
set terminal pngcairo enhanced font "arial,10" fontscale 1.0 size 600, 600
set datafile separator ","
set output '2dh$i.png'
set hidden3d
set dgrid3d 50,50 qnorm 2
set view 30, 30, 1, 1
set xyplane at 0
set zrange [0:3]
set zlabel "height" offset 0
set xlabel "x axis" offset 0
set ylabel "y axis" offset 0
splot "file_$i" u 2:3:4 t 'U at t=$t' with lines
EOF

done

rm *.avi
ffmpeg -r 12 -i h_midY%d.png h_midY$file.avi
ffmpeg -r 12 -i h_midX%d.png h_midX$file.avi
ffmpeg -r 12 -i 2dh%d.png 2h$file.avi

mkdir picsNfiles
mv *png file_* picsNfiles
}

rm *.r *.u[x,y,2]
c++ *.cpp
./a.out
ls *.u2>2dFileNames
sed -i 's/...$//' 2dFileNames
while IFS= read -r line; do
    rm -r $line
    mkdir $line
    cp $line.* $line
    cd $line
    if [ ! -d picsNfiles ]; then
        postProcess;
    fi
    cd ..
done < 2dFileNames

