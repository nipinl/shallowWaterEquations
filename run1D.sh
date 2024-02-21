
postProcess () {

gridSize=`tail -n 1 *.r|cut -f 1`
tSteps=`tail -n 1 *.r|cut -f 2`
gridSizep=$gridSize"p"
gridSized=$gridSize"d"
echo "tSteps : " $tSteps
echo "gridSize : " $gridSize
echo "gridSizep : " $gridSizep
echo "gridSized : " $gridSized
file=`ls *.u`
cp $file $file.bkp
for i in $(seq 1 $tSteps)
do sed -n 1,$gridSizep <$file > file_$i
sed -i 1,$gridSized $file
t=`tail -n 1 file_$i|cut -d "," -f 1`
cat<<EOF|gnuplot
set term postscript eps lw 3 solid 20 colour
set datafile separator ","
set size 1.0,0.5
set output 'h$i.eps'
set yrange [0.5:3.5]
set pointintervalbox 0.2
set grid
plot "file_$i" u 2:6 t 'u at t=$t' with linespoints
EOF
cat<<EOF|gnuplot
set term postscript eps lw 3 solid 20 colour
set datafile separator ","
set size 1.0,0.5
set output 'hu$i.eps'
set yrange [-0.5:0.5]
set pointintervalbox 0.2
set grid
plot "file_$i" u 2:8  t 'hu at t=$t' with linespoints lc rgb "forest-green"
EOF
done
for i in $(seq 1 $tSteps); do 
	convert -density 1000 h$i.eps -quality 100 h$i.jpg;
	convert -density 1000 hu$i.eps -quality 100 hu$i.jpg;
done
rm $file.avi
ffmpeg -r 5 -i h%d.jpg h$file.avi
ffmpeg -r 5 -i hu%d.jpg hu$file.avi

mkdir picsNfiles
mv *jpg *.eps file_* picsNfiles
}

rm *.r *.u
c++ *.cpp
./a.out
ls *.u>1dFileNames
sed -i 's/..$//' 1dFileNames
while IFS= read -r line; do
    rm -r $line
    mkdir $line
    cp $line.* $line
    cd $line
    if [ ! -d picsNfiles ]; then
        postProcess;
    fi
    cd ..
done < 1dFileNames

