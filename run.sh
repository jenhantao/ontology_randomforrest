treeFiles=()
for treeFileSplit in $1*txt
do
treeFiles+=($treeFileSplit)
done

# loop over tree files and compute distances
length=${#treeFiles[*]}
#for treeFile in ${treeFiles[@]}
for (( i=0; i<$(( $length )); i++ ))
do
if [ $(expr $i % 8) = 0 ]
then
sleep 150
fi
treeFile="${treeFiles[$i]}"
pythonOut=${treeFile/.txt/_distance.npy}
if [ ! -f $pythonOut ]
then
echo python treeDistanceCounter.py $treeFile $2
python treeDistanceCounter.py $treeFile $2 &
fi
done

# call script to sum individual distances and compute final matrix

