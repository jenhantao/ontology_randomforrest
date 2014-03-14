treeFiles=()
for treeFileSplit in $1*txt
do
treeFiles+=($treeFileSplit)
done

# loop over tree files and compute distances
for treeFile in ${treeFiles[@]}
do
pythonOut=${treeFile/.txt/_distance.npy}
if [ ! -f $pythonOut ]
then
echo python treeDistanceCounter.py $treeFile $2
# submit jobs to cluster
qsub python treeDistanceCounter.py $treeFile
fi
done

# call script to sum individual distances and compute final matrix

