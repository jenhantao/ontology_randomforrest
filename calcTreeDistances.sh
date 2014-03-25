# arguments: path to trees, path to ontology (term_term)
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
treeFile="${treeFiles[$i]}"
pythonOut=${treeFile/.txt/.npz}
if [ ! -f $pythonOut ]
then
	if [ $(expr $i % 8) = 0 ]
	then
	echo "napping"
	sleep 120
	fi
echo python treeDistanceCounter.py $treeFile $2
python treeDistanceCounter.py $treeFile $2 &
fi
done


# move files to dropbox if all files are complete
echo "moving files"
# give time for lst runs to finish
#sleep 150
#for npzFile in $1*npz
#do
#echo mv $npzFile ~/Dropbox/random_forrest/$npzFile
#mv $npzFile ~/Dropbox/random_forrest/$npzFile
#done


# call script to sum individual distances and compute final matrix


