for gzFile in $1*gz
do
echo gzip -d $gzFile
gzip -d $gzFile
echo python converGzArrayToNpz.py ${gzFile/.gz/}
python convertGzArrayToNpz.py ${gzFile/.gz/}
echo rm  ${gzFile/.gz/}
rm  ${gzFile/.gz/}
done

for gzFile in $1*npy
do
echo python convertGzArrayToNpz.py ${gzFile/.gz/}
python convertGzArrayToNpz.py ${gzFile/.gz/}
done
