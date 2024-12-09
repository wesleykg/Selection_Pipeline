dirs=$1

while read dir; do
    cd $dir    
    echo "Working on $dir"
    python ../../scripts/1b_alignment_concatenator.py
    cd ..
done < $dirs