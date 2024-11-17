for dir in ../alignments/*/;
do
    cd $dir
    for alignment in *.fasta;
    do
    	echo "Working on" $alignment "in" $dir
        trimal -in $alignment -out $alignment.noallgaps -noallgaps -keepseqs
        filename="$alignment.noallgaps"
        if [ -e "$filename" ]; then
            sed 's/ [0-9]* bp//g' $alignment.noallgaps > $alignment
            rm $alignment.noallgaps
        fi
    done
    cd ../../scripts/
done