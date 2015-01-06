
for f in $@; do

    ext="${f##*.}"
    filename="${f%.*}"

    OUT="${filename}_1mer.${ext}"
    FEAT="featurize -k 1 -i $f -o $OUT"
    echo "$FEAT"
    $FEAT

    OUT="${filename}_1-2mer.${ext}"
    FEAT="featurize -k 1 2 -i $f -o $OUT"
    echo "$FEAT"
    $FEAT

    OUT="${filename}_1-3mer.${ext}"
    FEAT="featurize -k 1 2 3 -i $f -o $OUT"
    echo "$FEAT"
    $FEAT

done
