
for f in $@; do

    ext="${f##*.}"
    filename="${f%.*}"

    EXE="/home/tbepler/gordan_lab/ets/ETS_customarray_analysis/scripts/cpp/featurize"

    OUT="${filename}_1mer.${ext}"
    FEAT="${EXE} -k 1 -i $f -o $OUT"
    echo "$FEAT"
    $FEAT

    OUT="${filename}_1-2mer.${ext}"
    FEAT="${EXE} -k 1 2 -i $f -o $OUT"
    echo "$FEAT"
    $FEAT

    OUT="${filename}_1-3mer.${ext}"
    FEAT="${EXE} -k 1 2 3 -i $f -o $OUT"
    echo "$FEAT"
    $FEAT

done
