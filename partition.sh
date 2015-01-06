NUM=$1

mkdir -p "train"
mkdir -p "test"

for f in ${@:2}; do

    filename=$(basename "$f")
    ext="${filename##*.}"
    filename="${filename%.*}"

    TRAIN="train/${filename}_train.${ext}"
    TEST="test/${filename}_test.${ext}"

    SHUF="shuf -n ${NUM} ${f}"
    echo "${SHUF} > ${TEST}"
    $SHUF > "$TEST"

    GREP="comm -23"
    echo "$GREP  <(sort ${f}) <(sort ${TEST}) > ${TRAIN}"
    $GREP  <(sort ${f}) <(sort "${TEST}") > "${TRAIN}"

done
