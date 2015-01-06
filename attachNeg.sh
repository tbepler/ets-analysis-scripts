NEG_PATH="/home/tbepler/gordan_lab/ets/ETS_customarray_analysis/neg_ctrl/"

for f in $@; do

    bname=$(basename $f)
    tf=${bname%%_*_*}
    echo $tf


done
