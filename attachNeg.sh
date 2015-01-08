NEG_PATH="/home/tbepler/gordan_lab/ets/ETS_customarray_analysis/neg_ctrl/"
BOUND_PATH="/home/tbepler/gordan_lab/ets/ETS_customarray_analysis/bound1_oriented/"

for f in $@; do

    negtest="${NEG_PATH}test/${f}_negctrl_avgReps_lowAff_test.txt"
    negtrain="${NEG_PATH}train/${f}_negctrl_avgReps_lowAff_train.txt"
    boundtest="${BOUND_PATH}test/${f}_bound_log_combinedReps_fwd_test.txt"
    boundtrain="${BOUND_PATH}train/${f}_bound_log_combinedReps_fwd_train.txt"
    outtest="test/${f}_boundNeg_log_combinedReps_fwd_test.txt"
    outtrain="train/${f}_boundNeg_log_combinedReps_fwd_train.txt"

    cmd="cat ${negtest} ${boundtest}"
    echo "$cmd" ">" "$outtest"
    $cmd > $outtest

    cmd="cat ${negtrain} ${boundtrain}"
    echo "$cmd" ">" "$outtrain"
    $cmd > $outtrain


done
