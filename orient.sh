PWM=$1
FILES=${@:2}

for F in $FILES; do
    
    NAME=$(basename $F)
    EXT=${NAME##*.}
    NAME=${NAME%.*}

    ORNT="${NAME}_oriented.${EXT}"
    FWD="${NAME}_fwd.${EXT}"
    RVS="${NAME}_rvs.${EXT}"

    EXE="/home/tbepler/gordan_lab/ets/ETS_customarray_analysis/scripts/cpp/orient" 

    $EXE -l -e -p $PWM -s $F -f $FWD -r $RVS

    #grep 'F' $ORNT > $FWD
    #grep 'R' $ORNT > $RVS

done
