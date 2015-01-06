PWM=$1
FILES=${@:2}

for F in $FILES; do
    
    NAME=$(basename $F)
    EXT=${NAME##*.}
    NAME=${NAME%.*}

    ORNT="${NAME}_oriented.${EXT}"
    FWD="${NAME}_fwd.${EXT}"
    RVS="${NAME}_rvs.${EXT}"

    orient -l -e -p $PWM -s $F -f $FWD -r $RVS

    #grep 'F' $ORNT > $FWD
    #grep 'R' $ORNT > $RVS

done
