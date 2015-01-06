

KS="1mer
1-2mer
1-3mer"

TFS="ELK1_100nM
ELK1_50nM
ELK1_10nM"

for TF in $TFS
do

    for K in $KS
    do

        M=models/bayesian-regression/all_tfs_2014_12_04/${TF}_bound_log_fwd*${K}_weightsMean.txt
        P="../../chip-seq/idr_peaks/elk1AndNonOverlappingDnaseLabeled.seqs"
        O="roc/elk1DnaseLabeled_${TF}_Bayesian${K}Scored.txt"
        echo ./scoreSeq.sh $M $P ">" $O
        ./scoreSeq.sh $M $P > $O

    done

done

TFS="ETS1_100nM
ETS1_10nM"

for TF in $TFS
do

    for K in $KS
    do

        M=models/bayesian-regression/all_tfs_2014_12_04/${TF}_bound_log_fwd*${K}_weightsMean.txt
        P="../../chip-seq/idr_peaks/ets1AndNonOverlappingDnaseLabeled.seqs"
        O="roc/ets1DnaseLabeled_${TF}_Bayesian${K}Scored.txt"
        echo ./scoreSeq.sh $M $P ">" $O
        ./scoreSeq.sh $M $P > $O

    done

done

TFS="GABPA_100nM
GABPA_50nM
GABPA_25nM"

for TF in $TFS
do

    for K in $KS
    do

        M=models/bayesian-regression/all_tfs_2014_12_04/${TF}_bound_log_fwd*${K}_weightsMean.txt
        P="../../chip-seq/idr_peaks/gabpaAndNonOverlappingDnaseLabeled.seqs"
        O="roc/gabpaDnaseLabeled_${TF}_Bayesian${K}Scored.txt"
        echo ./scoreSeq.sh $M $P ">" $O
        ./scoreSeq.sh $M $P > $O

    done

done
