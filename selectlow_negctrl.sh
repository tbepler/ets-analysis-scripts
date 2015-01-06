for f in $@; do
    ext=${f##*.}
    name=${f%.*}
    out="${name}_lowAff.${ext}"

    cmd="tail -n 300 $f"
    echo "$cmd > $out"
    $cmd > "$out"
    
done
