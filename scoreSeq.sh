
tmp=.tmp$(date +"%s")
transpose "$1" > "$tmp"
./scoreSeq.out "$tmp" "$2"
rm "$tmp"
