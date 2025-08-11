set i=$1
set domain=$3
set assembly_summary=$4
set curr_date=$5
mkdir "$curr_date"_logs
while ($i <= $2)
csh download_and_bla_hmm.csh $i $domain $assembly_summary > "$curr_date"_logs/out.$i &
@ i++
end

wait


