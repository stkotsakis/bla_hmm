set i=$1
set class=$3
set results=$4

while ($i <= $2)
csh meta_analysis_multi.csh $i $class $results > "$results"_class_"$class"/meta_analysis_multi_files/out.$i &
@ i++
end

wait

