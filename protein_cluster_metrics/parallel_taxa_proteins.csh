set i=$1
set class=$3
set results=$4

while ($i <= $2)
csh taxa_multi_proteins.csh $i $class $results > "$results"_class_"$class"/cluster_proteins_taxa/out.$i &
@ i++
end

wait

cat "$results"_class_"$class"/cluster_proteins_taxa/assembly_taxa_proteins_clusters.dat.* > "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat


