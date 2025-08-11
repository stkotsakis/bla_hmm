set class=$1
set results=$2
set cpus=$3

set unranked=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta.dat | awk '{print $6}' | sed 's/-/ /g' | awk '{for(i=1;i<=NF;i++){if($i ~ "group"){print $i}}}' | sort -u`

cp "$results"_class_"$class"/"$class"_whole_proteins_positives_meta.dat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat

foreach i ($unranked)
 sed -i 's/'$i'-//g' "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat
end



mkdir "$results"_class_"$class"/cluster_proteins_taxa

set clusters=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk '{print $33}' | sort -u | wc -l`
cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk '{print $33}' | sort -u > "$results"_class_"$class"/cluster_proteins_taxa/cluster_proteins.dat

cd "$results"_class_"$class"/cluster_proteins_taxa
set per_cpu=`echo $clusters | awk -v cpus=$cpus '{print int($1/cpus)}'`
awk -v a=$per_cpu 'NR%a==1{x="ASS"++i;}{print > x}' cluster_proteins.dat
set file_no=`ls ASS* | wc -l`
cd ../..

csh parallel_taxa_proteins.csh 1 $file_no $class $results





 

     
