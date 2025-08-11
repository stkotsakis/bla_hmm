# PROTEINS AND CDS EXTRACTION FOR WHOLE GENES, CLUSTERING

setenv PATH $PATH\:$HOME/signalp-5.0b/bin

set classes=(A B1 B2 B3 C D)
set results=$1
set cpus=$2
foreach class ($classes)


mkdir "$results"_class_"$class"

cat $results | awk -v cl=$class '$8 ~ "class_"cl' | awk -v cl=$class '{print cl"_"NR,$0}' > "$results"_class_"$class"/"$results"_class_"$class"_wcodes.dat # FILE 1

cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes.dat | awk '$17 !~ "YYYYYY" && $17 !~ "PSEUDO" && $29 == 0 && $26 >= 60 && $19 !~ "X" && $20 !~ "PSEUDO" {print $1}' > "$results"_class_"$class"/"$class"_whole_proteins_code.dat # FILE 2

cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes.dat | awk '$17 !~ "YYYYYY" && $17 !~ "PSEUDO" && $29 == 0 && $26 >= 60 && $19 !~ "X" && $20 !~ "PSEUDO" {print ">"$1" ""\n"$19}' > "$results"_class_"$class"/"$class"_whole_proteins.fasta # FILE 3

cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes.dat | awk '$17 !~ "YYYYYY" && $17 !~ "PSEUDO" && $29 == 0 && $26 >= 60 && $19 !~ "X" && $20 !~ "PSEUDO" {print ">"$1" ""\n"$16}' > "$results"_class_"$class"/"$class"_whole_cds.fasta # FILE 4

# cds clustering and codon alignment

cd-hit -i "$results"_class_"$class"/"$class"_whole_cds.fasta -o "$results"_class_"$class"/clustered_class_"$class"_whole_cds -c 0.8 -s 0.9 -d 0 -T 8 # FILES 5 AND 6

# extract clustering information in usable format

mkdir "$results"_class_"$class"/temp

cat "$results"_class_"$class"/clustered_class_"$class"_whole_cds.clstr | sed 's/>Cluster />Cluster_/g' > "$results"_class_"$class"/temp/clustered_class_"$class"_whole_cds.clstr

set clusters=`cat "$results"_class_"$class"/temp/clustered_class_"$class"_whole_cds.clstr | sed 's/>Cluster />Cluster_/g'  | awk '$1 ~ "Cluster_" {print $1}' | sed 's/>//g'`

echo ">Cluster_" >> "$results"_class_"$class"/temp/clustered_class_"$class"_whole_cds.clstr

foreach i ($clusters)
 cat "$results"_class_"$class"/temp/clustered_class_"$class"_whole_cds.clstr | sed 's/\<'$i'\>/'$i' /g' | sed -n '/>'$i' /,/>Cluster_/p' | sed '/>Cluster_/d' | sed 's/\.\.\.//g' > "$results"_class_"$class"/temp/"$i".temp
 set repr=`cat "$results"_class_"$class"/temp/"$i".temp | awk '$4 ~ "*" {print $3}'`
 set members=`cat "$results"_class_"$class"/temp/"$i".temp | awk '{print $3}'`
 foreach j ($members)
  set j=`echo $j | sed 's/>//g'`
  echo $i "$repr" $j >> "$results"_class_"$class"/class_"$class"_clusters_cds.dat # FILE 7
 end
 sed -i 's/'$repr' />'$i' /g' "$results"_class_"$class"/clustered_class_"$class"_whole_cds # CHANGE ENTRY NAMES IN FILE 5
end
rm -rf "$results"_class_"$class"/temp

# For CDS only:check relation of each cluster to AMR ref db and update the _clusters_cds.dat file

cat "$results"_class_"$class"/clustered_class_"$class"_whole_cds  | sed 's/Cluster_/'$class'_clbardb_/g' > "$results"_class_"$class"/clustered_class_"$class"_whole_cds.temp
cat "$results"_class_"$class"/clustered_class_"$class"_whole_cds.temp ncbi_blas/cdhit_clustered_class_"$class"_cds > "$results"_class_"$class"/"$class"_combined_w_refdb.fasta
cd-hit -i "$results"_class_"$class"/"$class"_combined_w_refdb.fasta -o "$results"_class_"$class"/"$class"_combined_w_refdb_clustered -c 0.8 -s 0.9 -d 0 -T 8
mkdir "$results"_class_"$class"/temp
cat "$results"_class_"$class"/"$class"_combined_w_refdb_clustered.clstr | sed 's/>Cluster />Cluster_/g' > "$results"_class_"$class"/temp/"$class"_combined_w_refdb_clustered.clstr
set clusters=`cat "$results"_class_"$class"/temp/"$class"_combined_w_refdb_clustered.clstr | sed 's/>Cluster />Cluster_/g'  | awk '$1 ~ "Cluster_" {print $1}' | sed 's/>//g'`
echo ">Cluster_" >> "$results"_class_"$class"/temp/"$class"_combined_w_refdb_clustered.clstr


foreach m ($clusters)
 cat "$results"_class_"$class"/temp/"$class"_combined_w_refdb_clustered.clstr | sed 's/\<'$m'\>/'$m' /g' | sed -n '/>'$m' /,/>Cluster_/p' | sed '/>Cluster_/d' | sed 's/\.\.\.//g' > "$results"_class_"$class"/temp/"$m".temp
 set members=`cat "$results"_class_"$class"/temp/"$m".temp | awk '{print $1}' | wc -l`
 if ($members > 1) then
  set clus=`cat "$results"_class_"$class"/temp/"$m".temp | awk -v cl=$class '$3 ~ cl"_clbardb_" {print $3}' | sed 's/>//g' | sed 's/'$class'_clbardb_/Cluster_/g'`
  set ref=`cat "$results"_class_"$class"/temp/"$m".temp | awk -v cl=$class '$3 !~ cl"_clbardb_" {print $3}' | sed 's/>//g'`
  foreach l ($clus)
  cat "$results"_class_"$class"/class_"$class"_clusters_cds.dat | awk -v cl=$l -v ref="$ref" '$1 == cl {print $0,ref}' >> "$results"_class_"$class"/class_"$class"_clusters_cds_wrefs.dat
  end
 else
  set entry=`cat "$results"_class_"$class"/temp/"$m".temp | awk '{print $3}' | sed 's/>//g'`
  set cl_check=`echo $entry | grep _clbardb_ | wc -l`
  if ($cl_check != "0") then
   set clus=`echo $entry | sed 's/'$class'_clbardb_/Cluster_/g'`
   set ref="CLUSTER_NOT_IN_NCBI_AMR_REF_DB"
   cat "$results"_class_"$class"/class_"$class"_clusters_cds.dat | awk -v cl=$clus -v ref="$ref" '$1 == cl {print $0,ref}' >> "$results"_class_"$class"/class_"$class"_clusters_cds_wrefs.dat
  else
   echo $entry >> "$results"_class_"$class"/amr_rf_db_class_"$class"_unclustered.dat
  endif
 endif
end
  
rm -rf "$results"_class_"$class"/temp



#for codon alignment: remove stop codons from the end of cds. 

mkdir "$results"_class_"$class"/temp
set cds_entries=`cat "$results"_class_"$class"/clustered_class_"$class"_whole_cds | awk '$1 ~ ">" {print $1}' | sed 's/>//g'`
cp "$results"_class_"$class"/clustered_class_"$class"_whole_cds "$results"_class_"$class"/temp/
echo ">" >> "$results"_class_"$class"/temp/clustered_class_"$class"_whole_cds 
foreach o ($cds_entries)
 cat "$results"_class_"$class"/temp/clustered_class_"$class"_whole_cds | sed -n '/>'$o' /,/>/p' | sed '/>/d' | tr -d '\n' | rev | cut -c4- | rev | sed '1s/^/>'$o' \n/' > "$results"_class_"$class"/temp/"$o".temp
end

cat "$results"_class_"$class"/temp/*.temp > "$results"_class_"$class"/clustered_class_"$class"_whole_cds_no_stop

rm -rf "$results"_class_"$class"/temp


# protein clustering

cd-hit -i "$results"_class_"$class"/"$class"_whole_proteins.fasta -o "$results"_class_"$class"/clustered_class_"$class"_whole_proteins -c 0.8 -s 0.9 -d 0 -T 8 # FILES 8 AND 9

mkdir "$results"_class_"$class"/temp

cat "$results"_class_"$class"/clustered_class_"$class"_whole_proteins.clstr | sed 's/>Cluster />Cluster_/g' > "$results"_class_"$class"/temp/clustered_class_"$class"_whole_proteins.clstr

set clusters_p=`cat "$results"_class_"$class"/temp/clustered_class_"$class"_whole_proteins.clstr | sed 's/>Cluster />Cluster_/g'  | awk '$1 ~ "Cluster_" {print $1}' | sed 's/>//g'`

echo ">Cluster_" >> "$results"_class_"$class"/temp/clustered_class_"$class"_whole_proteins.clstr
# extract clustering information in usable format

foreach k ($clusters_p)
 cat "$results"_class_"$class"/temp/clustered_class_"$class"_whole_proteins.clstr | sed 's/\<'$k'\>/'$k' /g' | sed -n '/>'$k' /,/>Cluster_/p' | sed '/>Cluster_/d' | sed 's/\.\.\.//g' > "$results"_class_"$class"/temp/"$k".temp
 set repr_p=`cat "$results"_class_"$class"/temp/"$k".temp | awk '$4 ~ "*" {print $3}'`
 set members_p=`cat "$results"_class_"$class"/temp/"$k".temp | awk '{print $3}'`
 foreach l ($members_p)
  set l=`echo $l | sed 's/>//g'`
  echo $k "$repr_p" $l >> "$results"_class_"$class"/class_"$class"_clusters_proteins.dat # FILE 10
 end
 sed -i 's/'$repr_p' />'$k' /g' "$results"_class_"$class"/clustered_class_"$class"_whole_proteins # CHANGE ENTRY NAMES IN FILE 8
end

rm -rf "$results"_class_"$class"/temp


# signal peptide predictions

cd "$results"_class_"$class"

signalp -org gram- -prefix gramneg -fasta "$class"_whole_proteins.fasta

signalp -org gram+ -prefix grampos -fasta "$class"_whole_proteins.fasta

cd ..

# add protein lengths

cat "$results"_class_"$class"/"$results"_class_"$class"_wcodes.dat | awk '{print $0,length($19)}' > "$results"_class_"$class"/"$results"_class_"$class"_wcodes_plen.dat

end

foreach class ($classes)
set entries_whole=`cat "$results"_class_"$class"/"$class"_whole_proteins_code.dat | awk '{print $1}' | sort -u | wc -l`
# invoke meta_analysis script if entries >= 10000 run meta-analysis in parallel
if ($entries_whole >= 10000) then 
 mkdir "$results"_class_"$class"/meta_analysis_multi_files
 cp "$results"_class_"$class"/"$class"_whole_proteins_code.dat "$results"_class_"$class"/meta_analysis_multi_files/"$class"_whole_proteins_code.dat
 cd "$results"_class_"$class"/meta_analysis_multi_files
 set per_cpu=`echo $entries_whole | awk -v cpus=$cpus '{print int($1/cpus)}'`
 awk -v a=$per_cpu 'NR%a==1{x="ASS"++i;}{print > x}' "$class"_whole_proteins_code.dat
 set file_no=`ls ASS* | wc -l`
 cd ../..
 csh parallel_multi.csh 1 $file_no $class $results
 cat "$results"_class_"$class"/meta_analysis_multi_files/"$class"_whole_proteins_positives_meta.dat.* > "$results"_class_"$class"/"$class"_whole_proteins_positives_meta.dat
 rm -rf meta_analysis_multi_files
else
 csh meta_analysis.csh $class $results
endif
end





