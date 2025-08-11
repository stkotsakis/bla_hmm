set class=$1
set results=$2

set unranked=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta.dat | awk '{print $6}' | sed 's/-/ /g' | awk '{for(i=1;i<=NF;i++){if($i ~ "group"){print $i}}}' | sort -u`

cp "$results"_class_"$class"/"$class"_whole_proteins_positives_meta.dat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat

foreach i ($unranked)
 sed -i 's/'$i'-//g' "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat
end

set clusters_proteins=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk '{print $33}' | sort -u`

foreach i ($clusters_proteins)
 set assemblies=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $2}' | sort -u`
 foreach j ($assemblies)
  set genus=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $6}' | sort -u`
  set family_check=`echo $genus | awk '$1 ~ "ceae"'| wc -l`
  if ($family_check != "0") then
   set genus=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $7}' | sort -u `
  endif
  if ($#genus == 0) then
   set genus=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $5}' | sort -u `
   if ($#genus == 0) then
    set genus=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $4}' | sort -u `
    if ($#genus == 0) then
     set genus=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $3}' | sort -u `
     if ($#genus == 0) then
      set genus=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $2}' | sort -u`
      if ($#genus == 0) then
       set genus=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $1}' | sort -u`
      endif
     endif
    endif
   endif
  endif
  set family=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $5}' | sort -u`
  if ($#family == 0) then
   set family=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $4}' | sort -u`
   if ($#family == 0) then
    set family=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $3}' | sort -u`
    if ($#family == 0) then
     set family=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $2}' | sort -u`
     if ($#family == 0) then
      set family=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $1}' | sort -u`
     endif
    endif
   endif
  endif
  set order=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $4}' | sort -u`
   if ($#order == 0) then
    set order=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $3}' | sort -u`
    if ($#order == 0) then
     set order=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $2}' | sort -u`
      if ($#order == 0) then
       set order=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $1}' | sort -u`
      endif
     endif
    endif
  set classt=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $3}' | sort -u`
  if ($#classt == 0) then
   set classt=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $2}' | sort -u`
   if ($#classt == 0) then
    set classt=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $1}' | sort -u`
   endif
  endif
  set phylum=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $2}' | sort -u`
   if ($#phylum == 0) then
    set phylum=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v ass=$j '$2 == ass {print $6}' | sed 's/-/ /g' | awk '{print $1}' | sort -u`
   endif
  echo $i $j "$genus" "$family" "$order" "$classt" "$phylum" >> "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat
 end
end



 

     
