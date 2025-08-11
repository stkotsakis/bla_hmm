set class=$1
set results=$2

set clusters_proteins=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk '{print $33}' | sort -u`

foreach i ($clusters_proteins)
 set members=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $1}' | wc -l`
 set species=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $5}' | sort -u`
 set species=`echo $species | sed 's/ /,/g'`
 set spec_taxids_num=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $4}' | sort -u | wc -l`
 set spec_taxids=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $4}' | sort -u`
 set spec_taxids=`echo $spec_taxids | sed 's/ /,/g'`
 set genera=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $3}' | sort -u`
 set genera=`echo $genera | sed 's/ /,/g'`
 set genus_num=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $3}' | sort -u | wc -l`
 set families=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $4}' | sort -u`
 set families=`echo $families | sed 's/ /,/g'`
 set family_num=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $4}' | sort -u | wc -l`
 set orders=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $5}' | sort -u`
 set orders=`echo $orders | sed 's/ /,/g'`
 set order_num=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $5}' | sort -u | wc -l`
 set classes=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $6}' | sort -u`
 set classes=`echo $classes | sed 's/ /,/g'`
 set class_num=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $6}' | sort -u | wc -l`
 set phyla=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $7}' | sort -u`
 set phyla=`echo $phyla | sed 's/ /,/g'`
 set phylum_num=`cat "$results"_class_"$class"/assembly_taxa_proteins_clusters.dat | awk -v clu=$i '$1 == clu {print $7}' | sort -u | wc -l` 
 if ($spec_taxids_num == "1") then
  set taxon_label=$spec_taxids
 else 
  if ($genus_num == "1") then
   set taxon_label=$genera
  else
   if ($family_num == "1") then
    set taxon_label=$families
   else
    if ($order_num == "1") then
     set taxon_label=$orders
    else
     if ($class_num == "1") then
      set taxon_label=$classes
     else
      if ($phylum_num == "1") then
       set taxon_label=$phyla
      else
       set taxon_label="MULTI_PHYLA"
      endif
     endif
    endif
   endif
  endif
 endif
 set cds_clusters_num=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $31}' | sort -u |  wc -l`
 set cds_clusters=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $31}' | sort -u`
 set cds_clusters=`echo $cds_clusters | sed 's/ /,/g'`
 set motif1=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $38}' | sort -u`
 set motif1=`echo $motif1 | sed 's/ /,/g'`
 set motif2=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $39}' | sort -u`
 set motif2=`echo $motif2 | sed 's/ /,/g'` 
 set motif3=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $40}' | sort -u`
 set motif3=`echo $motif3 | sed 's/ /,/g'` 
 set motif4=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $41}' | sort -u`
 set motif4=`echo $motif4 | sed 's/ /,/g'` 
 set secretion_neg=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $34}' | sort -u`
 set secretion_neg=`echo $secretion_neg | sed 's/ /,/g'`
 set secretion_pos=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu {print $36}' | sort -u`
 set secretion_pos=`echo $secretion_pos | sed 's/ /,/g'`
 set ref_db_check=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu && $32 !~ "CLUSTER_NOT_IN_NCBI_AMR_REF_DB"' | wc -l`
 if ($ref_db_check != "0") then
  set experimental_evidence="YES"
 else
  set experimental_evidence="NO"
 endif
 set comment="."
 if ($class == "D") then
  set true_hit_check=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu && $21 ~ "BlaR1"' | wc -l`
  if ($true_hit_check != "0") then
   set comment="BlaR1"
  endif
 endif
 if ($class == "B3") then
  set true_hit_check=`cat "$results"_class_"$class"/"$class"_whole_proteins_positives_meta_ranked.dat | awk -v clu=$i '$33 == clu && $21 ~ "HARLDQ_MOTIF"' | wc -l`
  if ($true_hit_check != "0") then
   set comment="HARLDQ_MOTIF"
  endif
 endif
 echo $i $members $spec_taxids_num "$species" "$spec_taxids" $genus_num "$genera" $family_num "$families" $order_num "$orders" $class_num "$classes" $phylum_num "$phyla" $taxon_label $cds_clusters_num $cds_clusters $motif1 $motif2 $motif3 $motif4 $secretion_neg $secretion_pos $experimental_evidence $comment >> "$results"_class_"$class"/data_protein_clusters.dat
end


 

     
