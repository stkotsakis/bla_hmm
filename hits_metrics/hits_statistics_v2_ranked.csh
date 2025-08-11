#beta-lactamase statistics carriage foreach taxonomic level true hits/tru hits secreted

set results=$1 #file used as input for cluster analysis
set output=$2

rm -rf $output
mkdir $output

set blas=(A C D B1 B2 B3)

sed -i 's/(//g' $results
sed -i 's/)//g' $results
sed -i 's/cellular_organisms-//g' $results
sed -i 's/\//_/g' $results
sed -i 's/in:/in/g' $results
sed -i 's/,_/_/g' $results
sed -i 's/\._/_/g' $results
sed -i 's/Gram-positive/Gram_positive/g' $results
sed -i 's/NO_LINEAGE/NO_LINEAGE-/g' $results
sed -i 's/G+C/G_C/g' $results

set unranked=`cat $results | awk '{print $5}' | sed 's/-/ /g' | awk '{for(i=1;i<=NF;i++){if($i ~ "group"){print $i}}}' | sort -u`

cp $results "$results"_ranked.dat

foreach i ($unranked)
 sed -i 's/'$i'-//g' "$results"_ranked.dat
end

foreach r ( $blas )
 sed -i 's/(//g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/)//g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/cellular_organisms-//g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/\//_/g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/in:/in/g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/,_/_/g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/\._/_/g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/Gram-positive/Gram_positive/g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/NO_LINEAGE/NO_LINEAGE-/g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/G+C/G_C/g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat
 sed -i 's/(//g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/)//g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/cellular_organisms-//g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/\//_/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/in:/in/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/,_/_/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/\._/_/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/Gram-positive/Gram_positive/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/NO_LINEAGE/NO_LINEAGE-/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/G+C/G_C/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes.dat
 sed -i 's/(//g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 sed -i 's/)//g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 sed -i 's/cellular_organisms-//g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 sed -i 's/\//_/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 sed -i 's/in:/in/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 sed -i 's/,_/_/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 sed -i 's/\._/_/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 sed -i 's/Gram-positive/Gram_positive/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 sed -i 's/NO_LINEAGE/NO_LINEAGE-/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 sed -i 's/G+C/G_C/g' "$results"_class_"$r"/"$results"_class_"$r"_wcodes_plen.dat
 cp "$results"_class_"$r"/"$r"_whole_proteins_positives_meta.dat "$results"_class_"$r"/"$r"_whole_proteins_positives_meta_ranked.dat
 foreach b ( $unranked )
  sed -i 's/'$i'-//g' "$results"_class_"$r"/"$r"_whole_proteins_positives_meta_ranked.dat
 end
 cp "$results"_class_"$r"/"$r"_whole_proteins_positives_meta_ranked.dat $output/"$r"_whole_proteins_positives_meta_ranked.dat
end  



#taxonomic levels highest to lowest first = domain, second = phylum ....
set first_tl=`cat "$results"_ranked.dat | awk '{print $5}' | sort -u | awk '{print $1}' | sed 's/-/ /g' | awk '{print $1}' | sort -u`

foreach i ( $first_tl )
 mkdir $output/$i
 cat "$results"_ranked.dat | awk -v otu=$i '$5 ~ otu"-"' > "$output"/$i/"$i"_total.dat
 foreach r ( $blas )
  cat $output/"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$i '$6 ~ otu"-"' > "$output"/$i/"$i"_"$r"_whole_proteins_positives_meta_ranked.dat
  cat $output/"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$i '$6 ~ otu"-" {print $1}' > "$output"/$i/"$i"_"$r"_whole_proteins_positives_codes.dat
  csh evaluate_count_ranked.csh "$output"/$i/"$i"_total.dat $output $r $i
  cat "$output"/"$i"/"$i"_class_"$r"_statistics.out | awk -v otu=$i '{print "FIRST",otu,$0}' >> "$output"/class_"$r"_statistics.dat
 end  
#---------------------------------------------------END OF FIRST OTU BLA EVAL-------------------------------------------------------------------------------------------------
 set second_tl=`cat "$results"_ranked.dat | awk '{print $5}' | sort -u | awk '{print $1}' | sed 's/-/ /g' | awk -v first=$i '$1 == first {print $2}' | sort -u`
 foreach j ( $second_tl )
  mkdir "$output"/$i/$j
  cat "$output"/$i/"$i"_total.dat | awk -v otu=$j '$5 ~ otu"-"' > "$output"/$i/$j/"$j"_total.dat
  foreach r ( $blas )
   cat "$output"/$i/"$i"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$j '$6 ~ otu"-"' > "$output"/$i/$j/"$j"_"$r"_whole_proteins_positives_meta_ranked.dat
   cat "$output"/$i/$j/"$j"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$j '$6 ~ otu"-" {print $1}' > "$output"/$i/$j/"$j"_"$r"_whole_proteins_positives_codes.dat
   csh evaluate_count_ranked.csh "$output"/$i/$j/"$j"_total.dat $output $r $j $i/
   cat "$output"/$i/$j/"$j"_class_"$r"_statistics.out | awk -v otu=$j '{print "SECOND",otu,$0}' >> "$output"/class_"$r"_statistics.dat
  end  
#--------------------------------------------------END OF SECOND OTU BLA EVAL---------------------------------------------------------------------------------------------------
  set third_tl=`cat "$results"_ranked.dat | awk '{print $5}' | sort -u | awk '{print $1}' | sed 's/-/ /g' | awk -v first=$i -v second=$j '$1 == first && $2 == second {print $3}' | sort -u`
  foreach k ( $third_tl )
   mkdir "$output"/$i/$j/$k
   cat "$output"/$i/$j/"$j"_total.dat | awk -v otu=$k '$5 ~ otu"-"' > "$output"/$i/$j/$k/"$k"_total.dat
   foreach r ( $blas )
    cat "$output"/$i/$j/"$j"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$k '$6 ~ otu"-"' > "$output"/$i/$j/$k/"$k"_"$r"_whole_proteins_positives_meta_ranked.dat
    cat "$output"/$i/$j/$k/"$k"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$k '$6 ~ otu"-" {print $1}' > "$output"/$i/$j/$k/"$k"_"$r"_whole_proteins_positives_codes.dat
    csh evaluate_count_ranked.csh "$output"/$i/$j/$k/"$k"_total.dat $output $r $k $i/$j/
    cat "$output"/$i/$j/$k/"$k"_class_"$r"_statistics.out | awk -v otu=$k '{print "THIRD",otu,$0}' >> "$output"/class_"$r"_statistics.dat
   end  
#--------------------------------------------------END OF THIRD OTU BLA EVAL----------------------------------------------------------------------------------------------------
   set forth_tl=`cat "$results"_ranked.dat | awk '{print $5}' | sort -u | awk '{print $1}' | sed 's/-/ /g' | awk -v first=$i -v second=$j -v third=$k '$1 == first && $2 == second && $3 == third {print $4}' | sort -u`
   foreach l ( $forth_tl )
    mkdir "$output"/$i/$j/$k/$l
    cat "$output"/$i/$j/$k/"$k"_total.dat | awk -v otu=$l '$5 ~ otu"-"' > "$output"/$i/$j/$k/$l/"$l"_total.dat
    foreach r ( $blas )
     cat "$output"/$i/$j/$k/"$k"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$l '$6 ~ otu"-"' > "$output"/$i/$j/$k/$l/"$l"_"$r"_whole_proteins_positives_meta_ranked.dat
     cat "$output"/$i/$j/$k/$l/"$l"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$l '$6 ~ otu"-" {print $1}' > "$output"/$i/$j/$k/$l/"$l"_"$r"_whole_proteins_positives_codes.dat
     csh evaluate_count_ranked.csh "$output"/$i/$j/$k/$l/"$l"_total.dat $output $r $l $i/$j/$k/
     cat "$output"/$i/$j/$k/$l/"$l"_class_"$r"_statistics.out | awk -v otu=$l '{print "FORTH",otu,$0}' >> "$output"/class_"$r"_statistics.dat
    end  
#--------------------------------------------------END OF FORTH OTU BLA EVAL----------------------------------------------------------------------------------------------------
    set fifth_tl=`cat "$results"_ranked.dat | awk '{print $5}' | sort -u | awk '{print $1}' | sed 's/-/ /g' | awk -v first=$i -v second=$j -v third=$k -v forth=$l '$1 == first && $2 == second && $3 == third && $4 == forth {print $5}' | sort -u`
    foreach m ( $fifth_tl )
     mkdir "$output"/$i/$j/$k/$l/$m
     cat "$output"/$i/$j/$k/$l/"$l"_total.dat | awk -v otu=$m '$5 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/"$m"_total.dat
     foreach r ( $blas )
      cat "$output"/$i/$j/$k/$l/"$l"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$m '$6 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/"$m"_"$r"_whole_proteins_positives_meta_ranked.dat
      cat "$output"/$i/$j/$k/$l/$m/"$m"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$m '$6 ~ otu"-" {print $1}' > "$output"/$i/$j/$k/$l/$m/"$m"_"$r"_whole_proteins_positives_codes.dat
      csh evaluate_count_ranked.csh "$output"/$i/$j/$k/$l/$m/"$m"_total.dat $output $r $m $i/$j/$k/$l/
      cat "$output"/$i/$j/$k/$l/$m/"$m"_class_"$r"_statistics.out | awk -v otu=$m '{print "FIFTH",otu,$0}' >> "$output"/class_"$r"_statistics.dat
     end 
#--------------------------------------------------END OF FIFTH OTU BLA EVAL----------------------------------------------------------------------------------------------------
     set sixth_tl=`cat "$results"_ranked.dat | awk '{print $5}' | sort -u | awk '{print $1}' | sed 's/-/ /g' | awk -v first=$i -v second=$j -v third=$k -v forth=$l -v fifth=$m '$1 == first && $2 == second && $3 == third && $4 == forth && $5 == fifth {print $6}' | sort -u`
     foreach n ( $sixth_tl )
      mkdir "$output"/$i/$j/$k/$l/$m/$n
      cat "$output"/$i/$j/$k/$l/$m/"$m"_total.dat | awk -v otu=$n '$5 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/$n/"$n"_total.dat
      foreach r ( $blas )
       cat "$output"/$i/$j/$k/$l/$m/"$m"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$n '$6 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/$n/"$n"_"$r"_whole_proteins_positives_meta_ranked.dat
       cat "$output"/$i/$j/$k/$l/$m/$n/"$n"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$n '$6 ~ otu"-" {print $1}' > "$output"/$i/$j/$k/$l/$m/$n/"$n"_"$r"_whole_proteins_positives_codes.dat
       csh evaluate_count_ranked.csh "$output"/$i/$j/$k/$l/$m/$n/"$n"_total.dat $output $r $n $i/$j/$k/$l/$m/
       cat "$output"/$i/$j/$k/$l/$m/$n/"$n"_class_"$r"_statistics.out | awk -v otu=$n '{print "SIXTH",otu,$0}' >> "$output"/class_"$r"_statistics.dat
      end 
#--------------------------------------------------END OF SIXTH OTU BLA EVAL----------------------------------------------------------------------------------------------------
      set seventh_tl=`cat "$results"_ranked.dat | awk '{print $5}' | sort -u | awk '{print $1}' | sed 's/-/ /g' | awk -v first=$i -v second=$j -v third=$k -v forth=$l -v fifth=$m -v sixth=$n '$1 == first && $2 == second && $3 == third && $4 == forth && $5 == fifth && $6 == sixth {print $7}' | sort -u`
      foreach o ( $seventh_tl )
       mkdir "$output"/$i/$j/$k/$l/$m/$n/$o
       cat "$output"/$i/$j/$k/$l/$m/$n/"$n"_total.dat | awk -v otu=$o '$5 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/$n/$o/"$o"_total.dat       
       foreach r ( $blas )
        cat "$output"/$i/$j/$k/$l/$m/$n/"$n"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$o '$6 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/$n/$o/"$o"_"$r"_whole_proteins_positives_meta_ranked.dat        
        cat "$output"/$i/$j/$k/$l/$m/$n/$o/"$o"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$o '$6 ~ otu"-" {print $1}' > "$output"/$i/$j/$k/$l/$m/$n/$o/"$o"_"$r"_whole_proteins_positives_codes.dat
        csh evaluate_count_ranked.csh "$output"/$i/$j/$k/$l/$m/$n/$o/"$o"_total.dat $output $r $o $i/$j/$k/$l/$m/$n/
        cat "$output"/$i/$j/$k/$l/$m/$n/$o/"$o"_class_"$r"_statistics.out | awk -v otu=$o '{print "SEVENTH",otu,$0}' >> "$output"/class_"$r"_statistics.dat
       end 
#--------------------------------------------------END OF SEVENTH OTU BLA EVAL----------------------------------------------------------------------------------------------------
       set eighth_tl=`cat "$results"_ranked.dat | awk '{print $5}' | sort -u | awk '{print $1}' | sed 's/-/ /g' | awk -v first=$i -v second=$j -v third=$k -v forth=$l -v fifth=$m -v sixth=$n -v seventh=$o '$1 == first && $2 == second && $3 == third && $4 == forth && $5 == fifth && $6 == sixth && $7 == seventh {print $8}' | sort -u`
       foreach p ( $eighth_tl )
        mkdir "$output"/$i/$j/$k/$l/$m/$n/$o/$p
        cat "$output"/$i/$j/$k/$l/$m/$n/$o/"$o"_total.dat | awk -v otu=$p '$5 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/$n/$o/$p/"$p"_total.dat       
        foreach r ( $blas )
         cat "$output"/$i/$j/$k/$l/$m/$n/$o/"$o"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$p '$6 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/$n/$o/$p/"$p"_"$r"_whole_proteins_positives_meta_ranked.dat
         cat "$output"/$i/$j/$k/$l/$m/$n/$o/$p/"$p"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$p '$6 ~ otu"-" {print $1}' > "$output"/$i/$j/$k/$l/$m/$n/$o/$p/"$p"_"$r"_whole_proteins_positives_codes.dat 
         csh evaluate_count_ranked.csh "$output"/$i/$j/$k/$l/$m/$n/$o/$p/"$p"_total.dat $output $r $p $i/$j/$k/$l/$m/$n/$o/
         cat "$output"/$i/$j/$k/$l/$m/$n/$o/$p/"$p"_class_"$r"_statistics.out | awk -v otu=$p '{print "EIGHTH",otu,$0}' >> "$output"/class_"$r"_statistics.dat
        end 
#--------------------------------------------------END OF EIGHTH OTU BLA EVAL----------------------------------------------------------------------------------------------------
        set ninth_tl=`cat "$results"_ranked.dat | awk '{print $5}' | sort -u | awk '{print $1}' | sed 's/-/ /g' | awk -v first=$i -v second=$j -v third=$k -v forth=$l -v fifth=$m -v sixth=$n -v seventh=$o -v eighth=$p '$1 == first && $2 == second && $3 == third && $4 == forth && $5 == fifth && $6 == sixth && $7 == seventh && $8 == eighth {print $9}' | sort -u`
        foreach q ( $ninth_tl )
         mkdir "$output"/$i/$j/$k/$l/$m/$n/$o/$p/$q
         cat "$output"/$i/$j/$k/$l/$m/$n/$o/$p/"$p"_total.dat | awk -v otu=$q '$5 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/$n/$o/$p/$q/"$q"_total.dat    
         foreach r ( $blas )
          cat "$output"/$i/$j/$k/$l/$m/$n/$o/$p/"$p"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$q '$6 ~ otu"-"' > "$output"/$i/$j/$k/$l/$m/$n/$o/$p/$q/"$q"_"$r"_whole_proteins_positives_meta_ranked.dat
          cat "$output"/$i/$j/$k/$l/$m/$n/$o/$p/$q/"$q"_"$r"_whole_proteins_positives_meta_ranked.dat | awk -v otu=$q '$6 ~ otu"-" {print $1}' > "$output"/$i/$j/$k/$l/$m/$n/$o/$p/$q/"$q"_"$r"_whole_proteins_positives_codes.dat 
          csh evaluate_count_ranked.csh "$output"/$i/$j/$k/$l/$m/$n/$o/$p/$q/"$q"_total.dat $output $r $q $i/$j/$k/$l/$m/$n/$o/$p/
          cat "$output"/$i/$j/$k/$l/$m/$n/$o/$p/$q/"$q"_class_"$r"_statistics.out | awk -v otu=$q '{print "NINTH",otu,$0}' >> "$output"/class_"$r"_statistics.dat
         end 
#--------------------------------------------------END OF NINTH OTU BLA EVAL----------------------------------------------------------------------------------------------------
        rm "$output"/$i/$j/$k/$l/$m/$n/$o/$p/$q/*_whole_proteins_positives_meta_ranked.dat
        end
#--------------------------------------------------END OF NINTH OTU ----------------------------------------------------------------------------------------------------
       rm "$output"/$i/$j/$k/$l/$m/$n/$o/$p/*_whole_proteins_positives_meta_ranked.dat
       end
#--------------------------------------------------END OF EIGHTH OTU----------------------------------------------------------------------------------------------------
      rm "$output"/$i/$j/$k/$l/$m/$n/$o/*_whole_proteins_positives_meta_ranked.dat
      end
#--------------------------------------------------END OF SEVENTH OTU----------------------------------------------------------------------------------------------------
     rm "$output"/$i/$j/$k/$l/$m/$n/*_whole_proteins_positives_meta_ranked.dat
     end
#--------------------------------------------------END OF SIXTH OTU----------------------------------------------------------------------------------------------------
    rm "$output"/$i/$j/$k/$l/$m/*_whole_proteins_positives_meta_ranked.dat
    end
#--------------------------------------------------END OF FIFTH OTU----------------------------------------------------------------------------------------------------
   rm "$output"/$i/$j/$k/$l/*_whole_proteins_positives_meta_ranked.dat
   end
#--------------------------------------------------END OF FORTH OTU----------------------------------------------------------------------------------------------------
  rm "$output"/$i/$j/$k/*_whole_proteins_positives_meta_ranked.dat
  end
#--------------------------------------------------END OF THIRD OTU----------------------------------------------------------------------------------------------------
 rm "$output"/$i/$j/*_whole_proteins_positives_meta_ranked.dat
 end
#--------------------------------------------------END OF SECOND OTU---------------------------------------------------------------------------------------------------
 rm "$output"/$i/*_whole_proteins_positives_meta_ranked.dat
end
#---------------------------------------------------END OF FIRST OTU-------------------------------------------------------------------------------------------------

        
       
           

