set bla=$1
set PFs=(PF13354 PF02113 PF00905 PF00768 PF00144)
set cut_off=`sort -k2 -n scores.dat | head -1 | awk '{print int($2)}'`

foreach i ($PFs)
 hmmsearch class_"$bla"_model.hmm ../PFAM_serine_blas_clan/"$i"_clustered > "$i"_vs_class_"$bla"_model.dat
 cat "$i"_vs_class_"$bla"_model.dat | sed -n '/Scores for complete sequences/,/Domain annotation for each/p' | sed '/------ inclusion threshold ------/d' | sed '1,4d' | head -n -3 > "$i"_vs_class_"$bla"_model.tab
 set seqs=`cat ../PFAM_serine_blas_clan/"$i"_clustered | awk '$1 ~ ">" {print $1}' | wc -l`
 set hits=`cat "$i"_vs_class_"$bla"_model.tab | awk '{print $1}' | wc -l`
 set hits_co=`cat "$i"_vs_class_"$bla"_model.tab | awk -v c=$cut_off '$2 > c {print $2}' | wc -l` 
 echo "$i"" vs class_"$bla"_model cut-off=""$cut_off"" ""$seqs"" / ""$hits"" / ""$hits_co" >> class_"$bla"_model_vs_PFs.dat
end

