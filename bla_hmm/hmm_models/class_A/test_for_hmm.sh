
export protein_ids=`cat cdhit_clustered_class_A | awk '$1 ~ ">" {print $1}' | sed 's/>//g' | sed 's/ /\n/g'`
sed -i "\$a>" cdhit_clustered_class_A
for i in $protein_ids
 do
 mkdir $i
 cat cdhit_clustered_class_A | sed '/^>'$i'/,/^>/{/^>'$i'/!{/^>/!d}}' | sed '/>'$i'/d' > "$i"/class_A_w_o_cluster_"$i".fasta
 cat cdhit_clustered_class_A | sed -n '/>'$i'/,/>/p' | sed '$ d' > "$i"/"$i".fasta
 sed -i '$ d' "$i"/class_A_w_o_cluster_"$i".fasta
 clustalo -i "$i"/class_A_w_o_cluster_"$i".fasta -o "$i"/class_A_w_o_cluster_"$i"_aligned.fasta
 hmmbuild "$i"/class_A_w_o_cluster_"$i".hmm "$i"/class_A_w_o_cluster_"$i"_aligned.fasta
 hmmsearch "$i"/class_A_w_o_cluster_"$i".hmm "$i"/"$i".fasta > "$i"/"$i"_vs_model.out
 export score=`cat "$i"/"$i"_vs_model.out | sed -n 15p | awk '{print $2}'`
 echo $i $score >> scores.dat
done
sed -i '$ d' cdhit_clustered_class_A
clustalo -i cdhit_clustered_class_A -o cdhit_clustered_class_A_aligned.fasta
hmmbuild class_A_model.hmm cdhit_clustered_class_A_aligned.fasta


