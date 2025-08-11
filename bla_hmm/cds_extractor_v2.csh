#cds extractor v 2.0 for bla_hmm
set contig=$1 
set sequen_g=$2
set frag_start=$3
set frag_end=$4
set strand=$5
set sequen=$6
set protein_frag=$7
set annotation=$8
set tolerance=$9 # limit of nucleotides difference between hit and identified ORF in either direction in case of truncated ORF
set output=$10
set outdir=$11

rm -rf $outdir/cds_extraction
mkdir $outdir/cds_extraction

#check if genome is annotated and if not annotate the contig with PROKKA

set annot=`ls $outdir | awk '$0 ~ "annotation"' | wc -l`
if ($annot == 0) then
 echo "GENOME NOT ANNOTATED...PROKKA ANNOTATION OF GENOME HAS STARTED"
 prokka $sequen_g --outdir $outdir/annotation --addgenes --norrna --notrna --force --cpus 8 --quiet
 set gff=$outdir/annotation/*.gff
 set faa=$outdir/annotation/*.faa
else
 set gff=$outdir/annotation/*.gff
 set faa=$outdir/annotation/*.faa
endif

cat $faa > $outdir/cds_extraction/proteins.faa
echo ">" >> $outdir/cds_extraction/proteins.faa

echo $protein_frag | sed 's/_//g' | awk '{print ">""\n"$1}' > $outdir/cds_extraction/protein_fragment.fasta # _ deletion required in case of stop codon in hit

if ($strand == "plus") then
 set direction="+"
else
 set direction="-"
endif

cat $gff | sed '/#/d' | awk -F'\t' -v contig=$contig '$1 == contig' | awk -F'\t' '$3 == "CDS"' | awk -v direc=$direction -F'\t' '$7 == direc' > $outdir/cds_extraction/gff

cat $outdir/cds_extraction/gff | awk -v b=$frag_start -v e=$frag_end -F'\t' '$4 <= e && $5 >= b' | awk -v b=$frag_start -v e=$frag_end -v tol=$tolerance -F'\t' '$4-b <= tol && e-$5 <= tol' > $outdir/cds_extraction/target_annot.temp

set coord_targets=`cat $outdir/cds_extraction/target_annot.temp | wc -l`

if ($coord_targets != "0") then
 set protein_ids_c=`cat $outdir/cds_extraction/target_annot.temp | awk -F'\t' '$0 !~ "pseudo=true" {print $9}' | wc -l`
 if ($protein_ids_c == "0") then 
  set orf_check="pseudogene"
  set protein_id="PSEUDO"
  set cds_start=$frag_start
  set cds_end=$frag_end
  set protein_annot="annotated_as_pseudogene"
  set protein_seq=$protein_frag
 else
  set protein_ids=`cat $outdir/cds_extraction/target_annot.temp | awk -F'\t' '$0 !~ "pseudo=true" {print $9}' | sed 's/\;/ /g' | awk '{print $1}' | sed 's/ID=cds-//g' | sed 's/ID=//g' | sed 's/-/ /g' | awk '{print $1}'`
  foreach i ($protein_ids)
   cat $outdir/cds_extraction/proteins.faa | sed -n '/>'$i'/,/>/p' | sed '$d' >> $outdir/cds_extraction/target_proteins.fasta
  end
  #blast database
  makeblastdb -in $outdir/cds_extraction/target_proteins.fasta -input_type fasta -dbtype prot -title "proteins" -out  $outdir/cds_extraction/target_proteins_db >&-
  rm -
  blastp -evalue 0.001 -query $outdir/cds_extraction/protein_fragment.fasta -task blastp -max_hsps 1 -culling_limit 1 -db $outdir/cds_extraction/target_proteins_db -outfmt "7 qacc sallseqid evalue bitscore pident qstart qend sstart send mismatch gaps qframe qseq" -out $outdir/cds_extraction/blastp_results.dat -soft_masking true
  set best_hit=`cat $outdir/cds_extraction/blastp_results.dat | sed '/#/d' | sort -k4,4nr | head -n1 | awk '{print $2}'`
  if ($best_hit != "") then
   if ($annotation == "no") then
    set protein_id="PROKKA"
   else
    set protein_id="$best_hit"
   endif
   set orf_check="ORF_OK"
   set cds_start=`cat $outdir/cds_extraction/target_annot.temp | awk -v best=$best_hit -F'\t' '$9 ~ best {print $4}'`
   set cds_end=`cat $outdir/cds_extraction/target_annot.temp | awk -v best=$best_hit -F'\t' '$9 ~ best {print $5}'`
   set protein_annot=`cat $outdir/cds_extraction/target_annot.temp | awk -v best=$best_hit -F'\t' '$9 ~ best {print $9}' | grep -oP '(?<=product=).*' | awk -F';' '{print $1}' | sed 's/ /_/g'`
   set partial_check=`cat $outdir/cds_extraction/target_annot.temp | awk -v best=$best_hit -F'\t' '$9 ~ best && $9 ~ "partial=true" {print "_partial"}'`
   set protein_annot=`echo "$protein_annot""$partial_check"`
   set protein_seq=`cat $outdir/cds_extraction/proteins.faa | sed -n '/>'$best_hit'/,/>/p' | sed '/>/d' | tr -d '\n'`
  else
   set orf_check="pseudogene"
   set protein_id="YYYYYY"
   set cds_start=$frag_start
   set cds_end=$frag_end
   set protein_annot="not_annotated_partial"
   set protein_seq=$protein_frag
  endif
 endif
else
 set orf_check="pseudogene"
 set protein_id="YYYYYY"
 set cds_start=$frag_start
 set cds_end=$frag_end
 set protein_annot="not_annotated_partial"
 set protein_seq=$protein_frag
endif

if ($strand == "plus") then
 set cds_seq=`cat $sequen | tr -d '\n' | cut -c$cds_start-$cds_end`
else 
 set cds_seq=`cat $sequen | tr -d '\n' | cut -c$cds_start-$cds_end | rev | tr '[CGAT]' '[GCTA]'`
endif 

echo $cds_start $cds_end $cds_seq $protein_id $protein_annot $protein_seq $orf_check > $outdir/cds_extraction/"$output"_wseqs.dat




