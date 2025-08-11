# wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/archaea/assembly_summary.txt gia genbank. gia refseq antikatastash
#  cat archaea_genbank_assembly_summary.txt | awk -F'\t' '{print $1}' > archaea_genbank_assemblies.dat


setenv PATH $PATH\:$HOME/bla_hmm

set ass=$1

set domain=$2

set assembly_summary=$3

set curr_date=`date +%D | sed 's/\//_/g'`

set bacteria=`cat $domain/"$curr_date"_assemblies/ASS"$ass" | awk '{print $1}' | sort -u`


foreach j ($bacteria)
 mkdir $domain/$j
 #retrieve taxonomy
 set taxid=`cat $domain/$assembly_summary | awk -F'\t' -v as=$j '$1 == as {print $6}'`
 set spec_taxid=`cat $domain/$assembly_summary | awk -F'\t' -v as=$j '$1 == as {print $7}'`
 set species=`cat $domain/$assembly_summary | awk -F'\t' -v as=$j '$1 == as {print $8}' | sed 's/ /_/g' | sort -u | sed 's/\[/</g' | sed 's/\]/>/g' | sed 's/\*//g'`
 set lineage=`cat taxonomy/prokaryotic_lineages.dat | awk -F "|" -v tx=$spec_taxid '$1 == tx {print $3}' | sed 's/\t//g' | sed 's/\; /-/g' | sed 's/ /_/g' | sed 's/\[/</g' | sed 's/\]/>/g' | sed 's/\*//g'`

 if ($#taxid == 0) then
  set taxid="NO_TAXID"
 endif

 if ($#spec_taxid == 0) then
  set spec_taxid="NO_SPECIES_TAXID"
 endif

 if ($#species == 0) then
  set species="NO_SPECIES"
 endif

 if ($#lineage == 0) then
  set lineage="NO_LINEAGE"
 endif
 echo $j "$taxid" "$spec_taxid" "$species" "$lineage"

 #download fna, md5checksums, gff, faa
 set file_path_b=`cat $domain/$assembly_summary | awk -F'\t' -v ass=$j '$1 == ass {print $20}' | sed 's/https:/rsync:/g'`
 set preffix_b=`echo $file_path_b | grep -oP '(?='$j').*?(?=$)'`
 csh download_rsync_issue.csh $file_path_b/"$preffix_b"_genomic.fna.gz $file_path_b/md5checksums.txt $domain $j fna txt
 set annotation=`cat $domain/$j/md5checksums.txt | awk '$0 ~ "gff"' | wc -l`
  if ($annotation != "0") then
   csh download_rsync_issue.csh $file_path_b/"$preffix_b"_genomic.gff.gz $file_path_b/"$preffix_b"_protein.faa.gz $domain $j gff faa
  endif
 cd $domain/$j
 echo "y" | gunzip *.gz

 #bla_hmm analysis
 bla_hmm_v11.csh "$preffix_b"_genomic.fna all
 cat blas_hmm_results/positives/*_positives_total.dat | awk -v i=$j -v k=$taxid -v m=$spec_taxid -v l=$species -v o=$lineage '{print i,k,m,l,o,$0}' >> ../../$domain/"$curr_date"_bla_hmm_results/"$ass"_positives_total.tab
 cat blas_hmm_results/greys/*_greys_total.dat | awk -v i=$j -v k=$taxid -v m=$spec_taxid -v l=$species -v o=$lineage '{print i,k,m,l,o,$0}' >> ../../$domain/"$curr_date"_bla_hmm_results/"$ass"_greys_total.tab 
 cat blas_hmm_results/positives/*_pos_nums.dat | awk -v i=$j -v k=$taxid -v m=$spec_taxid -v l=$species -v o=$lineage '{print i,k,m,l,o,$0}' >> ../../$domain/"$curr_date"_bla_hmm_results/"$ass"_positives_stats.dat
 cat blas_hmm_results/greys/*_greys_nums.dat | awk -v i=$j -v k=$taxid -v m=$spec_taxid -v l=$species -v o=$lineage '{print i,k,m,l,o,$0}' >> ../../$domain/"$curr_date"_bla_hmm_results/"$ass"_greys_stats.dat
 cd ../..
 rm -rf $domain/$j
end  


