# script for getting prokaryotes lineage data

#1 download taxonomy dump files

mkdir temp

wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.Z --directory-prefix=temp

gunzip temp/*.Z

tar -xvf temp/*.tar -C temp/

cp temp/fullnamelineage.dmp fullnamelineage.dmp

rm -rf temp

cat fullnamelineage.dmp | awk -F "|" '$3 ~ "Bacteria" || $3 ~ "Archaea"' > prokaryotic_lineages.dat





