mkdir -p test/output
# Building from GenBank file and blast results
# ruby cgview_json_builder.rb -s test/input/sequence -o test/output/cgview.json -c test/input/config.yaml -b test/input/blast_results.txt

# Bulding from Genbank Contigs
ruby cgview_json_builder.rb -s test/input/chloroflexi_contigs.gbk -o test/output/cgview_contigs.json -c test/input/config.yaml -t test/input/contigs.csv




# BioRuby
# require 'bio'
# path = 'chloroflexi_contigs.gbk'
# flatfile = Bio::FlatFile.auto(path)
# so = flatfile.first
# bs = so.to_biosequence
# bs.na
# f1 = bs.features[2]
# l1 = Bio::Locations.new(f1.position).first
# bs.complement!
# f1c = bs.features[2]
# l1c = Bio::Locations.new(f1c.position).first

