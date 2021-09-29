DIR=`dirname $0`
rm -rf "${DIR}/output"
mkdir -p "${DIR}/output"

# Building from GenBank file and blast results
# ruby cgview_json_builder.rb -s test/input/sequence -o test/output/cgview.json -c test/input/config.yaml -b test/input/blast_results.txt

# Building from Genbank Contigs with no optional contig file
# ruby cgview_json_builder.rb -s test/input/chloroflexi_contigs.gbk -o test/output/cgview_contigs_no_file.json -c test/input/config.yaml
ruby cgview_builder_cli.rb \
  -s "${DIR}/input/chloroflexi_contigs.gbk" \
  -o "${DIR}/output/cgview_contigs.json" \
  -c "${DIR}/input/config.yaml"


