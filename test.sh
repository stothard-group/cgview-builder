mkdir -p test/output
ruby cgview_json_builder.rb -s test/input/sequence -o test/output/cgview.json -c test/input/config.yaml -b test/input/blast_results.txt
