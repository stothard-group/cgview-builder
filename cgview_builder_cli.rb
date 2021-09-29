require_relative 'cgview_builder'
require 'optparse'
require 'ostruct'

# Version should match latest Major.Minor version of CGView.js
VERSION='1.1.0'

# Command line options will be stored in *options*
options = OpenStruct.new

# Define and grab the options 
optparse = OptionParser.new do |opts|
  opts.summary_width = 30

  opts.banner = "Usage: cgview_builder_cli.rb [options]"
  opts.separator ""
  opts.separator "Required Arguments:"

  opts.on("-s", "--sequence FILE", "Sequence File [GenBank for now]") do |seqfile|
    options.seqfile = seqfile
  end

  opts.on("-m", "--map_id STRING", "ID for map [Default: Random 40 character hex string]") do |map_id|
    options.map_id = map_id
  end

  opts.on("-n", "--map_name STRING", "Name for map") do |map_name|
    options.map_name = map_name
  end

  opts.on("-o", "--outfile FILE", "Write JSON to this file") do |outfile|
    options.outfile = outfile
  end

  opts.separator ""
  opts.separator "General Options:"

  opts.on("-c", "--config FILE", "Config File") do |config|
    options.config = config
  end

  opts.on("-b", "--blasts FILEs", "One or more blast files (separated by ,)") do |blast_paths|
    options.blast_paths = blast_paths
  end

  # This will print an options summary.
  opts.on('-v', '', '--version', "Print version") do
    puts "CGViewBuilder #{VERSION}"
    exit!
  end

  # This will print an options summary.
  opts.on('-h', '-?', '--help', "Show this message") do
    puts opts
    exit!
  end
  opts.separator ""
end

# Parse command line arguments
begin
  optparse.parse!(ARGV)
rescue Exception => e
  puts e, "", optparse
  exit
end

# Check for required arguments
if !(options.seqfile && options.outfile) then
  puts "\nMissing Required Arguments!", "", optparse
  exit
end

cgview_options = {
  config: options.config,
  map_id: options.map_id,
  map_name: options.map_name,
  contigs: options.contigs,
}

if options.blast_paths
  cgview_options[:blasts] =  options.blast_paths.split(',')
end

cgview = CGViewBuilder.new(options.seqfile, cgview_options)
cgview.write_json(options.outfile)








