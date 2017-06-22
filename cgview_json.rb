require 'json'
require 'yaml'
require 'bio'
require 'csv'
require 'ostruct'


class CGViewJSON

  VERSION = '0.1'

  attr_accessor :config, :options, :seq_object, :sequence, :cgview, :seq_type, :features,
                :tracks, :debug, :captions

  def initialize(sequence_path, options={})
    @cgview = initialize_cgview
    @options = options
    @features = []
    @plots = []
    @tracks = []
    @blast_tracks = []
    @debug = options[:debug]
    @config = options[:config] ? read_config(options[:config]) : {}
    read_sequence(sequence_path)
    extract_features
    read_gff_analysis(options[:analysis_path]) if options[:analysis_path]
    read_blasts(options[:blasts]) if options[:blasts]
    build_legend
    build_captions
    build_tracks
    build_cgview
  end

  def initialize_cgview
    {
      version: VERSION,
      created: Time.now.strftime("%Y-%m-%d %H:%M:%S")
      settings: {},
      sequence: {},
      captions: [],
      legend: {},
      features: [],
      layout: { tracks: [] }
    }

  end

  def read_config(path)
    config = symbolize(YAML.load_file(path)['cgview'])
    @cgview[:settings] = config[:settings]
    @cgview[:sequence] = config[:sequence]
    @cgview[:legend] = config[:legend]
    @cgview[:layout][:tracks] = config[:layout] && config[:layout][:tracks] || []
    config
  end

  def symbolize(obj)
    if obj.is_a? Hash
      return obj.inject({}) do |memo, (k, v)|
        memo.tap { |m| m[k.to_sym] = symbolize(v) }
      end
    elsif obj.is_a? Array
      return obj.map { |memo| symbolize(memo) }
    end
    obj
  end

  def read_sequence(path)
    puts "Reading sequence file..."
    flatfile = Bio::FlatFile.auto(path)
    # Determine Sequence file type
    case flatfile.dbclass.to_s
    when 'Bio::GenBank'
      @seq_type = :genbank
    when 'Bio::EMBL'
      @seq_type = :embl
    when 'Bio::FastaFormat'
      @seq_type = :fasta
    else
      @seq_type = :raw
    end

    # Extract sequence
    if @seq_type == :raw
      @sequence = File.read(path).gsub(/[^A-Za-z]/, '').upcase
    else
      @seq_object = flatfile.first
      @sequence = @seq_object.to_biosequence.to_s.upcase
    end
  end

  def extract_features
    puts "Extracting features..."
    return unless [:embl, :genbank].include?(@seq_type)
    features_to_skip = ['source', 'gene', 'exon']
    # TODO: look into complex features from xml-builder
    @seq_object.features.each do |feature|
      next if features_to_skip.include?(feature.feature)
      next if feature.position.nil?
      locations = Bio::Locations.new(feature.position)
      unless locations.first == locations.last
        # FIXME Complex Feature...What Now
        next
      end
      # Feature Name
      # NOTE: This converts the array of qualifiers to a easily accessible hash.
      # However, there is a risk some information is lost when two or more qualifiers are the same.
      qualifiers = feature.assoc
      name = qualifiers['gene'] || qualifiers['locus_tag'] || qualifiers['note'] || feature.feature
      # Feature Location
      location = locations.first
      # Skip features with the same length as the sequence
      next if location.from == 1 && location.to == @seq_object.length
      # Create Feature
      @features.push({
        type: feature.feature,
        name: name,
        start: location.from,
        stop: location.to,
        strand: location.strand,
        source: "sequence-features"
      })
    end
  end

  def build_legend
    puts "Building legend..."
    config_items = {}
    default_legend_name = nil
    # Read config file legend items
    if @config[:legend] && @config[:legend][:items].is_a?(Array)
      @config[:legend][:items].each { |i| config_items[i[:name]] = i }
      default_legend_name =  @config[:legend][:default]
    end
    # FIXME: add default legend if one does not exist
    @features.each do |feature|
      if config_items[feature[:type]]
        feature[:legend] = config_items[feature[:type]][:name]
      elsif default_legend_name
        feature[:legend] = config_items[default_legend_name][:name]
      end
    end
    # Intersection of legend names (They will be in the same order as the config
    feature_legend_names = config_items.keys & @features.map { |f| f[:legend] }.uniq
    items = []
    feature_legend_names.each { |n| items.push config_items[n] }
    @cgview[:legend][:items] = items
  end

  def build_captions
    puts "Building captions..."
    @captions = []
    config_captions = @config[:captions].is_a?(Array) ? @config[:captions] : []
    config_captions.each do |caption|
      if caption[:items]
        @captions << caption
      elsif caption[:name].downcase == 'title' && map_title != ""
        caption[:items] = [ { name:  map_title}]
        @captions << caption
      end
    end
  end

  # Currently only reads blastn, blastx, tblastx results properly
  # Blast results are expected to have the typical format without a header
  # (i.e. option -outfmt 6)
  # Columns
  # 0: query_id
  # 1: match_id
  # 2: %_identity
  # 3: alignment_length
  # 4: mismatches
  # 5: gap_openings
  # 6: q_start
  # 7: q_end
  # 8: s_start
  # 9: s_end
  # 10: evalue
  # 11: bit_score
  # OPTIONAL
  #  - Custom formatting can be used to provide the strand information
  #    12: strand (optional) NOT IMPLEMENTED YET
  #  - The query_id may contain additional details such as the start, end of the query in
  #    relation to the genome sequence
  #    This format is produced from the sequence_to_multi_fasta.pl script
  #    (e.g. some_id._start=123;end=456)
  def read_blasts(paths)
    paths.each_with_index do |path, i|
      puts "Creating features for BLAST #{i+1} results..."
      num = i + 1
      # Create Features
      CSV.foreach(path, col_sep: "\t") do |row|
        query_id = row[0]
        match_id = row[1]
        start = row[6].to_i
        stop = row[7].to_i
        strand = 1
        offset = 0
        if query_id =~ /^([^\t]+)_start=(\d+);end=(\d+)/
          offset = $2.to_i - 1
        end

        # Collect meta data
        meta = {
          identity: row[2],
          mimatches: row[4],
          evalue: row[10],
          score: row[11]
        }

        if start > stop
          start,stop = stop,start
          strand = -1
        end

        @features.push({
          type: 'blast',
          meta: meta,
          start: offset + start,
          stop: offset + stop,
          strand: strand,
          source: "blast_#{num}"
        })
      end

      # Create Track
      @blast_tracks << {
        name: "blast_#{num}",
        position: 'inside',
        strand: 'combined',
        contents: {
          type: 'feature',
          from: 'source',
          extract: "blast_#{num}"
        }
      }
    end
  end

  def parse_query_id(id)
    query = OpenStruct.new
    if id =~ /^([^\t]+)_start=(\d+);end=(\d+)/
      query.start = $2
      query.end = $3
    else
      query.start = 1
    end
    query
  end

  def read_gff_analysis(path)
    starts = []
    stops = []
    raw_scores = []
    CSV.foreach(path, col_sep: "\t", headers: true) do |row|
      starts << row['start'].to_i
      stops << row['end'].to_i
      raw_scores << row['score'].to_f
    end
    max = raw_scores.max.to_f
    min = raw_scores.min.to_f
    baseline = 0
    if min < 0
      baseline = scale_score(0, min, max)
    end
    positions = []
    scores = []
    starts.each_with_index do |start, i|
      stop = stops[i]
      score = scale_score(raw_scores[i], min, max)
      positions << start
      scores << score
      positions << stop
      scores << baseline
    end
    puts(min)
    puts(max)
    puts(raw_scores[1..10])
    puts(scores[1..10])

    @plots.push({
      source: 'analysis',
      positions: positions,
      scores: scores,
      baseline: baseline
    })

  end

  def scale_score(score, min, max)
    (score - min) / (max - min)
  end

  def build_tracks
    puts "Building Tracks..."
    @tracks << {
      name: 'Features',
      readingFrame: 'combined',
      strand: 'separated',
      position: 'both',
      contents: {
        type: 'feature',
        from: 'source',
        extract: 'sequence-features'
      }
    }

    unless @blast_tracks.empty?
      @tracks += @blast_tracks
    end

    if @options[:analysis_path]
      @tracks << {
        name: 'Analysis',
        position: 'inside',
        contents: {
          type: 'plot',
          from: 'source',
          extract: 'analysis'
        }
      }
    end
  end

  def build_cgview
    puts "Creating CGView JSON"
    if @debug
      @cgview[:sequence][:seq] = "SEQUENCE WOULD GO HERE"
      @cgview[:features] += @features[1..5]
    else
      @cgview[:sequence][:seq] = @sequence
      @cgview[:features] += @features
    end
    # @cgview[:layout][:tracks] += @tracks
    @cgview[:layout][:tracks] = @tracks + @cgview[:layout][:tracks]
    unless @plots.empty?
      @cgview[:plots] = @plots
    end

    @cgview[:captions] += @captions
  end


  def map_title
    if @options[:mapTitle]
      @options[:mapTitle]
    elsif @seq_type == :raw
      ''
    else
      @seq_object.definition
    end
  end

  def to_json
    JSON.generate({ cgview: @cgview })
  end

  def write_json(path)
    File.open(path, 'w') { |f| f.write(self.to_json) }
  end

end

debug = false
# # debug = true
# file = "data/sequences/NC_001823.gbk" # 70 KB
file = "data/sequences/NC_000907.gbk" # 1.8 MB
# # file = "data/sequences/NC_000913.gbk" # 4.6 MB
config_path = 'scripts/cgview-builder/config_example.yaml'
cgview = CGViewJSON.new(file, config: config_path, debug: debug)
#
# # file = "data/sequences/B_pert_TahomaI.gbk" # 4 MB
# # config_path = 'scripts/cgview_json_builder/test_config.yaml'
# # analysis_path = '/Users/jason/Desktop/merged_hits_cov.gff'
# # cgview = CGViewJSON.new(file, config: config_path, debug: debug, analysis_path: analysis_path)
#
#
cgview.write_json("/Users/jason/workspace/stothard_group/cgview-js/data/tests/builder.json")

