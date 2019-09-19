require 'json'
require 'yaml'
require 'bio'
require 'csv'
require 'ostruct'
require 'securerandom'


class CGViewJSON

  VERSION = '0.1'

  attr_accessor :config, :options, :seq_object, :sequence, :cgview, :seq_type, :features,
                :tracks, :debug, :captions, :contigs

  def initialize(sequence_path, options={})
    @map_id = options[:map_id] || SecureRandom.hex(20)
    @map_name = options[:map_name]
    @cgview = initialize_cgview
    @options = options
    @features = []
    @plots = []
    @tracks = []
    @blast_tracks = []
    @debug = options[:debug]
    @config = options[:config] ? read_config(options[:config]) : {}
    @contigs = options[:contigs] ? read_contigs(options[:contigs]) : []
    read_sequence(sequence_path)
    # extract_features
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
      created: Time.now.strftime("%Y-%m-%d %H:%M:%S"),
      id: @map_id,
      name: @map_name,
      settings: {},
      backbone: {},
      ruler: {},
      dividers: {},
      sequence: {},
      captions: [],
      legend: {},
      features: [],
      tracks: []
    }

  end

  def read_config(path)
    config = symbolize(YAML.load_file(path)['cgview'])
    @cgview[:settings] = config[:settings] if config[:settings]
    @cgview[:backbone] = config[:backbone] if config[:backbone]
    @cgview[:ruler] = config[:ruler] if config[:ruler]
    @cgview[:dividers] = config[:dividers] if config[:dividers]
    @cgview[:sequence] = config[:sequence]
    @cgview[:legend] = config[:legend]
    @cgview[:tracks] = config[:tracks] || []
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
    # print "Reading sequence file..."
    puts "Extracting sequence and features..."
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
      print "\nCould not autodetect filetype. Trying different method..."
      @seq_type = detect_filetype(path)
      if @seq_type != :raw
        case @seq_type
        when :genbank
          flatfile = Bio::FlatFile.open(Bio::GenBank, path)
        when :embl
          flatfile = Bio::FlatFile.open(Bio::EMBL, path)
        when :fasta
          flatfile = Bio::FlatFile.open(Bio::FastaFormat, path)
        end
      end
      # @seq_type = :raw
    end

    # Extract sequence
    sequence_num = 0
    sequence_length = 0
    if @seq_type == :raw
      @sequence = File.read(path).gsub(/[^A-Za-z]/, '').upcase
      sequence_num += 1
      sequence_length = @sequence.length
    else
      if self.contigs?
        # puts @contigs
        print "Contigs:"
        flatfile.each_with_index do |seq_object, i|
          sequence_num += 1
          print "."
          # find contig for this seq_object using seq_object.entry_id
          # add sequence to contig object
          #   - if it's reversed, then complement the sequence
          id = seq_object.entry_id
          contig = @contigs.find { |c| c[:name] == id }
          # TODO: complement if needed
          biosequence = seq_object.to_biosequence
          if contig[:orientation] == '-'
            biosequence.na
            complement_biosequence = biosequence.complement
            contig[:seq] = complement_biosequence.to_s.upcase
          else
            contig[:seq] = biosequence.to_s.upcase
          end
          # puts " #{i}:#{contig[:seq][0..10]}"
          # puts contig
          unless contig
            fail "Cound not find contig '#{id}'; sequence index (base-1): #{i+1}"
          end
          unless contig[:length] == seq_object.length
            fail "Contig '#{id}' length '#{contig[:length]}' != sequence length '#{seq_object.length}'"
          end
          sequence_length += seq_object.length

          extract_features(seq_object, contig)
          if i == 0
            @map_name ||= seq_object.entry_id
          end
        end
        puts ""
      else
        seq_object = flatfile.first
        @sequence = seq_object.to_biosequence.to_s.upcase
        extract_features(seq_object)
        sequence_num += 1
        sequence_length = @sequence.length
        @map_name ||= seq_object.entry_id
      end
    end
    # puts "Count: #{sequence_num}, Length: #{sequence_length} bp, Type: #{@seq_type}"
    puts "Extracted Sequences: #{sequence_num}, Length: #{sequence_length} bp, Type: #{@seq_type}"
    puts "Extracted Features: #{@features.count}"
  end

  # Simple method for detecting file type
  # Mirrors JavaScript method in CGView Server
  def detect_filetype(path)
    first_line = ''
    File.foreach(path) do |line|
      next if line.empty?
      next if line =~ /^#/
      first_line = line
      break
    end
    case first_line
    when /^LOCUS\s+/ then :genbank
    when /^ID\s+/    then :embl
    when /^>/     then :fasta
    else               :raw
    end
  end

  def extract_features(seq_object, contig=nil)
    # if contig
    #   puts(contig[:id], contig[:name], contig[:orientation])
    # end
    return unless [:embl, :genbank].include?(@seq_type)
    # print "Extracting features..."
    features_to_skip = ['source', 'gene', 'exon']
    # TODO: look into complex features from xml-builder
    seq_object.features.each do |feature|
      next if features_to_skip.include?(feature.feature)
      next if feature.position.nil?
      locations = Bio::Locations.new(feature.position)
      unless locations.first == locations.last
        # FIXME Complex Feature...What Now
        # - Each location should become a feature
        # - The features could have compoundfeatures attribute to link to join features
        next
      end
      # Feature Name
      # NOTE: This converts the array of qualifiers to a easily accessible hash.
      # However, there is a risk some information is lost when two or more qualifiers are the same.
      qualifiers = feature.assoc
      name = qualifiers['gene'] || qualifiers['locus_tag'] || qualifiers['note'] || feature.feature
      # FIXME: addtional qualifiers could become part of meta tag called 'qualifiers'
      # Feature Location
      location = locations.first
      # Skip features with the same length as the sequence
      next if location.from == 1 && location.to == seq_object.seq.length
      start = location.from
      stop = location.to
      strand = location.strand

      # Handle contig
      if contig && contig[:orientation] == '-'
        strand = (strand == 1) ? -1 : 1
        stop = contig[:length] - location.from + 1
        start = contig[:length] - location.to + 1
      end

      # Create Feature
      feature = {
        type: feature.feature,
        name: name,
        start: start,
        stop: stop,
        strand: strand,
        # FIXME: source becomes collection
        # FIXME: more collections: genbank-cds, genbank-trna, ensembl-cds, etc
        source: "sequence-features"
      }
      if contig
        feature[:contig] = contig[:id]
      end
      @features.push(feature)
    end
    # puts @features.count
  end

  def build_legend
    print "Building legend..."
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
      else
        feature[:legend] = feature[:type]
      end
    end
    # Intersection of legend names (They will be in the same order as the config
    feature_legend_names = config_items.keys & @features.map { |f| f[:legend] }.uniq
    items = []
    feature_legend_names.each { |n| items.push config_items[n] }
    @cgview[:legend][:items] = items
    puts feature_legend_names.count
  end

  def build_captions
    print "Building captions..."
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
    puts @captions.count
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
          identity: row[2].to_f,
          mimatches: row[4].to_i,
          evalue: row[10].to_f,
          score: row[11].to_i
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
          score: (meta[:identity] / 100).round(3),
          # FIXME: source becomes collection
          source: "blast_#{num}"
        })
          puts (meta[:identity] / 100).round(3)
      end

      # Create Track
      @blast_tracks << {
        name: "blast_#{num}",
        position: 'inside',
        strand: 'combined',
        contents: {
          type: 'feature',
          # FIXME: source becomes collection
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

  # Config file format:
  # id,index,original_index,length,strand
  def read_contigs(path)
    puts "Reading contigs..."
    contigs = []
    CSV.foreach(path, headers: true) do |row|
      length = row['length'].to_i
      contigs.push({
        id: "cgv-contig-#{row['original_index']}",
        name: row['id'],
        length: length,
        orientation: row['strand']
      })
    end
    contigs
  end

  def contigs?
    @contigs && @contigs.length > 0
  end

  def read_contigs_OLD(path)
    puts "Creating contig features..."
    # Create Features
    start = 1
    CSV.foreach(path, headers: true) do |row|
      length = row['length'].to_i
      @features.push({
        type: 'Contig',
        name: row['id'],
        # meta: meta,
        start: start,
        stop: start + length - 1,
        strand: row['strand'],
        # FIXME: source becomes collection
        source: "contigs",
      })
      start += length.to_i
    end

    # Create Track
    @contig_track = {
      name: "Contigs",
      position: 'inside',
      strand: 'combined',
      thicknessRatio: 0.5,
      contents: {
        type: 'feature',
        # FIXME: source becomes collection
        from: 'source',
        extract: "contigs"
      }
    }
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
    print "Building Tracks..."
    @tracks << {
      name: 'Features',
      readingFrame: 'combined',
      strand: 'separated',
      position: 'both',
      contents: {
        type: 'feature',
        # FIXME: source becomes collection
        from: 'source',
        extract: 'sequence-features'
      }
    }

    # if @contig_track
    #   @tracks << @contig_track
    # end

    unless @blast_tracks.empty?
      @tracks += @blast_tracks
    end

    if @options[:analysis_path]
      @tracks << {
        name: 'Analysis',
        position: 'inside',
        contents: {
          type: 'plot',
          # FIXME: source becomes collection
          from: 'source',
          extract: 'analysis'
        }
      }
    end
    puts @tracks.count
  end

  def build_cgview
    puts "Creating CGView JSON"
    if @debug
      @cgview[:sequence][:seq] = "SEQUENCE WOULD GO HERE"
      @cgview[:features] += @features[1..5]
    else
      if contigs?
        @cgview[:sequence][:contigs] = @contigs
      else
        # FIXME: this will become "chromosome"
        @cgview[:sequence][:seq] = @sequence
      end
      @cgview[:features] += @features
    end
    @cgview[:tracks] = @tracks + @cgview[:tracks]
    unless @plots.empty?
      @cgview[:plots] = @plots
    end

    @cgview[:captions] += @captions
  end


  def map_title
    if @map_name
      @map_name
    elsif @seq_type == :raw
      ''
    else
      # FIXME: there will be no more seq_object
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

# class Contig
#   attr_accessor :id, :name, :orientation, :seq, :length
# end

def fail(msg)
  puts msg
  exit
end

# debug = false
# # # debug = true
# # file = "data/sequences/NC_001823.gbk" # 70 KB
# file = "data/sequences/NC_000907.gbk" # 1.8 MB
# # # file = "data/sequences/NC_000913.gbk" # 4.6 MB
# config_path = 'scripts/cgview-builder/config_example.yaml'
# cgview = CGViewJSON.new(file, config: config_path, debug: debug)
# #
# # # file = "data/sequences/B_pert_TahomaI.gbk" # 4 MB
# # # config_path = 'scripts/cgview_json_builder/test_config.yaml'
# # # analysis_path = '/Users/jason/Desktop/merged_hits_cov.gff'
# # # cgview = CGViewJSON.new(file, config: config_path, debug: debug, analysis_path: analysis_path)
# #
# #
# cgview.write_json("/Users/jason/workspace/stothard_group/cgview-js/data/tests/builder.json")

