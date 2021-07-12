# CGView Builder

CGViewBuilder is a ruby script that generates JSON files for
[CGView.js](http://cgview.ca) from sequence files (e.g. GenBank, EMBL, FASTA).

CGViewBuilder consists of two files:

- cgview_builder.rb contains the CGViewBuilder class which can be used by other ruby programs.
- cgview_builder_cli.rb is a command line interface to the CGViewBuilder class.

## Installation

- Download the git repository.
- Install CGViewBuilder's one dependency [BioRuby](https://bioruby.org):
```bash
gem install bio
```
- Confirm installation by running the test script:
```bash
./test.sh
```


## Running the CLI script:

```bash
ruby cgview_builder_cli.rb [OPTIONS]
```

## Command line options

```bash
Usage: cgview_builder_cli.rb [options]

Required Arguments:
    -s, --sequence FILE            Sequence File [GenBank for now]
    -m, --map_id STRING            ID for map [Default: Random 40 character hex string]
    -n, --map_name STRING          Name for map
    -o, --outfile FILE             Write JSON to this file

General Options:
    -c, --config FILE              Config File
    -b, --blasts FILEs             One or more blast files (separated by ,)
    -v, --version
                                   Print version
    -h, -?, --help                 Show this message
```


