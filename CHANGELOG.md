--------------------------------------------------------------------------------
# CGViewBuilder Changelog
--------------------------------------------------------------------------------

## Unreleased/Tagged
- Added product and db_xref as fallback for feature names
- Remove semilcolons from contigs names (replace with _)
- Fix error where empty qualifiers strings were converted to Boolean
- Fix contig names encoding issues ("\\xBB" from ASCII-8BIT to UTF-8 (Encoding::UndefinedConversionError))

## v1.1.1 - 2022-10-27
- Fixed: raw sequences now work

## v1.1.0 - 2021-09-29

- Remove 'id' attribute from Contig (Use 'name' instead)
- Duplicate Contig names will automically become unique by appending a number
- Default Annotation font changed from 'sans-serif' to 'monospace'

## v1.0.0 - 2021-08-24
- Initial Release
