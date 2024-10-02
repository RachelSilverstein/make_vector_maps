Creating Vector Maps 

- make amino a acid substitutions at certain positions in an annotation of a genbank file
- also creates annotations showing the substituted AAs
- make an inputs file following the format of the provided example

DEPENDENCIES:
- pandas
- biopython

PARAMS IN INPUT FILE:

name
- name of the new vector map you want to make

AAs
- amino acid 1 letter codes of the AAs you are changing to, comma separated

positions
- positions of the AAs you are changing, comma separated (numbering starts from the beginning of the reference annotation ex. "SpCas9"). You need to have an annotation in order to be able to make the substitutions.

reference_annotation
- you need to specify the name of the annotation in the reference map where you are trying to make mutations... make sure this name matches exactly the annotation in the gb file

new_annotation
- what you want to change the name of the annotation to in your new vector map

reference_map
- filename of the reference vector map (needs to be a gb file)
- also provide the path to the file here if it is not in the directory where you're running the script
- Choosing a clean gb file to start with can be a little finicky because the script can't recognize annotations if they span the origin. If your vector map has this, import it to Geneious, change the origin to part of the plasmid that does not have an annotation, and export as gb file again.
- the script also cannot change amino acid positions in a annotation that is on the reverse strand. If this is the case, reverse complement the whole map in Geneious and then re-export it as gb file.

codon_table
- txt file with the mapping of the AAs to the codons used to encode them. The current table provided uses the "Rustymax" codon usage which is the most used codon for each AA in humans. If you use a different codon usage for your AA subs, then you need to make a different table.

RUNNING THE SCRIPT:
python3 make_variant_maps.py