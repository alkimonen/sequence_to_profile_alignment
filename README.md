# Sequence to Profile Alignment
This project implements sequence to profile alignment in **C**. It takes n DNA sequences, 2 <= n <= 10, in a single aln-formatted
(multiple alignment formatted) file and another DNA sequence. It utilizes naive implementation of Needleman-Wunsh alignment algorithm. It takes gap penalty, match score, and mismatch penalty via parameters.

### Parameters:
- **--aln:** Only one aln-formatted file containing all given alignments "aligned sequences.aln",
which contains n DNA sequences line-by-line.
- **--fasta:** Sequence fasta file to be aligned to the given profile.
- **--gap:** gap penalty score.
- **--match:** matching score.
- **--mismatch:** mismatch penalty score.

### Output
- **--out:** sequence.aln

### Compilation
Makefile uses **Makefile** and **gcc** as compiler. These packages should be installed into the system in order to run the project.
After navigating to the projects directory to compile the project execute
```
make
```
For test purposes, there are examples sequences in *test* folder. After compilation, an example to run the project can be
```
alignSeqToProfile \
  --fasta seq. fasta \
  --aln aligned_sequences .aln \
  --out seq. aln \
  --gap 1 \
  --match 2 \
  --mismatch 3
```

