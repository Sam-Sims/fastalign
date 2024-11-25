# fastalign

Super fast multiple sequence alignments against a reference sequence. 

Uses minimap2 (specifically rust bindings to minimap2 via the [minimap2-rs](https://github.com/jguhlin/minimap2-rs) crate) to get alignments between a reference sequence and a set of query sequences, and outputs this as a single aligned fasta file.

You might use this when you have several thousand seqeunces to align, but expect low diverstiy between sequences (i.e an outbreak). In these situations minimap2 does a good enough job, and is much quicker than traditional MSA methods.

## Usage
```
Quick multiple sequnce alignment using minimap2

Usage: fastalign [OPTIONS] --input <Unaligned FASTA> --reference <Reference FASTA> --output <Output FASTA>

Options:
  -i, --input <Unaligned FASTA>      Input (unaligned) FASTA file
  -r, --reference <Reference FASTA>  Input reference FASTA file
  -o, --output <Output FASTA>        Output alignment file
  -t, --threads <Threads>            Number of threads to use. Default: 1 [default: 1]
  -h, --help                         Print help
  -V, --version                      Print version
```

## How it works

1. Reads a reference sequence from the specified FASTA file.
2. Processes the input FASTA file containing the query sequences.
3. Using minimap2, each query sequence is aligned to the reference sequence.
4. The alignment process is parallelised across as many threads as you can give it.
5. Insertions are omitted from the output alignment in order to preserve the reference sequence length, deletions are kept as `-`
6. The aligned sequences are written to the output FASTA file.
