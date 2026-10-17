# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**Codon Harmony** is a bioinformatics tool for amino acid reverse translation and DNA sequence optimization based on species-specific codon usage profiles. It optimizes DNA sequences to match host organism codon usage while removing unwanted sequence features.

## Commands

### Installation
```bash
pip install -e .
```

### Running the CLI
```bash
codon_harmony --input input.fasta --output output.fasta --host "Escherichia coli"
```

Key arguments:
- `--input`: Input FASTA file with amino acid sequences
- `--output`: Output FASTA file for DNA sequences
- `--host`: Host organism name or NCBI Taxonomy ID (default: "413997" = E. coli B)
- `--host-threshold`: Lowest codon fraction per AA allowed (default: 0.10)
- `--cycles`: Number of independent codon samples (default: 10)
- `--inner-cycles`: Iterations per sample (default: 10)
- `--max-relax`: Maximum deviation from host profile (default: 0.1)
- `--restriction-enzymes`: List of restriction enzyme sites to remove
- `--remove-splice-sites` / `--no-remove-splice-sites`: Toggle splice site removal
- `--remove-start-sites` / `--no-remove-start-sites`: Toggle alternate start site removal

### Running Tests
```bash
pytest tests/
```

Run a single test file:
```bash
pytest tests/test_codon_harmony.py
```

Run a specific test:
```bash
pytest tests/test_codon_harmony.py::test_function_name
```

### Linting
```bash
flake8
```

## Architecture

The project has a simple, flat structure:

```
codon_harmony/
├── codon_harmony.py    # Main entry point, CLI parser, harmonization workflow
├── util/               # Utility modules
│   ├── codon_use.py    # Codon usage table handling, host profile generation
│   ├── seq_opt.py     # DNA sequence optimization (GC content, hairpins, etc.)
│   ├── seq.py         # Sequence utilities
│   ├── codon_table_downloader.py
│   └── logging.py     # Custom logging with verbosity levels
└── data/
    └── __init__.py    # Data constants (GC thresholds, ribosome binding sites, restriction enzymes)
```

### Core Workflow

1. **Reverse Translation** (`codon_harmony.py:170`): Convert amino acids to DNA using Biopython's back_translate
2. **Codon Resampling** (`seq_opt.py:49`): Replace codons based on host species' codon usage frequencies
3. **Profile Harmonization** (`seq_opt.py:175`): Adjust codon usage to match host profile
4. **Sequence Cleanup** (`codon_harmony.py:196-215`): Iteratively remove:
   - High/low GC content regions
   - Alternate start sites
   - Repeating sequences (3+ consecutive identical codons, 9-mer repeats)
   - Local homopolymers
   - RNA hairpins
   - Restriction enzyme sites
   - RNA splice sites (optional)
5. **Ranking**: Select best sequence by Codon Adaptation Index (CAI)

### Key Dependencies

- `biopython >= 1.73`: Sequence handling, codon tables
- `python_codon_tables >= 0.1.8`: Species-specific codon usage data
