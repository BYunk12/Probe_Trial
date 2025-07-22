# NCBI-based Probe Designer

This program identifies unique 21-30 nucleotide segments from mRNA sequences using the NCBI database and generates reverse complement probe sequences for molecular biology applications.

## Features

- **NCBI Integration**: Fetches mRNA sequences directly from NCBI database
- **Flexible Input**: Search by gene name or NCBI accession number
- **Uniqueness Checking**: Multiple levels of uniqueness verification
- **Biochemical Filtering**: Filters probes by melting temperature and GC content
- **Reverse Complement Generation**: Automatically generates probe sequences
- **BLAST Integration**: Optional BLAST-based uniqueness verification
- **Customizable Parameters**: Adjustable probe length, temperature, and GC content ranges

## Installation

1. Install Python dependencies:
```bash
pip install -r requirements.txt
```

2. Make the script executable:
```bash
chmod +x probe_designer.py
```

## Usage

### Basic Usage

**Search by gene name:**
```bash
python probe_designer.py --gene "GAPDH" --email "your.email@example.com"
```

**Search by NCBI accession number:**
```bash
python probe_designer.py --accession "NM_002046" --email "your.email@example.com"
```

### Advanced Options

**Custom probe length range:**
```bash
python probe_designer.py --gene "ACTB" --email "your.email@example.com" --min-length 25 --max-length 35
```

**Enable BLAST uniqueness checking (slower but more thorough):**
```bash
python probe_designer.py --gene "TP53" --email "your.email@example.com" --use-blast
```

**Specify organism:**
```bash
python probe_designer.py --gene "GAPDH" --email "your.email@example.com" --organism "Mus musculus"
```

**Save results to file:**
```bash
python probe_designer.py --gene "GAPDH" --email "your.email@example.com" --output "gapdh_probes.txt"
```

**Limit number of probes:**
```bash
python probe_designer.py --gene "GAPDH" --email "your.email@example.com" --max-probes 5
```

### Complete Command Reference

```bash
python probe_designer.py [OPTIONS]

Required Arguments:
  --email EMAIL         Your email address (required by NCBI)

Input (choose one):
  --gene GENE_NAME      Gene name to search for
  --accession ACC_NUM   NCBI accession number

Optional Arguments:
  --organism ORGANISM   Target organism (default: "Homo sapiens")
  --min-length MIN      Minimum probe length (default: 21)
  --max-length MAX      Maximum probe length (default: 30)
  --use-blast          Use BLAST for uniqueness checking (slower)
  --max-probes NUM     Maximum number of probes to return (default: 10)
  --output FILE        Output file for probe sequences
```

## How It Works

### 1. Sequence Retrieval
The program searches NCBI's nucleotide database for mRNA sequences matching your gene name or fetches a specific sequence by accession number.

### 2. Segment Generation
All possible segments of length 21-30 nucleotides (configurable) are generated from the retrieved sequences.

### 3. Uniqueness Filtering
- **Internal uniqueness**: Removes segments that appear multiple times within the same sequence
- **Optional BLAST checking**: Verifies uniqueness against the entire NCBI database

### 4. Biochemical Filtering
Segments are filtered based on:
- **Melting Temperature (Tm)**: 50-65°C (optimal for PCR/hybridization)
- **GC Content**: 40-60% (balanced stability)

### 5. Probe Generation
For each qualifying segment, the program generates:
- **Target sequence**: The original mRNA segment
- **Probe sequence**: Reverse complement for hybridization

## Output Format

The program outputs detailed information for each probe:

```
Probe 1:
  Target Sequence (5' -> 3'): ATGGCACCGTCAAGGCTGAGAAC
  Probe Sequence (5' -> 3'):  GTTCTCAGCCTTGACGGTGCCAT
  Source: NM_002046.7
  Position: 123-146
  Length: 23 bp
  Melting Temperature: 58.0°C
  GC Content: 52.2%
```

## Important Notes

### NCBI Requirements
- **Email address**: Required by NCBI for API access
- **Rate limiting**: The program includes delays between BLAST queries to respect NCBI's servers
- **Internet connection**: Required for NCBI database access

### Performance Considerations
- **BLAST checking**: Significantly slower but more accurate for uniqueness verification
- **Large genes**: May generate many segments; consider using `--max-probes` to limit output
- **Network dependent**: Performance depends on internet connection and NCBI server load

### Limitations
- **Simple Tm calculation**: Uses basic formula; for critical applications, verify with specialized software
- **Uniqueness scope**: Internal uniqueness only checks within retrieved sequences unless BLAST is enabled
- **Species specificity**: Default searches human sequences; specify `--organism` for other species

## Example Workflows

### Quick probe design for common gene:
```bash
python probe_designer.py --gene "GAPDH" --email "researcher@university.edu"
```

### Comprehensive analysis with BLAST verification:
```bash
python probe_designer.py --gene "TP53" --email "researcher@university.edu" --use-blast --output "tp53_probes.txt"
```

### Mouse-specific probes:
```bash
python probe_designer.py --gene "Actb" --email "researcher@university.edu" --organism "Mus musculus"
```

### Custom probe specifications:
```bash
python probe_designer.py --gene "MYC" --email "researcher@university.edu" --min-length 28 --max-length 32 --max-probes 3
```

## Troubleshooting

### Common Issues:
1. **No sequences found**: Check gene name spelling and organism
2. **Network errors**: Verify internet connection and try again
3. **No suitable probes**: Adjust filtering criteria or try a different gene region
4. **BLAST timeouts**: Use `--use-blast` sparingly; it's slower but more thorough

### Error Messages:
- `"No mRNA sequences found"`: Gene name not found in NCBI database
- `"Error fetching sequences"`: Network or NCBI server issue
- `"No suitable probes found"`: All segments filtered out by criteria

## Contributing

This tool can be extended with additional features such as:
- Secondary structure prediction
- Off-target analysis
- Primer3 integration
- Custom filtering criteria
- Batch processing capabilities

## License

This program is provided as-is for research and educational purposes.