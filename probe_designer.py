#!/usr/bin/env python3
"""
NCBI-based Probe Designer

This program identifies unique 21-30 nucleotide segments from mRNA sequences
using the NCBI database and generates reverse complement probe sequences.
"""

import os
import sys
from typing import List, Tuple, Set, Dict
import argparse
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import time
import random
from collections import defaultdict

class ProbeDesigner:
    def __init__(self, email: str, min_length: int = 21, max_length: int = 30):
        """
        Initialize the probe designer.
        
        Args:
            email: Your email address (required by NCBI)
            min_length: Minimum probe length (default: 21)
            max_length: Maximum probe length (default: 30)
        """
        self.email = email
        self.min_length = min_length
        self.max_length = max_length
        Entrez.email = email
        
    def fetch_sequence_by_gene_name(self, gene_name: str, organism: str = "Homo sapiens") -> List[SeqRecord]:
        """
        Fetch mRNA sequences for a given gene name from NCBI.
        
        Args:
            gene_name: Name of the gene
            organism: Target organism (default: Homo sapiens)
            
        Returns:
            List of SeqRecord objects containing mRNA sequences
        """
        print(f"Searching for {gene_name} mRNA sequences in {organism}...")
        
        # Search for mRNA sequences
        search_term = f"{gene_name}[Gene Name] AND {organism}[Organism] AND mRNA[Filter]"
        
        try:
            handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=10)
            search_results = Entrez.read(handle)
            handle.close()
            
            if not search_results["IdList"]:
                print(f"No mRNA sequences found for {gene_name}")
                return []
            
            print(f"Found {len(search_results['IdList'])} sequences")
            
            # Fetch the actual sequences
            handle = Entrez.efetch(db="nucleotide", id=search_results["IdList"], rettype="fasta", retmode="text")
            sequences = list(SeqIO.parse(handle, "fasta"))
            handle.close()
            
            return sequences
            
        except Exception as e:
            print(f"Error fetching sequences: {e}")
            return []
    
    def fetch_sequence_by_accession(self, accession: str) -> SeqRecord:
        """
        Fetch a sequence by its accession number.
        
        Args:
            accession: NCBI accession number
            
        Returns:
            SeqRecord object
        """
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            sequence = SeqIO.read(handle, "fasta")
            handle.close()
            return sequence
        except Exception as e:
            print(f"Error fetching sequence {accession}: {e}")
            return None
    
    def generate_segments(self, sequence: str) -> List[Tuple[str, int, int]]:
        """
        Generate all possible segments of specified length from the sequence.
        
        Args:
            sequence: Input DNA/RNA sequence
            
        Returns:
            List of tuples (segment, start_pos, end_pos)
        """
        segments = []
        
        for length in range(self.min_length, self.max_length + 1):
            for i in range(len(sequence) - length + 1):
                segment = sequence[i:i + length]
                segments.append((segment, i, i + length))
        
        return segments
    
    def check_uniqueness_blast(self, segment: str, max_hits: int = 5) -> bool:
        """
        Check if a segment is unique using NCBI BLAST.
        
        Args:
            segment: DNA/RNA segment to check
            max_hits: Maximum number of BLAST hits to consider unique
            
        Returns:
            True if segment appears to be unique (few BLAST hits)
        """
        try:
            print(f"  Checking uniqueness of segment: {segment[:20]}...")
            
            # Perform BLAST search
            result_handle = NCBIWWW.qblast("blastn", "nt", segment, expect=1000, hitlist_size=max_hits + 5)
            blast_records = NCBIXML.parse(result_handle)
            
            hit_count = 0
            for blast_record in blast_records:
                hit_count += len(blast_record.alignments)
                
                # Check if we have too many hits
                if hit_count > max_hits:
                    print(f"    Too many hits ({hit_count}), not unique")
                    return False
            
            print(f"    Found {hit_count} hits, appears unique")
            return True
            
        except Exception as e:
            print(f"    Error in BLAST search: {e}")
            return False
    
    def simple_uniqueness_check(self, segments: List[Tuple[str, int, int]]) -> List[Tuple[str, int, int]]:
        """
        Perform a simple uniqueness check by looking for duplicates within the input.
        
        Args:
            segments: List of segments to check
            
        Returns:
            List of segments that appear only once in the input
        """
        segment_counts = defaultdict(int)
        segment_info = {}
        
        # Count occurrences
        for segment, start, end in segments:
            segment_counts[segment] += 1
            if segment not in segment_info:
                segment_info[segment] = (start, end)
        
        # Return only unique segments
        unique_segments = []
        for segment, count in segment_counts.items():
            if count == 1:
                start, end = segment_info[segment]
                unique_segments.append((segment, start, end))
        
        return unique_segments
    
    def calculate_tm(self, sequence: str) -> float:
        """
        Calculate melting temperature using a simple approximation.
        
        Args:
            sequence: DNA/RNA sequence
            
        Returns:
            Estimated melting temperature in Celsius
        """
        # Simple Tm calculation: Tm = 2*(A+T) + 4*(G+C)
        sequence = sequence.upper()
        at_count = sequence.count('A') + sequence.count('T') + sequence.count('U')
        gc_count = sequence.count('G') + sequence.count('C')
        
        return 2 * at_count + 4 * gc_count
    
    def calculate_gc_content(self, sequence: str) -> float:
        """
        Calculate GC content as a percentage.
        
        Args:
            sequence: DNA/RNA sequence
            
        Returns:
            GC content as percentage
        """
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        return (gc_count / len(sequence)) * 100 if sequence else 0
    
    def get_reverse_complement(self, sequence: str) -> str:
        """
        Generate reverse complement of a DNA/RNA sequence.
        
        Args:
            sequence: Input DNA/RNA sequence
            
        Returns:
            Reverse complement sequence
        """
        # Convert RNA to DNA for reverse complement calculation
        dna_sequence = sequence.replace('U', 'T')
        seq_obj = Seq(dna_sequence)
        return str(seq_obj.reverse_complement())
    
    def filter_segments_by_criteria(self, segments: List[Tuple[str, int, int]], 
                                  min_tm: float = 50.0, max_tm: float = 65.0,
                                  min_gc: float = 40.0, max_gc: float = 60.0) -> List[Tuple[str, int, int, Dict]]:
        """
        Filter segments based on melting temperature and GC content.
        
        Args:
            segments: List of segments to filter
            min_tm: Minimum melting temperature
            max_tm: Maximum melting temperature
            min_gc: Minimum GC content percentage
            max_gc: Maximum GC content percentage
            
        Returns:
            List of filtered segments with their properties
        """
        filtered_segments = []
        
        for segment, start, end in segments:
            tm = self.calculate_tm(segment)
            gc_content = self.calculate_gc_content(segment)
            
            if min_tm <= tm <= max_tm and min_gc <= gc_content <= max_gc:
                properties = {
                    'tm': tm,
                    'gc_content': gc_content,
                    'length': len(segment)
                }
                filtered_segments.append((segment, start, end, properties))
        
        return filtered_segments
    
    def design_probes(self, gene_name: str = None, accession: str = None, 
                     organism: str = "Homo sapiens", use_blast: bool = False,
                     max_probes: int = 10) -> List[Dict]:
        """
        Main function to design probes for a given gene or sequence.
        
        Args:
            gene_name: Name of the gene to design probes for
            accession: NCBI accession number (alternative to gene_name)
            organism: Target organism
            use_blast: Whether to use BLAST for uniqueness checking (slower but more accurate)
            max_probes: Maximum number of probes to return
            
        Returns:
            List of probe dictionaries with sequences and properties
        """
        # Fetch sequences
        if gene_name:
            sequences = self.fetch_sequence_by_gene_name(gene_name, organism)
        elif accession:
            seq = self.fetch_sequence_by_accession(accession)
            sequences = [seq] if seq else []
        else:
            raise ValueError("Either gene_name or accession must be provided")
        
        if not sequences:
            return []
        
        all_probes = []
        
        for seq_record in sequences:
            print(f"\nProcessing sequence: {seq_record.id}")
            print(f"Length: {len(seq_record.seq)} bp")
            
            # Generate all possible segments
            segments = self.generate_segments(str(seq_record.seq))
            print(f"Generated {len(segments)} segments")
            
            # Check for uniqueness within the sequence
            unique_segments = self.simple_uniqueness_check(segments)
            print(f"Found {len(unique_segments)} internally unique segments")
            
            # Filter by biochemical criteria
            filtered_segments = self.filter_segments_by_criteria(unique_segments)
            print(f"Found {len(filtered_segments)} segments passing biochemical criteria")
            
            # Optionally check uniqueness with BLAST
            if use_blast and filtered_segments:
                print("Performing BLAST uniqueness checks...")
                blast_unique = []
                for segment, start, end, props in filtered_segments[:min(5, len(filtered_segments))]:  # Limit BLAST queries
                    if self.check_uniqueness_blast(segment):
                        blast_unique.append((segment, start, end, props))
                    time.sleep(2)  # Be nice to NCBI servers
                
                filtered_segments = blast_unique
                print(f"Found {len(filtered_segments)} BLAST-verified unique segments")
            
            # Convert to probe format
            for segment, start, end, props in filtered_segments[:max_probes]:
                probe_sequence = self.get_reverse_complement(segment)
                
                probe = {
                    'target_sequence': segment,
                    'probe_sequence': probe_sequence,
                    'source_sequence': seq_record.id,
                    'start_position': start,
                    'end_position': end,
                    'length': props['length'],
                    'melting_temperature': props['tm'],
                    'gc_content': props['gc_content']
                }
                all_probes.append(probe)
        
        # Sort by melting temperature for consistency
        all_probes.sort(key=lambda x: x['melting_temperature'])
        
        return all_probes[:max_probes]

def main():
    parser = argparse.ArgumentParser(description="Design specific probes for mRNA sequences")
    parser.add_argument("--gene", type=str, help="Gene name to search for")
    parser.add_argument("--accession", type=str, help="NCBI accession number")
    parser.add_argument("--email", type=str, required=True, help="Your email address (required by NCBI)")
    parser.add_argument("--organism", type=str, default="Homo sapiens", help="Target organism")
    parser.add_argument("--min-length", type=int, default=21, help="Minimum probe length")
    parser.add_argument("--max-length", type=int, default=30, help="Maximum probe length")
    parser.add_argument("--use-blast", action="store_true", help="Use BLAST for uniqueness checking (slower)")
    parser.add_argument("--max-probes", type=int, default=10, help="Maximum number of probes to return")
    parser.add_argument("--output", type=str, help="Output file for probe sequences")
    
    args = parser.parse_args()
    
    if not args.gene and not args.accession:
        parser.error("Either --gene or --accession must be provided")
    
    # Create probe designer
    designer = ProbeDesigner(args.email, args.min_length, args.max_length)
    
    # Design probes
    print("Starting probe design...")
    probes = designer.design_probes(
        gene_name=args.gene,
        accession=args.accession,
        organism=args.organism,
        use_blast=args.use_blast,
        max_probes=args.max_probes
    )
    
    if not probes:
        print("No suitable probes found.")
        return
    
    # Display results
    print(f"\n{'='*80}")
    print(f"PROBE DESIGN RESULTS")
    print(f"{'='*80}")
    print(f"Found {len(probes)} suitable probes:\n")
    
    output_lines = []
    for i, probe in enumerate(probes, 1):
        probe_info = f"""
Probe {i}:
  Target Sequence (5' -> 3'): {probe['target_sequence']}
  Probe Sequence (5' -> 3'):  {probe['probe_sequence']}
  Source: {probe['source_sequence']}
  Position: {probe['start_position']}-{probe['end_position']}
  Length: {probe['length']} bp
  Melting Temperature: {probe['melting_temperature']:.1f}°C
  GC Content: {probe['gc_content']:.1f}%
"""
        print(probe_info)
        output_lines.append(probe_info)
    
    # Save to file if requested
    if args.output:
        with open(args.output, 'w') as f:
            f.write(f"Probe Design Results for {args.gene or args.accession}\n")
            f.write("="*80 + "\n")
            for line in output_lines:
                f.write(line + "\n")
        print(f"\nResults saved to {args.output}")

if __name__ == "__main__":
    main()