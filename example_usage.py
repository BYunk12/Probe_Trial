#!/usr/bin/env python3
"""
Example usage of the ProbeDesigner class

This script demonstrates how to use the probe designer programmatically
without command-line arguments.
"""

from probe_designer import ProbeDesigner

def main():
    # Initialize the probe designer
    # Replace with your actual email address
    email = "your.email@example.com"
    designer = ProbeDesigner(email, min_length=21, max_length=30)
    
    # Example 1: Design probes for GAPDH gene
    print("Example 1: Designing probes for GAPDH gene")
    print("=" * 50)
    
    probes = designer.design_probes(
        gene_name="GAPDH",
        organism="Homo sapiens",
        use_blast=False,  # Set to True for more thorough checking (slower)
        max_probes=5
    )
    
    if probes:
        for i, probe in enumerate(probes, 1):
            print(f"Probe {i}:")
            print(f"  Target: {probe['target_sequence']}")
            print(f"  Probe:  {probe['probe_sequence']}")
            print(f"  Tm: {probe['melting_temperature']:.1f}°C")
            print(f"  GC: {probe['gc_content']:.1f}%")
            print()
    else:
        print("No suitable probes found for GAPDH")
    
    # Example 2: Design probes using accession number
    print("\nExample 2: Designing probes using accession number")
    print("=" * 50)
    
    probes = designer.design_probes(
        accession="NM_002046",  # GAPDH mRNA
        use_blast=False,
        max_probes=3
    )
    
    if probes:
        for i, probe in enumerate(probes, 1):
            print(f"Probe {i}:")
            print(f"  Target: {probe['target_sequence']}")
            print(f"  Probe:  {probe['probe_sequence']}")
            print(f"  Position: {probe['start_position']}-{probe['end_position']}")
            print()
    else:
        print("No suitable probes found for NM_002046")

if __name__ == "__main__":
    main()