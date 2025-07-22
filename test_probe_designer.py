#!/usr/bin/env python3
"""
Test script for the probe designer

This demonstrates the key functionality without requiring NCBI access.
"""

from probe_designer import ProbeDesigner

def test_basic_functions():
    """Test the basic functions of the probe designer."""
    
    # Initialize probe designer
    designer = ProbeDesigner("test@example.com", min_length=21, max_length=25)
    
    print("Testing ProbeDesigner basic functions...")
    print("=" * 50)
    
    # Test sequence processing functions
    test_sequence = "ATCGATCGATCGATCGAAAAAATTTTTCGCGCGATATCGATCGATCGAAATTTCCCGGGAAATTTCGCGCG"
    
    print(f"Test sequence: {test_sequence}")
    print(f"Length: {len(test_sequence)} bp\n")
    
    # Test segment generation
    segments = designer.generate_segments(test_sequence)
    print(f"Generated {len(segments)} segments (lengths 21-25)")
    print("First 5 segments:")
    for i, (segment, start, end) in enumerate(segments[:5]):
        print(f"  {i+1}. {segment} (pos {start}-{end})")
    print()
    
    # Test uniqueness checking
    unique_segments = designer.simple_uniqueness_check(segments)
    print(f"Found {len(unique_segments)} unique segments")
    print("First 3 unique segments:")
    for i, (segment, start, end) in enumerate(unique_segments[:3]):
        print(f"  {i+1}. {segment} (pos {start}-{end})")
    print()
    
    # Test biochemical filtering
    filtered_segments = designer.filter_segments_by_criteria(unique_segments)
    print(f"Found {len(filtered_segments)} segments passing biochemical criteria")
    
    if filtered_segments:
        print("Best candidate segments:")
        for i, (segment, start, end, props) in enumerate(filtered_segments[:3]):
            tm = props['tm']
            gc = props['gc_content']
            length = props['length']
            probe = designer.get_reverse_complement(segment)
            
            print(f"\nProbe {i+1}:")
            print(f"  Target:  {segment}")
            print(f"  Probe:   {probe}")
            print(f"  Length:  {length} bp")
            print(f"  Tm:      {tm:.1f}°C")
            print(f"  GC:      {gc:.1f}%")
            print(f"  Position: {start}-{end}")
    else:
        print("No segments passed the biochemical criteria.")
    
    print("\n" + "=" * 50)
    print("Basic functionality test completed successfully!")

if __name__ == "__main__":
    test_basic_functions()