import re
import sys

def compare(filename1, filename2):
    def parse_gfa(file):
        sequences = {}
        links = []
        overlaps = {}
        is_overlap = False

        with open(file) as f:
            for line in f:
                if line.startswith("S"):
                    parts = line.strip().split("\t")
                    sequences[parts[1]] = parts[2]
                elif line[0] == 'L':
                    parts = line.strip().split("\t")
                    src = parts[1]
                    tgt = parts[3]
                    links.append((src, tgt))
                    overlap = int(parts[5].strip()[:-1])
                    if tgt not in overlaps:
                        if 0 < overlap:
                            is_overlap = True 
                        overlaps[tgt] = overlap # assign overlap
                    elif overlaps[tgt] != overlap:
                        print(f"Warning: Incoming overlap mismatch. {overlaps[tgt]} {overlap}")
        if is_overlap:
            for tgt, overlap in overlaps.items():
                sequences[tgt] = sequences[tgt][overlap:]
        return sequences, links
    
    def process_links(seqs, links, file):
        links_content = set()
        for src, tgt in links:
            seq1 = seqs.get(src, "")
            seq2 = seqs.get(tgt, "")
            if seq1 == "" or seq2 == "":
                print(f"Sequences {src}/{tgt} in links not defined in {file}.")
            else:
                links_content.add((seq1, seq2))  
        return links_content      
    
    seqs1, links1 = parse_gfa(filename1)
    seqs2, links2 = parse_gfa(filename2)

    common_seqs = len(set(seqs1.values()) & set(seqs2.values()))
    unique_seqs1 = len(set(seqs1.values()) - set(seqs2.values()))
    unique_seqs2 = len(set(seqs2.values()) - set(seqs1.values()))

    print("=== Sequence Comparison ===")
    print("Sequences in file1:", len(seqs1))
    print("Sequences in file2:", len(seqs2))
    print("Common sequences:", common_seqs)
    print("Unique to file1:", unique_seqs1)
    print("Unique to file2:", unique_seqs2)

    links1_content = process_links(seqs1, links1, filename1)
    links2_content = process_links(seqs2, links2, filename2)

    common_links = len(links1_content & links2_content)
    unique_links1 = len(links1_content - links2_content)
    unique_links2 = len(links2_content - links1_content)
    
    print("\n=== Link Comparison ===")
    print("Links in file1:", len(links1_content))
    print("Links in file2:", len(links2_content))
    print("Common links:", common_links)
    print("Unique to file1:", unique_links1)
    print("Unique to file2:", unique_links2)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare.py <gfa> <gfa>")
    else:
        compare(sys.argv[1], sys.argv[2])