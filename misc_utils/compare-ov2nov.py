import re
import sys

def compare(filename1, filename2):
    def parse_rgfa_nov(file):
        sequences = {}
        links = []
        overlaps = {}

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
        return sequences, links

    def parse_rgfa_ov(file):
        sequences = {}
        links = []
        overlaps = {}

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
                        overlaps[tgt] = overlap # assign overlap
                    elif overlaps[tgt] != overlap:
                        print(f"Warning: Incoming overlap mismatch. {overlaps[tgt]} {overlap}")
        return sequences, links, overlaps
    
    seqs1, links1, overlaps1 = parse_rgfa_ov(filename1)
    seqs2, links2 = parse_rgfa_nov(filename2)

    links1_content = [
        (seqs1[src][overlaps1.get(src, 0):], seqs1[tgt][overlaps1.get(tgt, 0):])
        for src, tgt in links1
    ]
    
    links2_content = [
        (seqs2.get(src, src), seqs2.get(tgt, tgt))
        for src, tgt in links2
    ]

    links1_set = set(links1_content)
    links2_set = set(links2_content)

    common_links = links1_set & links2_set
    unique_links1 = links1_set - links2_set
    unique_links2 = links2_set - links1_set

    print("\n=== Link Comparison ===")
    print("Links in file1:", len(links1_set))
    print("Links to file2:", len(links2_set))
    print("Common links:", len(common_links))
    print("Unique to file1:", len(unique_links1))
    print("Unique to file2:", len(unique_links2))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare-ov2nov.py <overlaping> <non-overlapping>")
    else:
        compare(sys.argv[1], sys.argv[2])