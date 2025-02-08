import re
import sys

def compare(filename1, filename2):
    def parse_rgfa(file):
        sequences = {}
        links = []
        with open(file) as f:
            for line in f:
                if line.startswith("S"):
                    parts = line.strip().split("\t")
                    sequences[parts[1]] = parts[2]
                elif line[0] == 'L':
                    parts = line.strip().split("\t")
                    links.append((parts[1], parts[3], parts[5].strip()))
        return sequences, links
    
    seqs1, links1 = parse_rgfa(filename1)
    seqs2, links2 = parse_rgfa(filename2)

    common_seqs = len(set(seqs1.values()) & set(seqs2.values()))
    unique_seqs1 = len(set(seqs1.values()) - set(seqs2.values()))
    unique_seqs2 = len(set(seqs2.values()) - set(seqs1.values()))

    links1_content = [
        (seqs1.get(src, src), seqs1.get(tgt, tgt), overlap)
        for src, tgt, overlap in links1
    ]
    links2_content = [
        (seqs2.get(src, src), seqs2.get(tgt, tgt), overlap)
        for src, tgt, overlap in links2
    ]

    links1_set = set(links1_content)
    links2_set = set(links2_content)

    common_links = len(links1_set & links2_set)
    unique_links1 = len(links1_set - links2_set)
    unique_links2 = len(links2_set - links1_set)

    print("=== Sequence Comparison ===")
    print("Links in file1:", len(seqs1))
    print("Links to file2:", len(seqs2))
    print("Common sequences:", common_seqs)
    print("Unique to file1:", unique_seqs1)
    print("Unique to file2:", unique_seqs2)

    print("\n=== Link Comparison ===")
    print("Links in file1:", len(links1_set))
    print("Links to file2:", len(links2_set))
    print("Common links:", common_links)
    print("Unique to file1:", unique_links1)
    print("Unique to file2:", unique_links2)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare-ov2ov.py <input_file_1> <input_file_2>")
    else:
        compare(sys.argv[1], sys.argv[2])