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

    seqs1_content_map = {v: k for k, v in seqs1.items()}
    seqs2_content_map = {v: k for k, v in seqs2.items()}

    common_seqs = set(seqs1.values()) & set(seqs2.values())
    unique_seqs1 = set(seqs1.values()) - set(seqs2.values())
    unique_seqs2 = set(seqs2.values()) - set(seqs1.values())

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

    common_links = links1_set & links2_set
    unique_links1 = links1_set - links2_set
    unique_links2 = links2_set - links1_set

    print("=== Sequence Comparison ===")
    print("Common sequences:", len(common_seqs))
    print("Unique to file1:", len(unique_seqs1))
    print("Unique to file2:", len(unique_seqs2))

    print("\n=== Link Comparison ===")
    print("Common links:", len(common_links))
    print("Unique to file1:", len(unique_links1))
    print("Unique to file2:", len(unique_links2))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python compare.py <input_file_1> <input_file_2>")
    else:
        compare(sys.argv[1], sys.argv[2])