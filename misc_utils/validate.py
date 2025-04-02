'''
This scripts validates the overlaps' correctness in rGFA/GFA file
'''

import re
import sys

class RGFA:
    def __init__(self):
        self.segments = {}
        self.links = []

    def add_segment(self, line):
        fields = line.split("\t")
        segment_name = fields[1]
        sequence = fields[2].strip()
        if len(sequence) == 0:
            print(f"Error: Sequence {segment_name} is empty")
        if segment_name in self.segments:
            print(f"Error: Duplicate id for sequence {segment_name} detected")
        self.segments[segment_name] = sequence

    def add_link(self, line):
        fields = line.split("\t")
        from_segment = fields[1]
        to_segment = fields[3]
        cigar_string = fields[5].strip()

        overlap_length = int(re.match(r'(\d+)M$', cigar_string).group(1))

        self.links.append({
            'from': from_segment,
            'to': to_segment,
            'overlap_length': overlap_length
        })

    def check_link_overlap(self):
        for link in self.links:
            from_segment = self.segments.get(link['from'])
            to_segment = self.segments.get(link['to'])

            if not from_segment or not to_segment:
                print(f"Warning: Missing segments for link between {link['from']} and {link['to']}")
                continue

            from_len = len(from_segment)
            to_len = len(to_segment)

            overlap_correct = self.check_overlap(from_segment, to_segment, link['overlap_length'])
            if not overlap_correct:
                print(f"Error: Overlap mismatch for link {link['from']} -> {link['to']} with length {link['overlap_length']}")

    def check_overlap(self, from_segment, to_segment, overlap_length):
        if overlap_length == 0 or from_segment[-overlap_length:] == to_segment[:overlap_length]:
            return True
        return False

    def process_rgfa(self, rgfa_file):
        with open(rgfa_file, 'r') as file:
            for line in file:
                if line.startswith("S"):
                    self.add_segment(line)
                elif line.startswith("L"):
                    self.add_link(line)

        self.check_link_overlap()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python validate.py <input_file>")
    else:
        rgfa = RGFA()
        rgfa.process_rgfa(sys.argv[1])