import sys
import math

def read_line(file_path):
    try:
        with open(file_path, 'r') as file:
            while True:
                line = file.readline()
                if not line:
                    break
                yield line.strip().split('\t')
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

class stats:
    def __init__(self):
        self.values = []

    def add_value(self, value):
        self.values.append(value)

    def calculate_average(self):
        if not self.values:
            return 0.0
        return sum(self.values) / len(self.values)

    def calculate_std_dev(self):
        if not self.values:
            return 0.0
        avg = self.calculate_average()
        variance = sum((x - avg) ** 2 for x in self.values) / len(self.values)
        return math.sqrt(variance)

def read_until_colon(string):
    result = []
    for char in string:
        if char == ':':
            break
        result.append(char)
    return ''.join(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 analyze.py out.rgfa")
        sys.exit(1)

    file_path = sys.argv[1]
    
    segment_count = 0
    link_count = 0
    segment_length_stats = stats()
    overlapping_count_stats = stats()
    
    with open(file_path, 'r') as file:
        for line in file:
            line_parts = line.strip().split('\t')
            if line_parts[0] == 'S': 
                sequence = line_parts[2]
                segment_length_stats.add_value(len(sequence))
                segment_count += 1
            elif line_parts[0] == 'L':
                overlapping_count = int(line_parts[5][:-1])
                overlapping_count_stats.add_value(overlapping_count)
                link_count += 1

        segment_avg_length = segment_length_stats.calculate_average()
        segment_std_dev = segment_length_stats.calculate_std_dev()

        link_avg_overlap = overlapping_count_stats.calculate_average()
        link_std_dev = overlapping_count_stats.calculate_std_dev()

        print(f"Segment count: {segment_count}")
        print(f"Segment avg length: {segment_avg_length:.5f}")
        print(f"Segment std length: {segment_std_dev:.5f}")
        print(f"Link count: {link_count}")
        print(f"Link overlapping avg length: {link_avg_overlap:.5f}")
        print(f"Link overlapping std length: {link_std_dev:.5f}")
