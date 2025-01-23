import sys

def read_line(file_path):
    """
    Reads lines from the specified file one by one and splits each by tab characters.
    
    Parameters:
        file_path (str): The path to the file to be read.
    
    Yields:
        list: A list of strings obtained by splitting the current line by tabs.
    """
    try:
        with open(file_path, 'r') as file:
            while True:
                line = file.readline()  # Read one line at a time
                if not line:  # End of file
                    break
                yield line.strip().split('\t')  # Split by tabs and remove whitespace
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


average = 0.0
count = 0

def calculate_avg(value):
    """
    Updates the running average with a new value.

    Parameters:
        value (float): The new value to include in the average.
    """
    global average, count
    count += 1
    average += (value - average) / count  # Update the running average incrementally


def read_until_colon(string):
    """
    Reads a string until the ':' character is reached.
    
    Parameters:
        string (str): The input string to process.
    
    Returns:
        str: The substring up to but not including the ':' character.
    """
    result = []
    for char in string:
        if char == ':':
            break
        result.append(char)
    return ''.join(result)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("No file path provided.")
        sys.exit(1)

    file_path = sys.argv[1]

    segment_count = 0
    link_count = 0

    try:
        for idx, line_parts in enumerate(read_line(file_path)):
            if line_parts[0] == 'S': 
                sequence = read_until_colon(line_parts[2])
                calculate_avg(len(sequence))  # Use len() to get the size of the sequence
                segment_count += 1
            elif line_parts[0] == 'L':
                link_count += 1

        print(f"Total segment count: {segment_count}, AVG segment size: {average:.2f}, Total link count: {link_count}")

    except Exception as e:
        print(f"Error during processing: {e}")
