import difflib
filename1="C:/Users/Artlab2/Documents/amir/phd_code/data/yasuda-data/qiime/yasuda-sample.tsv"
filename2= "C:/Users/Artlab2/Documents/amir/phd_code/data/biagi-data/qiime/biagi-sample.tsv"

with open(filename1) as file1, open(filename2) as file2:
    file1_lines = file1.readlines()
    file2_lines = file2.readlines()

    # Find unique lines
    unique_lines = set([line.split()[1] for line in file1_lines]).symmetric_difference([line.split()[1] for line in file2_lines])
    with open('C:/Users/Artlab2/Documents/amir/phd_code/data/unique_lines.txt', 'w') as f:
        f.write('\n'.join(unique_lines))

    # Find common lines
    common_lines = set([line.split()[1] for line in file1_lines]).intersection([line.split()[1] for line in file2_lines])
    with open('C:/Users/Artlab2/Documents/amir/phd_code/data/common_lines.txt', 'w') as f:
        f.write('\n'.join(common_lines))

    # Find longest contiguous matching subsequence
    s = difflib.SequenceMatcher(None, ''.join([line.split()[1] for line in file1_lines]), ''.join([line.split()[1] for line in file2_lines]))
    match = s.find_longest_match(0, len([line.split()[1] for line in file1_lines]), 0, len([line.split()[1] for line in file2_lines]))
    with open('C:/Users/Artlab2/Documents/amir/phd_code/data/longest_contiguous_matching_subsequence.txt', 'w') as f:
        f.write(''.join([line.split()[1] for line in file1_lines])[match.a: match.a + match.size])
