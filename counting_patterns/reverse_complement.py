def reverse_complement(text, pattern):
    pattern_len = len(pattern)
    text_len = len(text)
    found_at_positions = []
    replacement_map = {'A':'T', 'T':'A', 'C':'G', 'G': 'C'}
    complement_content = pattern[::-1]
    complement_content = ''.join(replacement_map[c] for c in complement_content)

    for pos in range(0, text_len - pattern_len + 1):
        current = text[pos:pos + pattern_len]
        if complement_content == current:
            found_at_positions.append(pos)

    return found_at_positions


print(reverse_complement('ACGAACGTGCTGA', 'CGT'))