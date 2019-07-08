def most_frequent_pattern_of_length_in(k, text):
    occurances = {}
    max_frequency = 0
    most_frequent_patterns = None 
    for pos in range(0, len(text) + 1 - k):
        pattern = text[pos:pos + k]
        pattern_occured = occurances.get(pattern, 0) + 1
        occurances[pattern] = pattern_occured
        
        if pattern_occured > max_frequency:
            most_frequent_patterns = [pattern]
            max_frequency = pattern_occured
        elif pattern_occured == max_frequency:
            most_frequent_patterns = [*most_frequent_patterns, pattern]
    return [most_frequent_patterns, max_frequency]

print(most_frequent_pattern_of_length_in(4, '''GATCAGCATAAGGGTCCCTGCAATGCATGACAAGCCTGCAGTtGTTTTAC''' ))
