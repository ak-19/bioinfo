def countPatternOccurence(pattern, text):
    k = len(pattern)
    t = len(text)
    found = []
    for pos in range(0, t - k + 1):
        if pattern == text[pos:pos + k]:
            found.append(pos)
        # print(text[pos:pos + k])
    return [found, len(found)]

[positins, count] = countPatternOccurence('an', 'ante goes an home an an')
print('Found', count, ' matches')
print('At positions', positins)