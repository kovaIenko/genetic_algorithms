chr1 = [0, 1, 1]
chr2 = [0, 1, 0]
result = []

for i in range(len(chr1)):
    result.append(chr1[i] ^ chr2[i])

print(result)