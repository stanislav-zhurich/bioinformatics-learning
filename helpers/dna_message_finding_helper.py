import numpy as np
import random

def find_hamming_distance(string1, string2):
    distance = 0
    for i in range(0, len(string1)):
        if string1[i] != string2[i]:
            distance += 1         
    return distance

def find_approximate_pattern_matching(genome, pattern, d):
    pattern_len = len(pattern)
    positions = []
    for i in range(0, len(genome) - (pattern_len - 1)):
        window = genome[i: i + pattern_len]
        if find_hamming_distance(pattern, window) <= d:
            positions.append(i)
    return positions

def find_approximate_pattern_count(genome, pattern, d):
    positions = find_approximate_pattern_matching(genome, pattern, d)
    return len(positions)

def find_immediate_neighbors(pattern, d):
    neighbors = set([pattern])
    nucleotides = ['A', "T", "C", "G"]
    if d == 0:
        return neighbors
    
    for j in range(0, d):
        for i in range(0, len(pattern)):
            for nucleotide in nucleotides:
                if nucleotide != pattern[i]:
                    neighbor = pattern[:i] + nucleotide + pattern[i+1:]
                    result = find_immediate_neighbors(neighbor, d - 1)
                    neighbors.update(result)
    return neighbors

def get_motiffs(dna_list, k, d):
    candidates = set()
    result = set()
    for dna in dna_list:
        for i in range(0, len(dna) - (k-1)):
            pattern = dna[i: i + k]
            neighbors = find_immediate_neighbors(pattern, d)
            candidates.update(neighbors)
    for candidate in candidates:
        occurancy = True
        for dna in dna_list:
            count = find_approximate_pattern_count(dna, candidate, d)
            if count == 0:
                occurancy = False
                break
        if occurancy == True:
            result.add(candidate)
    return result

def find_patterns(dna_list, k):
    joined_dna = "".join(dna_list)
    patterns = set()
    for i in range(0, len(joined_dna) - (k - 1)):
        pattern = joined_dna[i: i + k]
        patterns.add(pattern)
    return patterns

def find_min_distance(dna, pattern):
    minimum = 1000000
    for i in range(0, len(dna) - (len(pattern) - 1)):
        occurance = dna[i: i + len(pattern)]
        distance = find_hamming_distance(pattern, occurance)
        if distance < minimum:
            minimum = distance
    return minimum

def find_distance_between_pattern_and_strings(pattern, dna_list):
    distance = 0
    for dna in dna_list:
        distance = distance + find_min_distance(dna, pattern)
    return distance

def find_median_string(dna_list, k):
    patterns = find_patterns(dna_list, k)
    pattern_d_map = {}
    for pattern in patterns:
        d_pattern = 0
        for dna in dna_list:
            current_min = find_min_distance(dna, pattern)
            d_pattern = d_pattern + current_min
        pattern_d_map[pattern] = d_pattern
    return dict(sorted(pattern_d_map.items(), key=lambda item: item[1]))

def find_pattern_probability(array, pattern):
    index_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    prob = 1
    for i in range(0, len(pattern)):
        prob = prob * array[index_map[pattern[i]]][i]
    return prob

def find_most_probable_k_mer(dna, matrix, k):
    array = np.asarray(matrix)
    best_prob = 0
    most_prob_pattern = None
    for i in range(0, len(dna) - (k - 1)):
        pattern = dna[i: i + k]
        prob = find_pattern_probability(array, pattern)
       
        if prob > best_prob:
            most_prob_pattern = pattern
            best_prob= prob
    return most_prob_pattern

def get_profile_of_pattern(k_mer_list):
    
    k_mer_len = len(k_mer_list[0])
    count = np.zeros((4, k_mer_len))
    for i in range(0, k_mer_len):
        column = [s[i] for s in k_mer_list]
        count[0, i] = column.count("A") + 1
        count[1, i] = column.count("C") + 1
        count[2, i] = column.count("G") + 1
        count[3, i] = column.count("T") + 1
    profile = count/np.sum(count, axis = 0)
    return profile

def find_score(k_mer_list, profile):
    index_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    score = 0
    for k_mer in k_mer_list:
        for i in range(0, len(k_mer)):
            nucleotide_index = index_map[k_mer[i]]
            cell_score = profile[nucleotide_index][i]
            score = score + cell_score
    return score

def greedy_motifs_search(dna_list, k, t):
    first_motif = dna_list[0]
    best_score = 0
    best_motifes = None
    for i in range(0, len(first_motif) - (k - 1)):
        k_mer = first_motif[i: i + k]
        motifes = [k_mer]
        profile = get_profile_of_pattern([k_mer])
        
        for j in range(1, len(dna_list)):
            best_k_mer = find_most_probable_k_mer(dna_list[j], profile, k)
            if best_k_mer is None:
                best_k_mer = dna_list[j][0: k]
            motifes.append(best_k_mer)   
            profile = get_profile_of_pattern(motifes)
        
        
        score = find_score(motifes, profile)
        print(f"Motifs: {motifes} with score: {score} ")
        if score > best_score:
            best_motifes = motifes
            best_score = score
    print(best_score)
    return best_motifes

def get_random_k_mers(dna_list, k):
    #random.seed(42)
    random_k_mers = []
    for i in range(0, len(dna_list)):   
        string_len = len(dna_list[i])
        random_number = random.randint(0, string_len - k)
        k_mer = dna_list[i][random_number: random_number + k]
        random_k_mers.append(k_mer)
    return random_k_mers

def score_motifs(motifs):
    score = 0
    for i in range(len(motifs[0])):
        nucleotide_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            nucleotide_counts[motif[i]] += 1
        score += len(motifs) - max(nucleotide_counts.values())
    return score