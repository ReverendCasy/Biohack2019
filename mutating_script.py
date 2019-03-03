from Bio import SeqIO
from Bio import SeqUtils
import os
from pprint import pprint
from collections import defaultdict
import random

os.chdir('/home/reverend_casy/Biohack/Mutator')
toxin_file = 'toxin_cites.txt'
rec_file = 'rec_cites.txt'


def converter(file, verse, adjust):
    op = open(file, 'r')
    cites = str()
    for line in op:
        if line.find(verse) == 0:
            cites = line
    op.close()
    cites = cites.replace(verse, '').replace('[', '').replace(']', '').replace('\'', '').replace('(', '').replace(')', '').replace(' ', '').split(',')
    loc_dict = dict()
    for char in range(0, len(cites), 2):
        loc_dict[int(cites[char]) + adjust] = SeqUtils.seq1(cites[char + 1])
    return loc_dict

def grouping_parsing(file, adjust1, adjust2):
    grouping_dict = defaultdict(list)
    op = open(file, 'r')
    for line in op:
        line = line.rstrip()
        ultra_key = int(line[0:3])
        newline = line[3:].replace('[', '').replace(']', '').replace('\'', '').replace('(', '').replace(')', '').replace(' ', '').split(',')
        # print(newline)
        for char in range(0, len(newline), 2):
            if newline[0] != '':
                grouping_dict[ultra_key + adjust1].append(int(newline[char]) + adjust2)
            else:
                grouping_dict[ultra_key + adjust1].append(None)
    op.close()

    return grouping_dict


def please_mutate_me(dictionary, indices, sequence, opposite_dictionary):

    def amino_type(char):
        nonpolar = ['G', 'A', 'V', 'L', 'I', 'M', 'W', 'F', 'P']
        polar = ['S', 'T', 'C', 'Y', 'N', 'Q']
        positive = ['K', 'R', 'H']
        negative = ['D', 'E']
        for i in [nonpolar, polar, positive, negative]:
            if char in i:
                i.remove(char)
        return [nonpolar, polar, positive, negative]

    def random_mutation(char):
        substitutions = list()
        for nested_list in amino_type(char):
            substitutions.append(random.choice(nested_list))
        return substitutions

    ind_list = list(indices.keys())
    # print(ind_list)
    if max(ind_list) + 16 <= len(seq):
        borders = list(range(min(ind_list) - 15, max(ind_list) + 16))
    else:
        borders = list(range(min(ind_list) - 15, len(seq)))
    opp_ind_list = list(opposite_dictionary.keys())
    if max(opp_ind_list) + 16 <= len(seq):
        opp_borders = list(range(min(opp_ind_list) - 15, max(opp_ind_list) + 16))
    else:
        opp_borders = list(range(min(opp_ind_list) - 15, len(seq)))

    dat_fukkin_list = []

    ''' Unary query '''
    for i in ind_list:
        for j in random_mutation(dictionary[i]):
            query_targeted = f"{dictionary[i]}{i}{j}"
            # print(query_targeted)
            dat_fukkin_list.append(query_targeted)
        exclusion_list = borders.copy()
        exclusion_list.remove(i)
        rand_num = random.choice(exclusion_list)
        for k in random_mutation(sequence[rand_num]):
            query_randomized = f"{sequence[rand_num]}{rand_num}{k}"
            # print(query_randomized)
            dat_fukkin_list.append(query_randomized)

    ''' Binary query '''
    for i in ind_list:
        query_stranded_list = []
        query_opposite_list = []
        for j in random_mutation(dictionary[i]):
            query_stranded = f"{dictionary[i]}{i}{j}"
            query_stranded_list.append(query_stranded)
        if indices[i][0]:
            random_opposite = random.choice(indices[i])
            for k in random_mutation(sequence[random_opposite]):
                query_opposite = f"{sequence[random_opposite]}{random_opposite}{k}"
                query_opposite_list.append(query_opposite)
        for x in range(0, len(query_stranded_list)):
            if len(query_opposite_list) > 1:
                result = f"{query_stranded_list[x]},{query_opposite_list[x]}"
            else:
                result = query_stranded_list[x]
            # print(result)
            dat_fukkin_list.append(result)

    ''' Tertiary query '''
    for i in ind_list:
        query_stranded_list = []
        query_opposite_list_1 = []
        query_opposite_list_2 = []
        for j in random_mutation(dictionary[i]):
            query_stranded = f"{dictionary[i]}{i}{j}"
            query_stranded_list.append(query_stranded)
        if indices[i][0]:
            opposite = random.choice(indices[i])
            for l in random_mutation(sequence[opposite]):
                query_opposite = f"{sequence[opposite]}{opposite}{l}"
                query_opposite_list_1.append(query_opposite)
            exclusion_list = opp_borders.copy()
            exclusion_list.remove(opposite)
            print(opposite)
            rand_num = random.choice(exclusion_list)
            for m in random_mutation(sequence[rand_num]):
                query_opposite_random = f"{sequence[rand_num]}{rand_num}{m}"
                query_opposite_list_2.append(query_opposite_random)
        for x in range(0, len(query_stranded_list)):
            if len(query_opposite_list_1) > 1:
                result = f"{query_stranded_list[x]},{query_opposite_list_1[x]},{query_opposite_list_2[x]}"
            else:
                result = query_stranded_list[x]
            # print(result)
            dat_fukkin_list.append(result)

    ''' Quaternary query '''
    for i in ind_list:
        query_stranded_list_1 = []
        query_stranded_list_2 = []
        query_opposite_list_1 = []
        query_opposite_list_2 = []
        for j in random_mutation(dictionary[i]):
            query_stranded = f"{dictionary[i]}{i}{j}"
            query_stranded_list_1.append(query_stranded)
        exclusion_list = borders.copy()
        exclusion_list.remove(i)
        rand_num = random.choice(exclusion_list)
        for k in random_mutation(sequence[rand_num]):
            query_randomized = f"{sequence[rand_num]}{rand_num}{k}"
            query_stranded_list_2.append(query_randomized)
        if indices[i][0]:
            opposite = random.choice(indices[i])
            for l in random_mutation(sequence[opposite]):
                query_opposite = f"{sequence[opposite]}{opposite}{l}"
                query_opposite_list_1.append(query_opposite)
            exclusion_list = opp_borders
            exclusion_list.copy().remove(opposite)
            rand_num = random.choice(exclusion_list)
            for m in random_mutation(sequence[rand_num]):
                query_opposite_random = f"{sequence[rand_num]}{rand_num}{m}"
                query_opposite_list_2.append(query_opposite_random)
        for x in range(0, len(query_stranded_list)):
            if len(query_opposite_list_1) > 1:
                result = f"{query_stranded_list_1[x]},{query_stranded_list_2[x]},{query_opposite_list_1[x]},{query_opposite_list_2[x]}"
            else:
                result = query_stranded_list[x]
            # print(result)
            dat_fukkin_list.append(result)

    return  dat_fukkin_list


tox_dict = converter(toxin_file, '00002', -33)
rec_dict = converter(rec_file, '', 436)

toxin = list(SeqIO.parse('target_sequences.fasta', 'fasta'))[0].seq
rec = list(SeqIO.parse('target_sequences.fasta', 'fasta'))[1].seq

seq = str(toxin + rec)
toxin_grouping = grouping_parsing('toxin_surroundings.txt', -33, 436)
rec_grouping = grouping_parsing('ligand_surroundings.txt', 436, -33)
# print(toxin_grouping)

print(tox_dict)
print(toxin_grouping)
print(len(seq))

tox_fukk = please_mutate_me(tox_dict, toxin_grouping, seq, rec_dict)
rec_fukk = please_mutate_me(rec_dict, rec_grouping, seq, tox_dict)
print(tox_fukk)


for i in range(0, 10):
    with open('delta_gibbs_queries_toxin_target.txt', 'a') as handle:
        handle.write('\n'.join(tox_fukk))

    with open('delta_gibbs_queries_rec_target.txt', 'a') as handle:
        handle.write('\n'.join(rec_fukk))
