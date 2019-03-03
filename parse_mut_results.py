machine_result = open('/home/lavrentydanilov/Documents/Documents/Scientific_work/Additional_projects/BioHack_2019/data/parse_mutantion/real_data.txt', 'w')
machine_result.write('toxin\treceptor\tG_en\n')

tox = str('YTPIDISLSLTQFLLSEFVPGAGFVLGLVDIIWGIFGPSQWDAFLVQIEQLINQRIEEFA'
           'RNQAISRLEGLSNLYQIYAESFREWEADPTNPALREEMRIQFNDMNSALTTAIPLLAVQN'
            'YQVPLLSVYVQAANLHLSVLRDVSVFGQRWGFDAATINSRYNDLTRLIGNYTDYAVRWYN'
            'TGLERVWGPDSRDWVRYNQFRRELTLTVLDIVALFSNYDSRRYPIRTVSQLTREIYTNPV'
            'LENFDGSFRGMAQRIEQNIRQPHLMDILNSITIYTDVHRGFNYWSGHQITASPVGFSGPE'
            'FAFPLFGNAGNAAPPVLVSLTGLGIFRTLSSPLYRRIILGSGPNNQELFVLDGTEFSFAS'
            'LTTNLPSTIYRQRGTVDSLDVIPPQDNSVPPRAGFSHRLSHVTMLSQAAGAVYTLRAPTF'
            'SWQHRSAEFNNIIPSSQITQIPLTKSTNLGSGTSVVKGPGFTGGDILRRTSPGQISTLRV'
            'NITAPLSQRYRVRIRYASTTNLQFHTSIDGRPINQGNFSATMSSGSNLQSGSFRTVGFTT'
            'PFNFSNGSSVFTLSAHVFNSGNEVYIDRIEFVPAEVT')
rec = str('VGLPYTVVIHYAGNLSETFHGFYKSTYRTKEGELRILASTQFEPTAARMAFPCFDEPAFK')

print(len(tox), len(rec))
with open('/home/lavrentydanilov/Documents/Documents/Scientific_work/Additional_projects/BioHack_2019/data/parse_mutantion/p_res.txt', 'r') as res:
    for line in res:
        if line.startswith('0'):
            a = line.split(' ')
            subst = a[-2]
            g_energy = a[-1]
            list_sub = subst.split(',')
            rec_mut = rec
            tox_mut = tox
            for el in list_sub:
                sub = el[-1]
                pos = int(el[1:-1])
                if pos > 577:
                    new_pos = int(pos - 578)
                    if new_pos != 0:
                        rec_mut = rec_mut[:new_pos] + sub + rec_mut[new_pos+1:]
                    else:
                        rec_mut = sub + rec_mut[1:]
                else:
                    new_pos_1 = int(pos - 1)
                    tox_mut = tox_mut[:new_pos_1] + sub + tox_mut[new_pos_1+1:]
            st_w = str(tox_mut) + str('\t') + str(rec_mut) + str('\t') + str(g_energy) + str('\n')
            print(len(tox_mut), len(rec_mut), g_energy)
            machine_result.write(st_w)
        else:
            continue
machine_result.close()
