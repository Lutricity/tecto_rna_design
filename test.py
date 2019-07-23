import pandas as pd

# Make three cutting/tiling programs- where the internal
# region shortens from the top, from the bottom,
# and depending on the even and odd cuts up one way, then down from the other


# I'll have to define the hairpin region variably? in each cut program?
"""def seq_1cut_bottom(long_seq, long_ss):
    count = 0

    start_pos = 11
    end_pos = -13
    end_3variable = -21
    end_5variable = end_3variable -7 #(position 20 on varna)


    while 1:

    start_5prime_cut = long_seq[0:start_pos]
    ss_start_5prime_cut = long_ss[0:start_pos]
    start_3prime_cut = long_seq[-1:end_pos]
    ss_start_3prime_cut = long_seq[-1:end_pos]
    variabreg_5prime = long_seq[start_pos:]
    variable_3prime = long_seq[end_3variable:]

    start_pos += 1
    end_pos -= 1
    end_3variable

    if prime5_helix < 1:
        exit()
"""


# every time make end_5variable smaller ex: -20, -19...

# def seq_2cut_top(long_seq, long_ss):

# def seq_2cut_bottom(long_seq, long_ss):

def flip_motif(seq_1, seq_2, ss_1, ss_2):
    flip_seq1 = seq_2
    flip_seq2 = seq_1
    flip_ss1 = ""
    flip_ss2 = ""

    for i in range(0, len(ss_2)):
        if ss_2[i] == ")":
            flip_ss1 += "("
        else:
            flip_ss1 += "."
    #print(flip_ss1)

    for i in range(0,len(ss_1)):
        if ss_1[i] == "(":
            flip_ss2 += ")"
        else:
            flip_ss2 += "."
    #print(flip_ss2)

    return [flip_seq1, flip_seq2, flip_ss1, flip_ss2]


# Next step- work flip script into the order of operations to it runs quickly and seemlessly

def tile_motif_in_chip_sequence(chip_seq, chip_ss, seq_1, seq_2, ss_1, ss_2):
    count = 0
    # the first while loop i'll have to make should contain flip scrip and while the count = 0 then no flip script, else flip script, exit at third count
    # make another while loop up here that edits the original sequence between certain positions, cutting in pairs or adding in pairs

    # elif len(chip_seq)
    start_pos = 11
    end_pos = -12

    length_offset = (len(chip_seq) - 47) / 2

    start_hair = 20 + length_offset
    end_hair = 26 + length_offset
    # print (length_offset)

    tecto_variants = []

    while 1:

        start_5prime = chip_seq[0:start_pos]
        ss_start_5prime = chip_ss[0:start_pos]
        start_3prime = chip_seq[end_pos:]
        ss_start_3prime = chip_ss[end_pos:]
        hairpin = chip_seq[start_hair:end_hair]
        ss_hairpin = chip_ss[start_hair:end_hair]
        helix_5prime = chip_seq[start_pos:start_hair]
        ss_helix_5prime = chip_ss[start_pos:start_hair]
        helix_3prime = chip_seq[end_hair:end_pos]
        ss_helix_3prime = chip_ss[end_hair:end_pos]
        motif_length = len(seq_1)
        if motif_length < len(seq_2):
            motif_length = len(seq_2)

        if len(helix_5prime) < motif_length:
            break

        # rint start_5prime, start_3prime

        new_helix_5prime = seq_1 + helix_5prime[motif_length:start_hair + 1]
        new_helix_3prime = helix_3prime[0:-motif_length] + seq_2

        ss_new_helix_5prime = ss_1 + ss_helix_5prime[motif_length:start_hair + 1]
        ss_new_helix_3prime = ss_helix_3prime[0:-motif_length] + ss_2

        final_sequence = start_5prime + new_helix_5prime + hairpin + new_helix_3prime + start_3prime
        ss_final = ss_start_5prime + ss_new_helix_5prime + ss_hairpin + ss_new_helix_3prime + ss_start_3prime

        #print(final_sequence)
        #print(ss_final)

        tecto_variants.append([final_sequence, ss_final])




        start_pos += 1
        end_pos -= 1
        count += 1
        # if count > 1:
        #    exit()
        # exit()
    return tecto_variants


def main():

    # The idea is that if you give the program the longest number RNA for it to cut down then
    # It will sequentially tile it for you. Make a version of this program where you pair down from the top,
    # and make a version where you pair down from the bottom
    tar_seq_1 = "GAUCUG"
    tar_seq_2 = "CUC"

    tar_ss_1 = "((...("
    tar_ss_2 = ")))"

    df = pd.read_csv("chip_pieces.csv")
    df_result = pd.DataFrame(columns='project_name chip_type chip_sequence chip_structure insert_sequence insert_structure'.split())

    #print(df)
    pos = 0
    for i, row in df.iterrows():
        motif_tile_output = tile_motif_in_chip_sequence(row.sequence, row.structure, tar_seq_1, tar_seq_2, tar_ss_1, tar_ss_2)
        for l in motif_tile_output:
            df_result.loc[pos] = ["project_TAR", row['name'], l[0], l[1], tar_seq_1+"+"+tar_seq_2, tar_ss_1+"+"+tar_ss_2]
            pos += 1


    df_result.to_csv("motif_sequences.csv",index=False)
    #print(df_result)
    exit()


   # f = open("chip_pieces.csv")
    #lines = f.readlines()
    #f.close()

    #lines.pop(0)

    #for l in lines:
     #   l_strip = l.rstrip()
      #  spl = l_strip.split(",")
        #name, sequence, structure = spl
       # name = spl[0]
        #sequence = spl[1]
        #structure = spl[2]
        #print(name, sequence, structure)


        #name, sequence, structure = l.rstrip().split(",")
    #print sequence, structure

     #   final_sequences.append(tile_motif_in_chip_sequence(sequence, structure, tar_seq_1, tar_seq_2, tar_ss_1, tar_ss_2))

    #print(tile_motif_in_chip_sequence(seq_2, ss_2, tar_seq_1, tar_seq_2, tar_ss_1, tar_ss_2))
    # pum2 = "MNHDFQALALESRGMGELLPTKKFWEPDDSTKDGQKGIFLGDDEWRETAWGASHHSMSQPIMVQRRSGQGFHGNSEVNAILSPRSESGGLGVSMVEYVLSSSPADKLDSRFRKGNFGTRDAETDGPEKGDQKGKASPFEEDQNRDLKQGDDDDSKINGRGLPNGMDADCKDFNRTPGSRQASPTEVVERLGPNTNPSEGLGPLPNPTANKPLVEEFSNPETQNLDAMEQVGLESLQFDYPGNQVPMDSSGATVGLFDYNSQQQLFQRTNALTVQQLTAAQQQQYALAAAQQPHIAGVFSAGLAPAAFVPNPYIISAAPPGTDPYTAAGLAAAATLAGPAVVPPQYYGVPWGVYPANLFQQQAAAAANNTASQQAASQAQPGQQQVLRAGAGQRPLTPNQGQQGQQAESLAAAAAANPTLAFGQGLATGMPGYQVLAPTAYYDQTGALVVGPGARTGLGAPVRLMAPTPVLISSAAAQAAAAAAAGGTASSLTGSTNGLFRPIGTQPPQQQQQQPSTNLQSNSFYGSSSLTNSSQSSSLFSHGPGQPGSTSLGFGSGNSLGAAIGSALSGFGSSVGSSASSSATRRESLSTSSDLYKRSSSSLAPIGQPFYNSLGFSSSPSPIGMPLPSQTPGHSLTPPPSLSSHGSSSSLHLGGLTNGSGRYISAAPGAEAKYRSASSTSSLFSSSSQLFPPSRLRYNRSDIMPSGRSRLLEDFRNNRFPNLQLRDLIGHIVEFSQDQHGSRFIQQKLERATPAERQMVFNEILQAAYQLMTDVFGNYVIQKFFEFGSLDQKLALATRIRGHVLPLALQMYGCRVIQKALESISSDQQVISEMVKELDGHVLKCVKDQNGNHVVQKCIECVQPQSLQFIIDAFKGQVFVLSTHPYGCRVIQRILEHCTAEQTLPILEELHQHTEQLVQDQYGNYVIQHVLEHGRPEDKSKIVSEIRGKVLALSQHKFASNVVEKCVTHASRAERALLIDEVCCQNDGPHSALYTMMKDQYANYVVQKMIDMAEPAQRKIIMHKIRPHITTLRKYTYGKHILAKLEKYYLKNSPDLGPIGGPPNGML"
    # print len(pum2)

    # print tar_seq_1, tar_seq_2
    #flipped_motif= flip_motif(tar_seq_1, tar_seq_2, tar_ss_1, tar_ss_2)
    #print flipped_motif
    #tile_motif_in_chip_sequence(seq_2, ss_2, flipped_motif[0], flipped_motif[1], flipped_motif[2], flipped_motif[3])



if __name__ == "__main__":
    main()

#final_sequences.to_csv("motif_sequences.csv")

