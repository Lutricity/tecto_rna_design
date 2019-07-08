seq = "CTAGGATATGGAAGATCCTCGGGAACGAGGATCTTCCTAAGTCCTAG"
ss  = "(((((((..((((((((((((....))))))))))))...)))))))"


print(seq[-1])
print(seq[-2])

tar_seq_1 = "GAUCUG"
tar_seq_2 = "CUC"

tar_ss_1 = "((...("
tar_ss_2 = ")))"

print(seq[-1])

def tile_motif_in_chip_sequence(chip_seq, chip_ss, seq_1, seq_2, ss_1, ss_2):
    start_pos = 11
    end_pos = 35
    count = 0
    while 1:
        start_5prime = chip_seq[0:start_pos]
        ss_start_5prime = chip_ss[0:start_pos]
        start_3prime = chip_seq[end_pos:]
        ss_start_3prime = chip_ss[end_pos:]
        hairpin = chip_seq[20:26]
        ss_hairpin = chip_ss[20:26]
        helix_5prime = chip_seq[start_pos:20]
        ss_helix_5prime = chip_ss[start_pos:20]
        helix_3prime = chip_seq[26:end_pos]
        ss_helix_3prime = chip_ss[26:end_pos]
        motif_length = len(seq_1)
        if motif_length < len(seq_2):
            motif_length = len(seq_2)
#what's up
    

        if len(helix_5prime) < motif_length:
            break

        new_helix_5prime = seq_1 + helix_5prime[motif_length:21]
        new_helix_3prime = helix_3prime[0:-motif_length] + seq_2


        ss_new_helix_5prime = ss_1 + ss_helix_5prime[motif_length:21]
        ss_new_helix_3prime = ss_helix_3prime[0:-motif_length] + ss_2

        final_sequence = start_5prime + new_helix_5prime + hairpin + new_helix_3prime + start_3prime
        ss_final = ss_start_5prime + ss_new_helix_5prime + ss_hairpin + ss_new_helix_3prime + ss_start_3prime

        print(final_sequence)
        print(ss_final)
        print


        start_pos += 1
        end_pos -= 1
        count += 1
        #if count > 1:
        #    exit()

tile_motif_in_chip_sequence(seq, ss, tar_seq_1, tar_seq_2, tar_ss_1, tar_ss_2)
