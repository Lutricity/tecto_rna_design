# input an RNA sequence and outputs the molecular weight in Da (g/mol)

#RNA_seq = input("write the sequence of your RNA ")


#Only use global variables with constants (usually)
A_WEIGHT = 507.2
U_WEIGHT = 484.2
G_WEIGHT = 523.2
C_WEIGHT = 483.2

#this function takes the RNA sequences and returns the number of each nucleotide present.
def count_nucleotides(RNA_seq):
    x = len(RNA_seq)

    pos = 0
    A_count = 0
    U_count = 0
    G_count = 0
    C_count = 0
    for i in range(x):
        if pos <= len(RNA_seq):
            if RNA_seq[pos] == "A":
                A_count = A_count + 1
                pos += 1
            elif RNA_seq[pos] == "U":
                U_count = U_count + 1
                pos += 1
            elif RNA_seq[pos] == "G":
                G_count = G_count + 1
                pos += 1
            elif RNA_seq[pos] == "C":
                C_count = C_count + 1
                pos += 1

        elif pos == len(RNA_seq):
            if RNA_seq[pos] == "A":
                A_count = A_count + 1
                pos += 1
            elif RNA_seq[pos] == "U":
                U_count = U_count + 1
                pos += 1
            elif RNA_seq[pos] == "G":
                G_count = G_count + 1
                pos += 1
            elif RNA_seq[pos] == "C":
                C_count = C_count + 1
                pos += 1

    return(A_count, U_count, G_count, C_count)


#This function takes the number of each nucleotide present and adds together the weight of all these in Daltons.
def calc_RNA_Daltons(A_count, U_count, G_count, C_count):
    global A_WEIGHT
    global U_WEIGHT
    global G_WEIGHT
    global C_WEIGHT

    output_A = A_count * A_WEIGHT
    output_U = U_count * U_WEIGHT
    output_G = G_count * G_WEIGHT
    output_C = C_count * C_WEIGHT

    final_seq_weight = output_A + output_U + output_G + output_C

    return final_seq_weight

#This prints the final weight of the sequence with a descriptor.
def print_final_seq_weight(final_seq_weight):
    print("This is your RNA sequence weight: ", final_seq_weight)

#My main function to string all these together and return the value I want (weight of RNA seq in Daltons)
def main():
    RNA_seq = input("write the sequence of your RNA: ")


    pass_counts = count_nucleotides(RNA_seq.strip())
    print(pass_counts)

    #print (A_count, U_count, G_count, C_count)

    final_weight =  calc_RNA_Daltons(*pass_counts)


    print_final_seq_weight(final_weight)

if __name__ == "__main__":
     main()

