def mask():
    """This function takes the sequence and masks it"""
    with open("ITS2_collapsed_custom_insilico.fasta","r") as file:
        seq=[]
        for line in file:
            if line.startswith('>'):
                seq.append(line)
            else:
                length= len(line)
                if length< 150:
                    seq.append(line)
                elif length> 150:
                    line = line[:150]+ (20*"N")
                    seq.append(line)
    for line in seq:
        print(line)

if __name__ == "__main__":
    mask()
