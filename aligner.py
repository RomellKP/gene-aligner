import numpy as np
import argparse

#scores for match. mismatch, and gap


#reads and returns data array and data header
def read_file(input_filename):
    #opens file and reads
    f = open(input_filename)
    data = f.readlines()
    f.close()
    #gets header and data
    data_header = data[0]
    data = ''.join([line.strip() for line in data[1:] if not line.startswith(">")])
    print("read")
    return data, data_header

#outputs data to output file in correct format
def write_file(output, al_score, seq1_head, al1, alvis, al2, seq2_head):
    f = open(output, "a")
    f.write(str(int(al_score)) + "\n")
    f.write(seq1_head)
    f.write(al1 + "\n")
    f.write(alvis + "\n")
    f.write(al2 + "\n")
    f.write(seq2_head)
    f.close()
    print("written")

#the algorithm
def needleman_wunsch(sequence1, sequence2, match_value, gap_penalty, ignore_outer_gaps, mismatch_penalty, output_filename):
    print("started")
    #chcek if init_gap_pen
    init_gap_pen = 0 if ignore_outer_gaps else -2
    #penalties
    gap_pen = gap_penalty
    mismatch_pen = mismatch_penalty
    match_bon = match_value
    #data and data header
    seq1, seq1_h = read_file(sequence1)
    seq2, seq2_h  = read_file(sequence2)
    #lengths of data to access easily 
    len1 = len(seq1)
    len2 = len(seq2)
    #initialize alignment score
    al_score = 0

    #fills matrices with zeros
    #matrix is the matrix with values and ref_matrix is a matrix to determine matches and mismatches
    matrix = np.zeros((len1+1,len2+1))
    matrix_ref = np.zeros((len1,len2))

    #fills in reference matrix with either mismatch_pen or match_bon so we can reference it to matrix
    for i in range(len1):
        for j in range(len2):
            if seq1[i] == seq2[j]:
                matrix_ref[i][j] = match_bon
            else:
                matrix_ref[i][j] = mismatch_pen

    #fills in all of the gap penalties depending on if init_gap_pen is 0 or -2
    for i in range(len1):
        matrix[i][0] = i*init_gap_pen
    for i in range(len2):
        matrix[0][i] = i*init_gap_pen

    #matrix filling
    #takes the max value of the upper+gap_pen, left+gap_pen, and diagonal+corrosponding number in matrix_ref
    for i in range(1, len1+1):
        for j in range(1, len2+1):
            matrix[i][j] = max(matrix[i-1][j-1]+matrix_ref[i-1][j-1], matrix[i][j-1]+gap_pen, matrix[i-1][j]+gap_pen)
            #when the matrix finishes, the alignment score will be matrix[i][j]
            al_score = matrix[i][j]

    #alignment and alignment visualization strings
    al1 = ""
    al2 = ""
    alvis = ""
    #editable lengths for traceback
    edlen1 = len1
    edlen2 = len2
    #traceback
    while edlen1>0 and edlen2>0:
        #checks to see if it is a diagonal traceback
        if (matrix[edlen1][edlen2] == matrix[edlen1-1][edlen2-1]+matrix_ref[edlen1-1][edlen2-1]):
            al1 = seq1[edlen1-1] + al1
            al2 = seq2[edlen2-1] + al2
            #chceks to see if they the same letter and edits visualization accordingly
            if seq1[edlen1-1] == seq2[edlen2-1]:
                alvis = "|" + alvis
            else:
                alvis = "x" + alvis
            edlen1 -= 1
            edlen2 -= 1
        #checks to see if traceback should go left on the row
        elif (matrix[edlen1][edlen2] == matrix[edlen1-1][edlen2]+gap_pen):
            al1 = seq1[edlen1-1] + al1
            al2 = "_" + al2
            alvis = " " + alvis
            edlen1 -= 1
            
        #if the traceback doesnt go diagonal or left it goes up
        elif (matrix[edlen1][edlen2] == matrix[edlen1][edlen2-1]+gap_pen):
            al1 = "_" + al1
            al2 = seq2[edlen2-1] + al2
            alvis = " " + alvis
            edlen2 -= 1

    #writes the extra charcters to al1 and al2 because they are different sizes
    while edlen1>0:
        al1 = seq1[edlen1-1] + al1
        al2 = "_" + al2
        alvis = " " + alvis
        edlen1 -= 1
    while edlen2>0:
        al1 = "_" + al1
        al2 = seq2[edlen2-1] + al2
        alvis = " " + alvis
        edlen2 -= 1

    #calls write_file() to write the alignment
    write_file(output_filename, al_score, seq1_h, al1, alvis, al2, seq2_h)

    #returns alignment 1 and 2 and alignment score
    return(al1, al2, al_score)
    

def main():
    #parser stuff
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--query', required=True, type=str, dest = 'seq1_file')
    parser.add_argument('-r', '--reference', required=True, type=str, dest = 'seq2_file')
    parser.add_argument('-o', '--output', required=True, type=str, dest = 'output')
    parser.add_argument('-g', '--gap_penalty', required=True, type=int, default = -2)
    parser.add_argument('-p', '--mismatch_penalty', required=True, type=int, default = -1)
    parser.add_argument('-m', '--match_score', required=True, type=int, default = 1)
    parser.add_argument('--ignore_outer_gaps', action='store_true')
    parser.add_argument('-s', '--affine_gap_penalty', required=False, default=0, type=int)
    args = parser.parse_args()
    #algorithm function call
    needleman_wunsch(args.seq1_file, args.seq2_file, args.match_score, args.gap_penalty, args.ignore_outer_gaps, args.mismatch_penalty, args.output)

if __name__ == "__main__":
    main()