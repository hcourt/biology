#!/usr/bin/python

import sys

#A program to calculate a best alignment for two given ACTG sequences.

## a class which represents each position in the sequence matrix which will be calculated.  Position has a score, an up pointer (0 or 1), a left pointer (0 or 1), and a diagonal pointer (0 or 1).  Each pointer is 0 if the pointer is unused or 1 if the pointer is used.
class Position:
    def __init__(self, score=0):
        self.score=score
        self.up=0
        self.left=0
        self.diag=0

## a function to construct a matrix of length c and width r
def construct_matrix(r, c, type):
    matrix = []
    for i in range(0,r):
        list = []
        for j in range(0,c):
            if type=='p':
                list.append(Position())
            elif type=='s':
                list.append(0)
        matrix.append(list)
    return matrix

## a function which recieves an index for a sequence and returns the alphabet index of the letter
def get_index(seq_index, seq):
    letter = ''
    if seq == 'A':
        letter = left_side[seq_index]
    elif seq == 'B':
        letter = up_side[seq_index]
    for i in range (0, len(alphabet)):
        if alphabet[i]==letter:
            return i
    return -1

## a function which recieves three scores (diagonal, left, and up) and returns a list containing the value of the best score and the codes for the highest score directions
def get_best_point(up, left, diag):
    score = max([diag,left,up])
    return_list = []
    return_list.append(score)
    if (diag==score):
        return_list.append('D')
    if (left==score):
        return_list.append('L')
    if (up==score):
        return_list.append('U')
    return return_list

## the traceback function for all types of alignments
def traceback (cur_r, cur_c, seq_matrix, prev_align_A, prev_align_B, aligns, left_side, up_side):
    pos = seq_matrix[cur_r][cur_c]
    ## end case for 0,0
    if ((cur_r == 0 and cur_c == 0) or (align_type==2 and (cur_r==0 or cur_c==0))):
        aligns['A'].append(prev_align_A[::-1])
        aligns['B'].append(prev_align_B[::-1])
        return aligns
    ## resolve diagonal tree
    if (pos.diag == 1):
        working_align_A = prev_align_A
        working_align_B = prev_align_B
        working_align_A += left_side[cur_r]
        working_align_B += up_side[cur_c]
        aligns=traceback (cur_r-1, cur_c-1, seq_matrix, working_align_A, working_align_B, aligns, left_side, up_side)
    ## resolve left tree
    elif (pos.left == 1):
        working_align_A = prev_align_A
        working_align_B = prev_align_B
        working_align_A+='-'
        working_align_B+=up_side[cur_c]
        aligns=traceback (cur_r, cur_c-1, seq_matrix, working_align_A, working_align_B, aligns, left_side, up_side)
    ## resolve up tree
    elif (pos.up == 1):
        working_align_A= prev_align_A
        working_align_B= prev_align_B
        working_align_A+=left_side[cur_r]
        working_align_B+='-'
        aligns=traceback (cur_r-1, cur_c, seq_matrix, working_align_A, working_align_B, aligns, left_side, up_side)
    return aligns

## MAIN FUNCTION ##

file = open(sys.argv[1])
lines = file.readlines()
for l in range(0,len(lines)):
    lines[l]=lines[l].rstrip()
## sequence of letters indicating sequence A
seqA = lines[1]
left_side = "-"+seqA
## sequence of letters indicating sequence B
seqB = lines[0]
up_side = "-"+seqB
## an argument for the type of alignment to perform (0 is global, 1 is semi-global, 2 is local)
align_type = int(lines[2])
## gap penalty (asumming a constant gap penalty for both seqs)
gap_penalty = int(lines[3])
## sequence of letters indicating the alphabet for the strings
alphabet = lines[4]
n = len(alphabet)
## series of n lines, each consisting of n space-separated values, indicating the score of matching characters between the two sequences in matrix form.  Characters are in the same order as alphabet given above.  Matrix is symmetric.
score_matrix = construct_matrix(n, n, 's')
for i in range (0, n):
    line_list = lines[i+5].split()
    for j in range (0, n):
        score_matrix[i][j] = int(line_list[j])

best_score=0
align_As = []
align_Bs = []
## the matrix into which the sequence scores will be outputted.  Assuming seqB is on top and seqA is to left
seq_matrix = construct_matrix(len(left_side), len(up_side), 'p')
seq_matrix[0][0].score=0

if align_type == 0:
    ## fill first column and first row with gap penalty
    for i in range (1, len(up_side)):
        seq_matrix[0][i].score = gap_penalty + seq_matrix[0][i-1].score
        seq_matrix[0][i].left = 1
    for j in range (1, len(left_side)):
        seq_matrix[j][0].score = gap_penalty + seq_matrix[j-1][0].score
        seq_matrix[j][0].up = 1
elif align_type == 1 or align_type == 2:
    ## fill first column and first row with 0
    for i in range (1, len(up_side)):
        seq_matrix[0][i].score = 0
        seq_matrix[0][i].left = 1
    for j in range (1, len(left_side)):
        seq_matrix[j][0].score = 0
        seq_matrix[j][0].up = 1
## work through each row, calculating the resulting score for the position
for row in range (1, len(left_side)):
    for col in range (1, len(up_side)):
        p = seq_matrix[row][col]
        ## get the predicted scores for each pointer
        p_up_score = 0
        p_left_score = 0
        if align_type == 1:
            ##semi-global has no penalty for end gaps
            p_up_score = seq_matrix[row - 1][col].score if (col == len(up_side)-1) else seq_matrix[row - 1][col].score + gap_penalty
            p_left_score = seq_matrix[row][col - 1].score if (row == len(left_side)-1) else seq_matrix[row][col - 1].score + gap_penalty
        else:
            p_up_score = seq_matrix[row - 1][col].score + gap_penalty
            p_left_score = seq_matrix[row][col - 1].score + gap_penalty
        indexA = get_index(row, 'A')
        indexB = get_index(col, 'B')
        p_diag_score = seq_matrix[row - 1][col - 1].score + score_matrix[get_index(row, 'A')][get_index(col, 'B')]
        if indexA == -1 or indexB == -1:
            print "ERROR: letter combination not in score matrix"
        ## set any negative scores to 0, if doing local alignment
        if align_type == 2:
            if p_up_score < 0: p_up_score = 0
            if p_left_score < 0: p_left_score = 0
            if p_diag_score < 0: p_diag_score = 0
        ## get the best, and assign pointers appropriately
        best = get_best_point(p_up_score, p_left_score, p_diag_score)
        for i in range (1, len(best)):
            if (best[i]=='D'):
                p.diag=1
            elif (best[i]=='L'):
                p.left=1
            elif (best[i]=='U'):
                p.up=1
        p.score=int(best[0])
        seq_matrix[row][col] = p
## tracebacks for global and semi-global alignment
if align_type < 2:
    if len(up_side)>=len(left_side):
        end_score_list=[]
        best_score = seq_matrix[0][len(up_side)-1].score
        ## get best score and all positions to start from
        for i in range (0,len(left_side)):
            pos = seq_matrix[i][len(up_side)-1]
            s = pos.score
            if s > best_score:
                best_score = s
                end_score_list=[]
                end_score_list.append(i)
            elif s == best_score:
                end_score_list.append(i)
        ## begin tracebacks from each start position
        for start in end_score_list:
             ## start traceback call
             aligns={'A':[], 'B':[]}
             aligns=traceback(start, len(up_side) - 1, seq_matrix, "", "", aligns, left_side, up_side)
             align_As.append(aligns['A'])
             align_Bs.append(aligns['B'])
    else:
        end_score_list=[]
        best_score = seq_matrix[len(left_side) - 1][0].score
        ## get best score and all positions to start from
        for i in range (0,len(up_side)):
            pos = seq_matrix[len(left_side) - 1][i]
            s = pos.score
            if s > best_score:
                best_score = s
                end_score_list=[]
                end_score_list.append(i)
            elif s == best_score:
                end_score_list.append(i)
        ## begin tracebacks from each start position
        for start in end_score_list:
            ## start traceback call
            aligns={'A':[], 'B':[]}
            aligns=traceback(len(left_side) - 1, start, seq_matrix, "", "", aligns, left_side, up_side)
            align_As.append(aligns['A'])
            align_Bs.append(aligns['B'])    
##tracebacks for local alignment
elif align_type == 2:
    score_list=[]
    best_score = seq_matrix[0][0].score
    ## get best score and all associated positions to start from
    for row in range(0, len(left_side)):
        for col in range(0, len(up_side)):
            pos = seq_matrix[row][col]
            s = pos.score
            if s > best_score:
                best_score = s
                score_list=[]
                score_list.append((row,col))
            elif s == best_score:
                score_list.append((row,col))
    ## select the furthest to the right, and the furthest down
    furthest=(0,0)
    for (x,y) in score_list:
        if x>=furthest[0] and y>=furthest[1]:
            furthest=(x,y)
    ## traceback call
    aligns={'A':[], 'B':[]}
    aligns = traceback(furthest[0],furthest[1], seq_matrix, "", "", aligns, left_side, up_side)
    align_As.append(aligns['A'])
    align_Bs.append(aligns['B'])
## output, via print:
##     1. score of the best alignment
##     2. alignment of seqA
##     3. alignment of seqB

## toprint seq_matrix, remove both """
"""
sys.stdout.write(" ")
for i in range (0, len(up_side)):
    sys.stdout.write(" "+str(up_side[i]))
print ""
for row in range (0, len(left_side)):
    sys.stdout.write(left_side[row]+" ")
    for col in range(0, len(up_side)):
        sys.stdout.write(str(seq_matrix[row][col].score)+" ")
    print 
"""
    
print("Best score for all alignments: "+str(best_score))
for i in range (0, len(align_As)):
    print(align_Bs[i][0])
    print(align_As[i][0])
