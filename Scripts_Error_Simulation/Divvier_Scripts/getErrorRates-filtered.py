import sys


"""
Checks if line is the beginning of an alignment

Args: 
    line (str): the line to check

Return:
    bool: whether the string is the beginning of an alignment
"""
def isBeginAlignment( line ):
    if len( line.split(">") ) > 1:
        return True
    return False


"""
Gets the FP, FN, TP, and TN of the current line

Args:
    correct (str): the line from the correct file
    error   (str): the line from the error file
    output  (str): the line from the output file

Return:
    tuple: tuple of (fp, fn, tp, tn)

"""
def getLineErrorRates( correct, error, output ):
    fp = 0
    fn = 0
    tp = 0
    tn = 0
    for x in range( len(correct) ):
        output_char = output[x]
        error_char = error[x]
        correct_char = correct[x]
        if correct_char == '-':
            continue
        if correct_char == '.' and output_char != error_char:
            tp += 1
        elif correct_char == '.' and output_char == error_char:
            fn += 1
        elif correct_char != '.' and output_char != error_char:
            fp += 1
#            print("expected " + error_char + " but got " + output_char)
        else:
            tn += 1
            
    return (fp, fn, tp, tn)
        


USAGE = "python getErrorRates.py [correct file] [error file] [correction output file] [description]"

if not len(sys.argv) == 5:
    print()
    print("\tIncorrect Parameters")
    print("\tUSAGE: " + USAGE)
    sys.exit()

threshold = 0.05

correct = sys.argv[1]
error = sys.argv[2]
output = sys.argv[3]
description = sys.argv[4]

correct_f = open( correct, "r" )
error_f = open( error, "r" )
output_f = open( output, "r" )

false_positive = 0
false_negative = 0
true_positive = 0
true_negative = 0
correct_lines = correct_f.read().split('\n')[1::2]
error_lines = error_f.read().split('\n')[1::2]
output_lines = [list(s) for s in output_f.read().split('\n')[1::2]]
if len(correct_lines) != len(error_lines) or len(correct_lines) != len(output_lines):
    print("Error: Not all files have the same number of lines: " + description)
correct_f.close()
error_f.close()
output_f.close()
for i in range(len(output_lines[0])):
    if len([j for j in range(len(output_lines)) if output_lines[j][i] != '-']) < threshold * len(output_lines):
        for j in range(len(output_lines)):
            output_lines[j][i] = '-'
for i in range(len(output_lines)):
    line_stats = getLineErrorRates( correct_lines[i], error_lines[i], output_lines[i] )
    false_positive += line_stats[0]
    false_negative += line_stats[1]
    true_positive += line_stats[2]
    true_negative += line_stats[3]
print("FP: " + str(false_positive))
print("FN: " + str(false_negative))
print("TP: " + str(true_positive))
print("TN: " + str(true_negative))
#description = description.decode('unicode_escape') #python3
#description = description.decode('string_escape') #python2
sys.stderr.write(str(description) + "\t" + str(false_positive) + "\t" + str(false_negative) + "\t" + str(true_positive) + "\t" + str(true_negative) + "\n")
