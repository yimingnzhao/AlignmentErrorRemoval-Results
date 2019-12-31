import random
import numpy as np
import sys
import os

RNA_DATA = ["A", "U", "G", "C"];
DNA_DATA = ["A", "T", "G", "C"];
#AA_DATA = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
AA_DATA = ["F", "F",
        "L", "L", "L", "L", "L", "L", 
        "I", "I", "I",
        "M",
        "V", "V", "V", "V", 
        "S", "S", "S", "S", "S", "S",
        "P", "P", "P", "P",
        "T", "T", "T", "T",
        "A", "A", "A", "A",
        "Y", "Y",
        "H", "H",
        "Q", "Q",
        "N", "N",
        "K", "K",
        "D", "D",
        "E", "E", 
        "C", "C",
        "W",
        "R", "R", "R", "R", "R", "R",
        "G", "G", "G", "G",
]


"""
Checks if an input string is an integer

Args:
    num (str): the string to check

Return:
    bool: whether the string is an int or not
"""
def isInt( num ):
    try:
        int(num);
        return True;
    except ValueError:
        return False;


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
Generates the alignment to a single line

Args:
    path (str): the path of the file to read
    new_file (str): the file name of the reformatted file

Return:
    void
"""
def reformatFile( path, new_file ):
    f = open( path, "r" );
    result_f = open( new_file, "a" );
    current_line = f.readline();
    while current_line:
        if isBeginAlignment( current_line ):
            result_f.write( current_line );
        else:
            result_f.write( current_line[:-1] );
        current_line = f.readline();
        if isBeginAlignment( current_line ):
            result_f.write( "\n" );

    f.close();
    result_f.close();

    
"""
Uniformly chooses the alignments to contain errors

Args:
    length (int): the number of erroneous alignments
    size (int): the total number of alignments

Return:
    list: list of sequence indices to modify
"""
def getErrSequences( length, size ):
    sequence_errs = []
    for i in range( length ):
        while True:
            rand_seq_idx = random.randint( 1, size );
            if not rand_seq_idx in sequence_errs:
                sequence_errs.append( rand_seq_idx );
                break;
    sequence_errs.sort( reverse=True );
    return sequence_errs


"""
Gets the last possible char in range to randomly choose the start position of an error sequence
Skips over gaps ("-")

Args: 
    sequence (str): the sequence to modify
    length (int): the number of characters to modify

Return:
    int: the last possible index of the sequence to begin the error sequence of a given length
"""
def getLastCharInRange( sequence, length ):
    char_count = 0;
    current_index = len(sequence) - 2;
    while char_count < length:
        # Checks if current_index has become negative and exits program if true
        if current_index < 0:
            return 0
        current_char = sequence[current_index];
        # Only increments char_count if the char is not a gap
        if not current_char == "-":
            char_count += 1;
        current_index -= 1;
    return current_index + 1;



"""
Uniformly chooses the start position of error under the following conditions:
    - given a set length of the error
    - randomly chooses start position
    - randomly chooses the sequence character

Args:
    sequence (str): the sequence to modify
    length (int): the number of characters to modify
    data (list): list of characters

Return:
    tuple: (modified sequence with errors chars, modified sequence with error positions)
"""
def setErrSequence( sequence, length, data ):
    last_char = getLastCharInRange(sequence, length)
    if last_char == -1:
        return (0, 0)
    current_pos = random.randint( 0, last_char);
    count = 0;
    err_sequence = sequence
    pos_sequence = sequence
    while count < length:
        if not sequence[current_pos] == "-":
            count += 1;
            rand_segment = data[ random.randint( 0, len( data ) - 1 ) ];
            err_sequence = err_sequence[0:current_pos] + rand_segment + err_sequence[(current_pos + 1):];
            pos_sequence = pos_sequence[0:current_pos] + "." + err_sequence[(current_pos + 1):];
        current_pos += 1;
    return (err_sequence, pos_sequence);



"""
Checks if the last char is a newline and adds a newline if not

Args:
    path (str): the path to the file

Return:
    void
"""
def addNewlineToEOF( path ):
    with open(path, 'r+') as f:
        f.seek(0, os.SEEK_END)  # go at the end of the file
        f.seek(f.tell() - 1, os.SEEK_SET);
        if f.read(1) != '\n':
            # add missing newline if not already present
            f.write('\n')
            f.flush()
            f.seek(0)


# The dataset file
USAGE = "python generateErrorModel_MultipleErrorAreas.py [data file] [num of erroneous sequences] [value of k] [DNA/RNA/AA]"
DESCRIPTION = "Generates an error alignment with a random error length [2*k, 64*k] for each erroneous sequence. Outputs files:\n\t\treformat.fasta - alignment file with one sequence per line\n\t\terror.fasta - erroneous file\n\t\tposition.fasta - file where each inserted error character is represented by '.'"


if not len(sys.argv) == 5:
    print();
    print("\tError: Incorrect number of parameters\n");
    print("\tUSAGE: " + USAGE );
    print()
    print("DESCRIPTION: " + DESCRIPTION)
    print()
    sys.exit();
data_file = sys.argv[1];
num_erroneous_alignments = sys.argv[2];
k = sys.argv[3];
data_type = sys.argv[4];
if not isInt( num_erroneous_alignments ):
    print();
    print("\tError: Invalid parameter for [num of errneous alignments]");
    print("\tUSAGE: " + USAGE);
    sys.exit();
if not isInt( k ):
    print();
    print("\tError: Invalid parameter for [k]");
    print("\tUSAGE: " + USAGE);
    sys.exit();

num_erroneous_alignments = int(num_erroneous_alignments);
k = int(k)
reformat_file = "reformat.fasta";
error_file= "error.fasta";
position_file = "position.fasta"



# Creates a reformated file
reformatFile( data_file, reformat_file );

# Gets data about reformatted alignment file
num_alignments = 0;
chars_in_alignment = 0;
with open( reformat_file, "r" ) as file_object:
    for line in file_object:
        if line[0] == ">":
            continue;
        if chars_in_alignment == 0:
            chars_in_alignment = len( line );
        num_alignments += 1;
sys.stderr.write("Number of Alignments: " + str(num_alignments) + "\n");
sys.stderr.write("Chars in Alignment: " + str(chars_in_alignment) + "\n");

if ( num_alignments < num_erroneous_alignments ):
    num_erroneous_alignments = num_alignments;

# Creates file with alignment sequence errors
f = open( reformat_file, "r" );
error_f = open( error_file, "a" );
pos_f = open(position_file, "a")
sequence_errs = getErrSequences( num_erroneous_alignments, num_alignments );
count = 0;
if  data_type == "RNA":
    data_type = RNA_DATA
elif data_type == "AA":
    data_type = AA_DATA
else:
    data_type = DNA_DATA
with open( reformat_file, "r" ) as file_object:
    for line in file_object:
        if isBeginAlignment( line ):
            error_f.write(line);
            pos_f.write(line)
            count += 1;
            continue;
        if len(sequence_errs) > 0 and sequence_errs[-1] == count:
            
            #HEREEEEE
            
            sequence_error_len = int(abs(np.random.normal(50, 10)))
            sequences = setErrSequence( line, sequence_error_len, data_type )
            if sequences[0] == 0 and sequences[1] == 0:
                f.close();
                error_f.close();
                pos_f.close();
                open(reformat_file, 'w').close()
                open(error_file, 'w').close()
                open(position_file, 'w').close()
                sys.exit()
            error_f.write( sequences[0] );
            pos_f.write( sequences[1] );
            sequence_errs.pop();
        else:
            error_f.write( line );
            pos_f.write(line)
f.close();
error_f.close();
pos_f.close();
addNewlineToEOF(reformat_file);
addNewlineToEOF(error_file);
addNewlineToEOF(position_file);
