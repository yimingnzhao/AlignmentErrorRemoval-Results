import random
import sys


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
Gets the number of sequences in the file

Args:
    path (str): the path of the file to read

Return:
    int: the number of sequences in the file 
"""
def getSequencesInFile( path ):
    count = 0;
    for line in open( path ):
        if ">" in line:
            count += 1;
    return count;

    
"""
Uniformly chooses sequences

Args:
    count (int): the number of sequences to choose
    size (int): the size of the file

Return:
    list: list that specifies which indices should be chosen
"""
def getRandomSequenceIndices( count, size ):
    indices = [];
    for i in range( count ):
        while True:
            rand_index = random.randint( 1, size );
            if not rand_index in indices:
                indices.append( rand_index );
                break;
    indices.sort( reverse=True );
    return indices;



USAGE = "python chooseSeqeunces.py [data file] [num of sequences to choose]"
DESCRIPTION = "Chooses random sequences in an alignment file and generates a new file 'chosen_sequences.fasta' as the output";
if not len(sys.argv) == 3:
    print();
    print("\tError: Incorrect number of parameters\n");
    print("\tUSAGE: " + USAGE);
    print()
    print("\tDESCRIPTION: " + DESCRIPTION)
    print()
    sys.exit();
if not isInt(sys.argv[2]):
    print();
    print("\tError: Invalid parameter for [num of sequences to choose]");
    print("\tUSAGE: " + USAGE);
    print()
    sys.exit();

print("Getting random sequences...")
file_size = getSequencesInFile(sys.argv[1]);


# if sequences to choose > sequences in file, then all sequences are chosen
if file_size < int(sys.argv[2]):
    f = open( sys.argv[1], "r" )
    result_f = open( "chosen_sequences.fasta", "a" )
    line = f.readline()
    while line:
        result_f.write(line)
        line = f.readline()
    f.close()
    result_f.close()
    sys.exit()

indices = getRandomSequenceIndices( int(sys.argv[2]), file_size );
f = open( sys.argv[1], "r" );
result_f = open( "chosen_sequences.fasta", "a" );
count = 1;
while len(indices) > 0:
    data = f.readline() + f.readline();
    if indices[-1] == count:
        result_f.write( data );
        indices.pop();
    count += 1;
f.close();
result_f.close();


    
