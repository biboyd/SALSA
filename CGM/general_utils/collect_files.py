from sys import argv
from os import listdir
from astropy.table import QTable, vstack


def collect_files(directory, file_ext='.h5', key_words=[], black_list=[]):
    """
    Finds and returns files in a given directory with given extension. Optional
    key_word and black_list requirements
    """

    all_files = listdir(directory)

    # get np files only
    files=[]
    for f in all_files:
        #check has file extension and not black listed
        if file_ext in f and check_file(f, key_words, black_list):
            files.append(f)

    return files
def check_file(file, key_words, black_list):
    """
    Check the file against a black list as well as check if it has keywords in
    it.

    Parameters:
        file : string : filename to check
        key_words: list : list of strings required to be in file name
        black_list :list : list of files that are black listed, that file can't be called
    Returns:
        bool : Returns True if file is not in black_list and contains all key_words
    """
    #check if ile in black list
    if file in black_list:
        return False
    else:
        #check if keyword not in file
        for k in key_words:
            if k not in file:
                return False
        #return true if passes through
        return True

def combine_astropy_files(directory, kw='ice', outfile=None):

    #get files
    files = collect_files(directory, key_words=['ray', kw])

    tables = []
    # open up tables
    for f in files:
        tables.append(QTable.read(f"{directory}/{f}"))

    #combine tables
    main_table = vstack(tables)

    #write table
    if outfile is not None:
        main_table.write(outfile, overwrite=True)
    return main_table

if __name__ == '__main__':
    directory= argv[1] 
    outfile= argv[2]  

    combine_astropy_files(directory, outfile=outfile)
