'''
Created on May 5, 2017

@author: Naux
'''

import sys, getopt, os
#from sets import Set

## Version control
def get_version():
    return "v1.0 (Alpha)"

## Help message
help_message = '''

Usage:
    python ConstructDatabase.py [options]

Inputs:
    [*.pro.txt] files
        (multiple protein assembling results can be processed)
        (search automatically in current directory)
    [*.pep.txt] files
    [*.psm.txt] files
    fasta protein database file

Options:
    -h/--help
    -v/--version
    -w/--working-directory ./path/       # Directory path containing pro.txt/pep.txt/psm.txt output files
    -d/--database-file database.fasta    # fasta protein database
    -o/--output-file                     # new database file
    -pro                                 # get proteins from *.pro.txt (default)
    -pep                                 # get proteins from *.pep.txt 
    -psm                                 # get proteins from *.psm.txt

Outputs:
    DatabaseName.fasta
'''

pro_int = 1
pep_int = 2
psm_int = 4

## Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], "hvVw:d:o:p:x:",
                                ["help",
                                 "version",
                                 "working-directory",
                                 "database-file",
                                 "output-file",
                                 "input-file",
                                 "xdebug"])

    working_directory = "./"
    database_file = ""
    output_file = ""
    input_file_format_int = 0
    output_protein_file = None

    # Basic options
    for option, value in opts:
        if option in ("-h", "--help"):
            print(help_message)
            sys.exit(0)
        if option in ("-v", "-V", "--version"):
            print("ConstructDatabase.py {}".format(get_version()))
            sys.exit(0)
        if option in ("-w", "--working-directory"):
            working_directory = value
            if working_directory[-1] != '/':
                working_directory = working_directory + '/'
        if option in ("-d", "--database-file"):
            database_file = value
        if option in ("-o", "--output-file"):
            output_file = value
        if option in ("-x", "--xdebug"):
            output_protein_file = value
        if option in ("-p"):
            if value == 'ro':
                input_file_format_int = input_file_format_int | pro_int
            elif value == 'ep':
                input_file_format_int = input_file_format_int | pep_int
            elif value == 'sm':
                input_file_format_int = input_file_format_int | psm_int
            else:
                print(help_message)
                sys.exit(0)
    if input_file_format_int == 0:
        input_file_format_int = pro_int
    
    if working_directory == "" or database_file == "" or output_file == "":
        print(help_message)
        sys.exit(0)

    return (working_directory, database_file, output_file, input_file_format_int, output_protein_file)

reverse_prefix_str = 'Rev'

def get_protein_in_pro_file(filename_str):
    protein_set = set()
    with open(filename_str, 'r') as f:
        for line_str in f:
            if not (line_str.startswith('#') or line_str.startswith('ProteinID')):
                if line_str.startswith('{'):
                    e_list = line_str.split('\t')[0].strip().replace('{', '').replace('}', '').split(',')
                    for v in e_list:
                        if not v.startswith('Rev'):
                            protein_set.add(v)
                else:
                    v = line_str.split('\t')[0].strip()
                    if not v.startswith('Rev'):
                        protein_set.add(v)
    return protein_set

def get_protein_in_pep_file(filename_str):
    protein_set = set()
    with open(filename_str, 'r') as f:
        for line_str in f:
            if not (line_str.startswith('#') or line_str.startswith('IdentifiedPeptide')):
                v = line_str.split('\t')[3].strip()
                if v.startswith('{'):
                    e_list = v.replace('{', '').replace('}', '').split(',')
                    for e in e_list:
                        if not e.startswith('Rev'):
                            protein_set.add(e)
                else:
                    e_list = v.split(',')
                    for e in e_list:
                        if not e.startswith('Rev'):
                            protein_set.add(e)
    return protein_set

def get_protein_in_psm_file(filename_str):
    protein_set = set()
    with open(filename_str, 'r') as f:
        for line_str in f:
            if not (line_str.startswith('#') or line_str.startswith('Filename')):
                v = line_str.split('\t')[15].strip()
                if v.startswith('{'):
                    e_list = v.replace('{', '').replace('}', '').split(',')
                    for e in e_list:
                        if not e.startswith('Rev'):
                            protein_set.add(e)
                else:
                    e_list = v.split(',')
                    for e in e_list:
                        if not e.startswith('Rev'):
                            protein_set.add(e)
    return protein_set

def create_database(protein_set, output_str, protein_database_str):
    
    fw = open(output_str, 'w')
    id_str = ''
    id_full_str = ''
    sequence_str = ''
    save_bool = False
    protein_count = 0
    protein_saved_set = set()
    with open(protein_database_str, 'r') as f:
        for line_str in f:
            if line_str.startswith('>'):
                if save_bool:
                    protein_count += 1
                    fw.write(id_full_str)
                    fw.write(sequence_str)
                    sequence_str = ""
                    fw.write('\n')
                save_bool = False
                id_str = line_str.strip().split(' ')[0][1:]
                id_full_str = line_str
                if id_str in protein_set:
                    save_bool = True
                    protein_saved_set.add(id_str)
                sequence_str = ""
            else:
                sequence_str += line_str.strip()
    
    if save_bool:
        protein_count += 1
        fw.write(id_full_str)
        fw.write(sequence_str)
        fw.write('\n')
    
    print(str(protein_count))
    print(str(len(protein_set)))
    
    fw.close()
    
def create_database_multiple_files(directory_str, output_str, protein_database_str, input_file_format_int, output_protein_id_str=None):
    
    file_list = []
    
    for root, _directories, filenames in os.walk(directory_str):
        for filename in filenames: 
            if filename.endswith('.pro.txt') and (input_file_format_int & pro_int) == pro_int:
                file_list.append(os.path.join(root,filename))
            elif filename.endswith('.pep.txt') and (input_file_format_int & pep_int) == pep_int:
                file_list.append(os.path.join(root,filename))
            elif filename.endswith('.psm.txt') and (input_file_format_int & psm_int) == psm_int:
                file_list.append(os.path.join(root,filename))
    
    for file_str in file_list:
        print(file_str)
    
    protein_set = set()
    for file_str in file_list:
        if file_str.endswith('.pro.txt'):
            protein_local_set = get_protein_in_pro_file(file_str)
            for protein in protein_local_set:
                protein_set.add(protein)
        elif file_str.endswith('.pep.txt'):
            protein_local_set = get_protein_in_pep_file(file_str)
            for protein in protein_local_set:
                protein_set.add(protein)
        elif file_str.endswith('.psm.txt'):
            protein_local_set = get_protein_in_psm_file(file_str)
            for protein in protein_local_set:
                protein_set.add(protein)                
    
    create_database(protein_set, output_str, protein_database_str)
    
    if output_protein_id_str:
        fw = open(output_protein_id_str, 'w')
        for protein in protein_set:
            fw.write(protein)
            fw.write("\n")
        fw.close()
            


# # +------+
# # | Main |
# # +------+
def main(argv=None):
    if argv is None:
        argv = sys.argv
    
    # parse options
    (working_dir, protein_database_str, output_str, input_file_format_int, output_protein_file) = parse_options(argv)
    
    create_database_multiple_files(working_dir, output_str, protein_database_str, input_file_format_int, output_protein_file)
    
    print("Done")

if __name__ == '__main__':
    sys.exit(main())

