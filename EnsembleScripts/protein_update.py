'''
scan the psm table
find the peptides in the forward protein database
add those proteins to psm entries
Created on Aug. 21, 2017

@author: Xuan Guo
'''

import sys, getopt

# # Version control
def get_version():
    return "protein_update 1.0"

# # Help message
help_message = '''
Usage:
    python protein_update.py [options]

Inputs:
    -d protein database file
    -t psm table file
    -p common name of all target forward proteins

Options:
    -h show help info
    -v show version info

Outputs:
    -o new psm table file
'''

# # Parse options
def parse_options(argv):

    opts, _args = getopt.getopt(argv[1:], 'hvd:t:o:p:')

    # Default working dir and config file
    database_file = ""
    psm_table_file = ""
    new_psm_table_file = ""
    common_name_str = ""

    # Basic options
    for option, value in opts:
        if option in ("-h"):
            print(help_message)
            sys.exit(0)
        elif option in ("-v"):
            print(get_version())
            sys.exit(0)
        elif option in ("-d"):
            database_file = value
        elif option in ("-t"):
            psm_table_file = value
        elif option in ("-o"):
            new_psm_table_file = value
        elif option in ("-p"):
            common_name_str = value
            
    if database_file == "" or psm_table_file == "" or new_psm_table_file == "" or common_name_str == "":
        print(help_message)
        sys.exit(0)

    return (database_file, psm_table_file, new_psm_table_file, common_name_str)

def get_database_into_sequence(database_str):
    database_list = []
    with open(database_str, 'r') as fr:
        for line_str in fr:
            if not line_str.startswith('>'):
                database_list.append(line_str.strip())
            else:
                split_list = line_str.split(' ')
                database_list.append(split_list[0]+',')
    big_database_str = ''.join(database_list)
    return big_database_str

def update_psm_table(big_database_str, psm_table_file, new_psm_table_file, common_name_str):
    psm_table_f = open(psm_table_file, 'r')
    new_psm_table_f = open(new_psm_table_file, 'w')
    
    line_str = psm_table_f.readline()
    split_list = line_str.split('\t')
    protein_name_index = split_list.index('ProteinNames')
    original_peptide_index = split_list.index('OriginalPeptide')
    new_psm_table_f.write(line_str)
    
    for line_str in psm_table_f:
        split_list = line_str.split('\t')
        protein_entry = split_list[protein_name_index]
        if common_name_str not in protein_entry:
            original_peptide_str= split_list[original_peptide_index][1:-1]
            pos = 0
            beg = 0
            protein_names_list = []
            if original_peptide_str in big_database_str:
                while True:
                    pos = big_database_str.find(original_peptide_str, beg)
                    if pos == -1:
                        break
                    st = big_database_str.rfind('>', pos)
                    ed = big_database_str.rfind(',', pos)
                    protein_name = big_database_str[st+1:ed]
                    protein_names_list.append(protein_name)
                    beg = pos + 1
                
                protein_list = protein_entry[1:-1].split(',')
                protein_list.extend(protein_names_list)
                protein_entry = '{' + ','.join(protein_list) + '}'
                split_list[protein_name_index] = protein_entry
                new_psm_table_f.write('\t'.join(split_list))
                continue
        new_psm_table_f.write(line_str)
    
    psm_table_f.close()
    new_psm_table_f.close()

def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse options
    (database_file, psm_table_file, new_psm_table_file, common_name_str) = parse_options(argv)
    
    # get the protein database
    big_database_str = get_database_into_sequence(database_file)
    
    # update psm tabel
    update_psm_table(big_database_str, psm_table_file, new_psm_table_file, common_name_str)
    
    print('Done.')
    
if __name__ == '__main__':
    sys.exit(main())
