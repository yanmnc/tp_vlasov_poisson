import os

from create_file_name import create_file_name

def create_file_list(acronym, extension, syntax):
    ifile = 0
    filename = create_file_name(acronym, 0, extension, syntax)
    nb_files = 0
    file_list = []
    while os.path.exists(filename):
        ifile += 1
        nb_files += 1
        file_list.append(filename)
        filename = create_file_name(acronym, ifile, extension, syntax)
    return file_list, nb_files