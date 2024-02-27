def create_file_name(acronym, ifile, extension, syntax):
    if syntax == 1:
        if ifile < 10:
            file_name = f"{acronym}_d0000{ifile}{extension}"
        elif ifile < 100:
            file_name = f"{acronym}_d000{ifile}{extension}"
        elif ifile < 1000:
            file_name = f"{acronym}_d00{ifile}{extension}"
        elif ifile < 10000:
            file_name = f"{acronym}_d0{ifile}{extension}"
        else:
            file_name = f"{acronym}_d{ifile}{extension}"
    else:
        if ifile < 10:
            file_name = f"{acronym}_000{ifile}{extension}"
        elif ifile < 100:
            file_name = f"{acronym}_00{ifile}{extension}"
        elif ifile < 1000:
            file_name = f"{acronym}_0{ifile}{extension}"
        else:
            file_name = f"{acronym}_{ifile}{extension}"
    return file_name
