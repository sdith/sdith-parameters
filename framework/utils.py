def print_title(title, size=56):
    nb_spaces = size - len(title) - 4
    nb_spaces_before = nb_spaces//2
    nb_spaces_after  = nb_spaces - nb_spaces_before
    print('*'*size)
    print('* '+' '*nb_spaces_before+title+' '*nb_spaces_after+' *')
    print('*'*size)
