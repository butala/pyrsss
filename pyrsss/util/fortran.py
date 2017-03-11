def descriptor2colspecs(descriptor_str):
    """
    Convert fortran format string, e.g., `2I4,4I3,3I4,2I7,F6.2,I7,
    8F8.2,4F8.1,F7.2,F9.0,F6.2,2F7.2,F6.1,6F8.2,7I6,F7.2,F5.1`, to the
    list of fixed columns tuples expected by :func:`pandas.read_fwf`.
    """
    colspecs = []
    i = 0
    for tok in descriptor_str.split(','):
        tok = tok.strip()
        if tok[0].isdigit():
            tok_chop_left = tok.lstrip('0123456789')
            tok_left = tok[:len(tok) - len(tok_chop_left)]
        else:
            tok_chop_left = tok
            tok_left = '1'
        tok_type = tok_chop_left[0]
        assert tok_chop_left[1].isdigit()
        tok_right = tok_chop_left[1:]
        multiplyer = int(tok_left)
        N = int(tok_right.split('.')[0])
        for n in range(multiplyer):
            colspecs.append((i, i + N))
            i += N
    return colspecs
