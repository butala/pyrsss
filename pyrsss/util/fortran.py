def descriptor2colspecs(descriptor_str):
    """
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
        # print('----')
        # print(tok)
        # print(tok_chop_left)
        # print(tok_left, '???')
        tok_type = tok_chop_left[0]
        assert tok_chop_left[1].isdigit()
        tok_right = tok_chop_left[1:]
        # print(tok_right)
        multiplyer = int(tok_left)
        N = int(tok_right.split('.')[0])
        #tok_len = multiplyer * N
        for n in range(multiplyer):
            colspecs.append((i, i + N))
            i += N
    return colspecs
