import dialout

dialout.load_pairs('./dialout/orthogonal_primers_pairs.csv')
print(dialout.get_pair('BBF10K_00001'))

dialout.clear_pairs()