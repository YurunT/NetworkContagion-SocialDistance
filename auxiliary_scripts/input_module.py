import argparse

def input():
    print("This is input.py")
    

"""
Everything changed/added by Yurun
"""
def parse_args(args):
    parser = argparse.ArgumentParser(description = 'Parameters')
    parser.add_argument('-n', type = int, default = 200000, help='200,000 (default); the number of nodes')
    parser.add_argument('-e', type = int, default = 50, help='50 (default); the number of experiments')
    parser.add_argument('-m', type = float, default = 0.6, help='0.6 (default); the prob of wearing a mask')
    parser.add_argument('-tm1', type = float, default = 0.5, help='0.5 (default); T_mask1')
    parser.add_argument('-tm2', type = float, default = 0.5, help='0.5 (default); T_mask2')
    parser.add_argument('-T', type = float, default = 0.6, help='0.6 (default); transmissibility of the orinigal virus')
    parser.add_argument('-th', type = float, default = 0.001, help='0.001 (default); the treshold to consider a component giant')
    parser.add_argument('-nc', type = int, default = 40, help='number of Cores')
    parser.add_argument('-mind', type = int, default = 0, help='[min_degree, max_degree]')
    parser.add_argument('-maxd', type = int, default = 10, help='[min_degree, max_degree]')
    parser.add_argument('-ns', type = int, default = 50, help='Num of sample within [min_degree, max_degree]')
    parser.add_argument('-cp', type = int, default = 500, help='Num of exps to save a checkpoint')
    parser.add_argument('-change', type = int, default = 0, help='Change m(0), T(1), tm(2)')
    return parser.parse_args(args)