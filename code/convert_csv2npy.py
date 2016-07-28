#!/usr/bin/env python
''' read csv file and save it to npy '''

import sys
import getopt
import numpy as np

def main(argv):
    ''' read csv file and save it to npy '''

    input_csv_file = ''

    try:
        opts, _ = getopt.getopt(argv, "hi:", "input_csv_file=")

    except getopt.GetoptError:
        print 'convert_csv2npy.py -i <input_csv_file>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'convert_csv2npy.py -i <input_csv_file>'
            sys.exit()
        elif opt in ("-i", "--input_csv_file"):
            print arg
            input_csv_file = arg

    tmp = np.loadtxt(input_csv_file, skiprows=1, delimiter=',')
    output_npy_file = input_csv_file.replace('.csv', '.npy')
    np.save(output_npy_file, tmp)

    print '%s is written' %output_npy_file

if __name__ == "__main__":
    main(sys.argv[1:])

