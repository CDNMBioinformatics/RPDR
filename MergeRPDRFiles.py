import math
import os
from os import walk
import os.path
from os import path
#from os import listdir
#from os.path import isfile, join
import argparse
import re

def mergefiles(input_path, output_folder_name, nSubfolders, select_files = []):
    output_path = "%s%s" % (input_path, output_folder_name)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    dirs = []
    for count in range(1, nSubfolders+1):
        dirs.append("%s%i of %i/" % (input_path, count, nSubfolders))
    #Get list of files to read
    _, _, all_filenames = next(walk(dirs[0]))
    if len(select_files) > 0:
        filenames_temp = []
        for sf in select_files:
            filenames_temp.append([fn for fn in all_filenames if re.search("_%s" % sf, fn, re.IGNORECASE)])
        # flatten the list depending on if it is based on single name or first letter which could have multiple items
        filenames = []
        for ele in filenames_temp:
            # Case when using single letter
            if len(ele) > 1:
                filenames.extend([item for item in ele])
            # Case when using 3-letter name or if single letter returns 1 file
            else:
                filenames.extend(ele)
    else:
        filenames = all_filenames
    print(filenames)
    for fn in filenames:
        # check if file exists in all directories
        exists = []
        for folder in dirs:
            exists.append(path.exists("%s%s" % (folder, fn)))
        if all(exists):
            fullnames = [folder + fn for folder in dirs]
            print("Merging %s..." % fn)
            fp_merge = open("%s%s" % (output_path, fn), "w")
            # Read in files from directory 1 of n
            with open("%s1 of %s/%s" % (input_path, nSubfolders, fn)) as fp_1_n:
                for line in fp_1_n:
                    fp_merge.write(line)
            print("1 of %s complete" % nSubfolders)
            # read in remaining file directories
            for count in range(2, nSubfolders + 1):
                with open("%s%i of %s/%s" % (input_path, count, nSubfolders, fn)) as fp_count_n:
                    #skip header line since it should be the same
                    next(fp_count_n)
                    for line in fp_count_n:
                        fp_merge.write(line)
                print("%s of %s complete" % (count, nSubfolders))
            fp_merge.close()
            print("Merge complete")
        else:
            print("File is not listed in all directories. Skipping")

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_path", type = str, help= "Provide directory path to folders with subset files")
    parser.add_argument("-o", "--output_folder_name", type = str, default = "Full_Set/", help="Name of merged file (default = Full_Set)")
    parser.add_argument("-n", "--number_subfolders", type = int, help="The number of subset files")
    parser.add_argument("-s", "--select_files", type = str, nargs = "+", help = "list of specific files (example Bib, Dem, or Med) or first letter of item")
    args = parser.parse_args()
    if not(args.input_path.endswith("/")):
        args.input_path = "%s/" % (args.input_path)
    if not (args.output_folder_name.endswith("/")):
        args.output_folder_name = "%s/" % (args.output_folder_name)
    if not args.number_subfolders:
        print("No subfolder count specified. Please run again")
    else:
        if not(args.select_files):
            mergefiles(args.input_path, args.output_folder_name, args.number_subfolders)
        else:
            mergefiles(args.input_path, args.output_folder_name, args.number_subfolders, args.select_files)

