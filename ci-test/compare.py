import sys
import filecmp

file1 = sys.argv[1]
file2 = sys.argv[2]

if not( filecmp.cmp(file1, file2, shallow=True) ):
    sys.exit("Files " + file1 + " and " + file2 + " are not the same")
else:
    print("Files " + file1 + " and " + file2 + " are the same")