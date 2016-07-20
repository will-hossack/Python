import os
import tio as t

def main():
    d = t.getString("Dir")
    d = t.getExpandedFilename(d)
    for filename in os.listdir(d):
        t.tprint(filename)

main()
