#!/usr/bin/env python3
# ConcatFasta.py

import argparse
import os
import re

def main():
    # parse arguments
    parser = argparse.ArgumentParser(prog="ConcatFasta.py",
        description="Concatenates FASTA files into a single file.",
        epilog="""Files can be sepcified one by one with --files
        or by specifying a directory --dir with the files.
        Additionally, a --suffix suffix can be combined with --dir to
        filter out other files.""")
    parser.add_argument(
    '--files', '-f', metavar='FASTA_FILE', nargs="*", type=str,
    help='file to concatenate in FASTA format.')
    parser.add_argument(
    '--dir', '-D', nargs="?", default=".", type=str,
    help='directory where FASTA files to concatenate are (default: %(default)s).')
    parser.add_argument(
    '--suffix', '-s', nargs="?", type=str,
    help='suffix for FASTA files.')
    parser.add_argument(
    '--outfile', '-o', default="concat.fasta", nargs="?", type=str,
    help='name for output file (default: %(default)s).')
    parser.add_argument(
    '--missing_character', '-mc', default="?", nargs="?", type=str,
    help='character used for missing data (default: %(default)s).')
    parser.add_argument(
    '--delim', '-d', nargs="?", type=str,
    help='use delimiter to split FASTA header. The first element will be kept.')
    parser.add_argument(
    '--part', '-q', action="store_true", default=False,
    help='will print a partition table.')
    parser.add_argument(
    '--silent', action="store_true", default=False,
    help='will finish without messages.')
    parser.add_argument(
    '--wrap', '-w', const=100, nargs="?", type=int, default=False, metavar="N",
    help='sequences will be wrapped every N characters. (default: 100)')
    parser.add_argument(
    '--nexus', '-n', action="store_true", default=False,
    help='export in NEXUS format.')
    parser.add_argument(
    '--phylip', '-p', action="store_true", default=False,
    help='export in PHYLIP format.')

    # Welcome message
    print("ConcatFasta v2.0")

    args = parser.parse_args()
    if (args.dir == '.' and args.files == None):
        proceed = input("Do you which to run ConcatFasta on all files in the current directory? [y|n]")
        if proceed == 'n':
            parser.error(message="use either --files/-f or --dir/-d")
    # store data in:
    all_labels = {}
    datalist = {}
    datalen = {}
    if args.nexus and args.outfile == "concat.fasta":
        args.outfile = "concat.nex"
    elif args.phylip and args.outfile == "concat.fasta":
        args.outfile = "concat.phy"
    if args.nexus and args.phylip:
        parser.error(message="pick either NEXUS or PHYLIP formats, but not both")
    # determine the way to read the files
    if args.files is None:
        # the directory way
        files = os.listdir(args.dir)
        if args.suffix is not None:
            files = filter(lambda x: args.suffix in x, files)
        for file in files:
            # read one file at a time, store in dictionary
            if not args.silent:
                print("Reading: ", file, "\r", end='', flush=True)
            datalist[file] = readfasta(args.dir+"/"+file, args.delim)
            # exit if file not read
            if datalist[file] == {}:
                print(file+" is not FASTA")
                parser.error(message=file+" is not FASTA")
            alnlen = list(map(lambda x: len(datalist[file][x]), datalist[file].keys()))
            if not all_same(alnlen):
                parser.error(message="sequences are not the same length")
            for head in datalist[file].keys(): all_labels[head] = ''
            #all_labels.append(list(datalist[file].keys()))
            datalen[file] = alnlen[0]
        if not args.silent:
            print("\nDone reading files.")
        # reduce labels
        #all_labels = sum(all_labels, [])
        #all_labels = reduce(lambda x,y: x+y,all_labels)
        #all_labels = list(set(all_labels))
        all_labels = list(all_labels.keys())
        if not args.silent:
            print("Sorting ",len(all_labels), " labels ... ", end='')
        all_labels.sort()
        if not args.silent:
            print("Done.")
        # do the concatenation
        if not args.silent:
            print("Preparing concatenation ... ", end='')
        catd = catdata(datalist, all_labels, datalen, files, args.missing_character)
        if not args.silent:
            print("Done.")
        # write to file
        if args.nexus:
            exportnexus(catd, args.outfile)
        elif args.phylip:
            exportphylip(catd, args.outfile)
        else:
            writefasta(catd, args.outfile, args.wrap)
        # status to screen
        if not args.silent:
            print("Your concatenated file is "+args.outfile)
        if args.part:
            if args.nexus:
                partblock(args.outfile, datalen, files)
                if not args.silent:
                    print("Partition block added to NEXUS file")
            else:
                printpartition(datalen, files)
    else:
        # the files way.. similar to the previous block
        files = list(map(lambda x: args.dir + "/" + x, args.files))
        for file in files:
            if not args.silent:
                print("Reading: ", file, "\r", end='', flush=True)
            datalist[file] = readfasta(args.dir+"/"+file, args.delim)
            if datalist[file] == {}:
                parser.error(message=file+" is not FASTA")
            alnlen = list(map(lambda x: len(datalist[file][x]), datalist[file].keys()))
            if not all_same(alnlen):
                parser.error(message="sequences are not the same length")
            for head in datalist[file].keys(): all_labels[head] = ''
            #all_labels.append(list(datalist[file].keys()))
            datalen[file] = alnlen[0]
        if not args.silent:
            print("\nDone reading files.")
        #all_labels = sum(all_labels, [])
        #all_labels = reduce(lambda x,y: x+y,all_labels)
        #all_labels = list(set(all_labels))
        all_labels = list(all_labels.keys())
        if not args.silent:
            print("Sorting",len(all_labels), "labels ... ", end='')
        all_labels.sort()
        if not args.silent:
            print("Done.")
        # do concatenation
        if not args.silent:
            print("Preparing concatenation ... ", end='')
        catd = catdata(datalist, all_labels, datalen, files, args.missing_character)
        if not args.silent:
            print("Done.")
        # write to file
        if args.nexus:
            exportnexus(catd, args.outfile)
        elif args.phylip:
            exportphylip(catd, args.outfile)
        else:
            writefasta(catd, args.outfile, args.wrap)
        # print status to screen
        if not args.silent:
            print("Your concatenated file is "+args.outfile)
        # print partitions
        if args.part:
            if args.nexus:
                partblock(args.outfile, datalen, files)
                if not args.silent:
                    print("Partition block added to NEXUS file")
            else:
                printpartition(datalen, files)

# functions

def readfasta(file, delim):
    data = {}
    if delim:
        with open(file, "r") as f:
            for line in f:
                line = line.rstrip()
                if len(line) != 0:
                    if line[0]:
                        head = line[1:]
                        data[head.split(delim)[0]] = ''
                    else:
                        data[head.split(delim)[0]] += re.sub(" ","",line)
            return data
    else:
        with open(file, "r") as f:
            for line in f:
                line = line.rstrip()
                if len(line) != 0:
                    if line[0] == ">":
                        head = line[1:]
                        data[head] = ''
                    else:
                        data[head] += re.sub(" ","",line)
            return data

def writefasta(catd, outf, wrap):
    with open(outf,"w") as o:
        for i in catd.keys():
            o.write(">"+i+"\n")
            if not wrap:
                o.write(catd[i]+"\n")
            else:
                o.write(wrapseq(catd[i], wrap)+"\n")

def exportnexus(data, outf):
    labels = list(data.keys())
    maxlen = max([ len(i) for i in labels ])
    spaced = [ i + ' ' * (maxlen - len(i) + 1) for i in labels ]
    ntax = len(labels)
    nchar = len(data[labels[0]])
    with open(outf,"w") as o:
        o.write("#NEXUS\n\n")
        o.write("Begin DATA;\n")
        o.write(f"\tDimensions ntax={ntax} nchar={nchar};\n")
        o.write("\tFormat Datatype=DNA gap=- missing=?;\n\tMatrix\n")
        for i in zip(labels,spaced):
            o.write("\t"+i[1]+data[i[0]]+"\n")
        o.write(";\nEnd;\n")

def exportphylip(data, outf):
    labels = list(data.keys())
    maxlen = max([ len(i) for i in labels ])
    spaced = [ i + ' ' * (maxlen - len(i) + 1) for i in labels ]
    ntax = len(labels)
    nchar = len(data[labels[0]])
    with open(outf,"w") as o:
        o.write(f"\t{ntax} {nchar}\n")
        for i in zip(labels,spaced):
            o.write("\t"+i[1]+data[i[0]]+"\n")

def catdata(data, labels, seqlen, files, mc):
    cdat = {}
    for lab in labels:
        cdat[lab] = ''
        for dat in files:
            if lab in data[dat].keys():
                cdat[lab] += data[dat][lab]
            else:
                cdat[lab] += mc * seqlen[dat]
    return cdat

def printpartition(seqlen, files):
    if len(seqlen) != len(files):
        print(seqlen)
        print(files)
        parser.error(message="not the same number of items")
    seqlen = list(map(lambda x: seqlen[x], files))
    gnames = [ re.sub(".[fF][aA]{0,1}[sS]{0,1}[tT]{0,1}[aA]{0,1}$","",i.split('/')[-1]) for i in files ]
    prev = 1
    with open("part.txt", "w") as o:
        for i in range(len(seqlen)):
            sumlen = sum(seqlen[0:i+1])
            L = gnames[i] + " = " + str(prev) + "-" + str(sumlen) + ";"
            print(L)
            o.write(L+"\n")
            prev = sumlen+1

def partblock(outf, seqlen, files):
    if len(seqlen) != len(files):
        print(seqlen)
        print(files)
        parser.error(message="not the same number of items")
    with open(outf, "a") as o:
        o.write("\nBegin Sets;\n")
        gnames = [ i.split('/')[-1].split(".")[0] for i in files ]
        seqlen = list(map(lambda x: seqlen[x], files))
        prev = 1
        for i in range(len(seqlen)):
            sumlen = sum(seqlen[0:i+1])
            o.write("\tcharset " + gnames[i] + " = " + str(prev) + "-" + str(sumlen) + ";\n")
            prev = sumlen+1
        o.write("\n\tpartition all = "+str(len(gnames))+":"+",".join(gnames)+";\n")
        o.write("End;\n")

def wrapseq(seq, w):
    chunks = []
    interval = list(map(lambda x: x*w, range((int(len(seq)/w))+2)))
    for i in interval:
        if i != interval[-1]:
            chunks.append(seq[i:interval[interval.index(i)+1]-1])
    return "\n".join(chunks)

def all_same(items):
    return all(x == items[0] for x in items)

if __name__ == '__main__':
    main()
