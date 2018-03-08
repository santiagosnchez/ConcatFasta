# ConcatFasta.py

import argparse
from os import listdir
from re import match

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
    '--dir', '-d', nargs="?", default=".", type=str,
    help='directory where FASTA files to concatenate are (default: %(default)s).')
    parser.add_argument(
    '--suffix', '-s', nargs="?", type=str,
    help='suffix for FASTA files.')
    parser.add_argument(
    '--outfile', '-o', default="concat.fasta", nargs="?", type=str,
    help='name for output file (default: %(default)s).')
    parser.add_argument(
    '--part', '-q', action="store_true", default=False,
    help='will print a partition table.')
    parser.add_argument(
    '--wrap', '-w', const=100, nargs="?", type=int, default=False, metavar="N",
    help='sequences will be wrapped every N characters. (default: 100)')
    parser.add_argument(
    '--nexus', '-n', action="store_true", default=False,
    help='export in NEXUS format.')
    parser.add_argument(
    '--phylip', '-p', action="store_true", default=False,
    help='export in PHYLIP format.')
    args = parser.parse_args()
    if (args.dir == '.' and args.files == None):
        proceed = raw_input("Do you which to run ConcatFasta on all files in the current directory? [y|n]")
        if proceed == 'n':
            parser.error(message="use either --files/-f or --dir/-d")
    # store data in:
    all_labels = []
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
        files = listdir(args.dir)
        if args.suffix is not None:
            files = filter(lambda x: args.suffix in x, files)
        for file in files:
             # read one file at a time, store in dictionary
            datalist[file] = readfasta(args.dir+"/"+file)
            # exit if file not read
            if datalist[file] == {}:
                print file+" is not FASTA\n"
                parser.error(message=file+" is not FASTA")
            alnlen = map(lambda x: len(datalist[file][x]), datalist[file].keys())
            if not all_same(alnlen):
                parser.error(message="sequences are not the same length")
            all_labels.append(datalist[file].keys())
            datalen[file] = alnlen[0]
        # reduce labels
        all_labels = reduce(lambda x,y: x+y,all_labels)
        all_labels = list(set(all_labels))
        all_labels.sort()
        # do the concatenation
        catd = catdata(datalist, all_labels, datalen)
        # write to file
        if args.nexus:
            exportnexus(catd, args.outfile)
        elif args.phylip:
            exportphylip(catd, args.outfile)
        else:
            writefasta(catd, args.outfile, args.wrap)
        # print status to screen
        print "Your concatenated file is "+args.outfile
        if args.part:
            if args.nexus:
                partblock(args.outfile, datalen, files)
                print "Partition block added to NEXUS file"
            else:
                printpartition(datalen, files)
    else:
        # the files way.. similar to the previous block
        files = map(lambda x: args.dir + "/" + x, args.files)
        for file in files:
            datalist[file] = readfasta(args.dir+"/"+file)
            if datalist[file] == {}:
                parser.error(message=file+" is not FASTA")
            alnlen = map(lambda x: len(datalist[file][x]), datalist[file].keys())
            if not all_same(alnlen):
                parser.error(message="sequences are not the same length")
            all_labels.append(datalist[file].keys())
            datalen[file] = alnlen[0]
        all_labels = reduce(lambda x,y: x+y,all_labels)
        all_labels = list(set(all_labels))
        all_labels.sort()
        catd = catdata(datalist, all_labels, datalen)
        # write to file
        if args.nexus:
            exportnexus(catd, args.outfile)
        elif args.phylip:
            exportphylip(catd, args.outfile)
        else:
            writefasta(catd, args.outfile, args.wrap)
        # print status to screen
        print "Your concatenated file is "+args.outfile
        # print partitions
        if args.part:
            if args.nexus:
                partblock(args.outfile, datalen, files)
                print "Partition block added to NEXUS file"
            else:
                printpartition(datalen, files)

# functions

def readfasta(file):
    data = {}
    with open(file, "r") as f:
        for line in f:
            line = line.rstrip()
            if match("^>",line):
                head = line[1:]
                data[head] = ''
            else:
                data[head] += line
        return data

def writefasta(catd, outf, wrap):
    o = open(outf,"w")
    for i in catd.keys():
        o.write(">"+i+"\n")
        if not wrap:
            o.write(catd[i]+"\n")
        else:
            o.write(wrapseq(catd[i], wrap)+"\n")
    o.close()

def exportnexus(data, outf):
    labels = data.keys()
    maxlen = max([ len(i) for i in labels ])
    spaced = [ i + ' ' * (maxlen - len(i) + 1) for i in labels ]
    ntax = len(labels)
    nchar = len(data[labels[0]])
    o = open(outf,"w")
    o.write("#NEXUS\n\n")
    o.write("Begin DATA;\n")
    o.write("\tDimensions ntax=%s nchar=%s;\n" % (ntax,nchar))
    o.write("\tFormat Datatype=DNA gap=- missing=?;\n\tMatrix\n")
    for i in zip(labels,spaced):
        o.write("\t"+i[1]+data[i[0]]+"\n")
    o.write(";\nEnd;\n")
    o.close()

def exportphylip(data, outf):
    labels = data.keys()
    maxlen = max([ len(i) for i in labels ])
    spaced = [ i + ' ' * (maxlen - len(i) + 1) for i in labels ]
    ntax = len(labels)
    nchar = len(data[labels[0]])
    o = open(outf,"w")
    o.write("\t%s %s\n" % (ntax,nchar))
    for i in zip(labels,spaced):
        o.write("\t"+i[1]+data[i[0]]+"\n")
    o.close()

def catdata(data, labels, seqlen):
    cdat = {}
    for lab in labels:
        cdat[lab] = ''
        for dat in data.keys():
            if lab in data[dat].keys():
                cdat[lab] += data[dat][lab]
            else:
                cdat[lab] += '?' * seqlen[dat]
    return cdat

def printpartition(seqlen, files):
    if len(seqlen) != len(files):
        print seqlen
        print files
        parser.error(message="not the same number of items")
    seqlen = map(lambda x: seqlen[x], files)
    gnames = [ i.split('.')[-2][1:] for i in files ]
    prev = 1
    o = open("part.txt", "w")
    for i in range(len(seqlen)):
        sumlen = sum(seqlen[0:i+1])
        L = gnames[i] + " = " + str(prev) + "-" + str(sumlen) + ";"
        print(L)
        o.write(L+"\n")
        prev = sumlen+1
    o.close()

def partblock(outf, seqlen, files):
    if len(seqlen) != len(files):
        print seqlen
        print files
        parser.error(message="not the same number of items")
    o = open(outf, "a")
    o.write("\nBegin Sets;\n")
    gnames = [ i.split('.')[-2][1:] for i in files ]
    seqlen = map(lambda x: seqlen[x], files)
    prev = 1
    for i in range(len(seqlen)):
        sumlen = sum(seqlen[0:i+1])
        o.write("\tcharset " + gnames[i] + " = " + str(prev) + "-" + str(sumlen) + ";\n")
        prev = sumlen+1
    o.write("\n\tpartition all = "+str(len(gnames))+":"+",".join(gnames)+";\n")
    o.write("End;\n")
    o.close()

def wrapseq(seq, w):
    chunks = []
    interval = map(lambda x: x*w, range((len(seq)/w)+2))
    for i in interval:
        if i != interval[-1]:
            chunks.append(seq[i:interval[interval.index(i)+1]-1])
    return("\n".join(chunks))

def all_same(items):
    return all(x == items[0] for x in items)


if __name__ == '__main__':
    main()


