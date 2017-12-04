# ConcatFasta.py

import argparse
from sys import exit
from os import listdir

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
    '--part', '-p', const=True, nargs="?", type=bool, default=False, metavar="",
    help='will print a partition table.')
    parser.add_argument(
    '--wrap', '-w', const=True, nargs="?", type=bool, default=False, metavar="",
    help='sequences will be wrapped every 100 characters.')
    args = parser.parse_args()
    # store data in:
    all_labels = []
    datalist = {}
    datalen = {}
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
                exit( parser.print_help())
            alnlen = map(lambda x: len(datalist[file][x]), datalist[file].keys())
            if not all_same(alnlen):
                exit("Sequences are not the same length")
            all_labels.append(datalist[file].keys())
            datalen[file] = alnlen[0]
        # reduce labels
        all_labels = reduce(lambda x,y: x+y,all_labels)
        all_labels = list(set(all_labels))
        all_labels.sort()
        # do the concatenation
        catd = catdata(datalist, all_labels, datalen)
        # write to file
        writefasta(catd, args.outfile, args.wrap)
        # print status to screen
        print "Your concatenated file is "+args.outfile
        if args.part == True:
            printpartition(datalen, files)
    else:
        # the files way.. similar to the previous block
        files = map(lambda x: args.dir + "/" + x, args.files)
        for file in files:
            datalist[file] = readfasta(args.dir+"/"+file)
            if datalist[file] == {}:
                print file+" is not FASTA\n"
                exit( parser.print_help())
            alnlen = map(lambda x: len(datalist[file][x]), datalist[file].keys())
            if not all_same(alnlen):
                exit("Sequences are not the same length")
            all_labels.append(datalist[file].keys())
            datalen[file] = alnlen[0]
        all_labels = reduce(lambda x,y: x+y,all_labels)
        all_labels = list(set(all_labels))
        all_labels.sort()
        catd = catdata(datalist, all_labels, datalen)
        writefasta(catd, args.outfile, args.wrap)
        print "Your concatenated file is "+args.outfile
        if args.part == True:
            printpartition(datalen, files)


def readfasta(file):
    from re import match
    f = open(file, "r")
    data = {}
    seq = ''
    head = ''
    for line in f:
        line = line.rstrip()
        if match("^>",line):
            if len(head) != 0:
                data[head] = seq
                head = line[1:]
            else:
                head = line[1:]
            if len(seq) != 0:
                seq = ''
        else:
            seq += line
    f.close()
    return data

def writefasta(catd, outf, wrap):
    o = open(outf,"w")
    for i in catd.keys():
        o.write(">"+i+"\n")
        if wrap == True:
            o.write(wrapseq(catd[i])+"\n")
        else:
            o.write(catd[i]+"\n")
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
        print "Not the same number of items"
        exit( parser.print_help() )
    seqlen = map(lambda x: seqlen[x], files)
    prev = 1
    for i in range(len(seqlen)):
        sumlen = sum(seqlen[0:i+1])
        print files[i] + " = " + str(prev) + "-" + str(sumlen) + ";"
        prev = sumlen+1

def wrapseq(seq):
    chunks = []
    interval = map(lambda x: x*100, range((len(seq)/100)+2))
    for i in interval:
        if i != interval[-1]:
            chunks.append(seq[i:interval[interval.index(i)+1]-1])
    return("\n".join(chunks))

def all_same(items):
    return all(x == items[0] for x in items)


if __name__ == '__main__':
    main()


