import Fasta_one_line as fol
import statistics as st

class NoFasta(Exception):
    pass

def fixedWrite(fastaDict, outHandle, width):
    for header in fastaDict.keys():
        if width !=0:
            seq = "\n".join([fastaDict[header][i:(i+width)] for i in range(0, len(fastaDict[header]), width)])
        else:
            seq = fastaDict[header]
        outHandle.write("{}\n{}\n".format(header,seq))

def combine_fastas(fastaList):
    fastaDict = {}
    for fastaFile in fastaList:
        fastaDict.update(fol.one_line_d(fastaFile))
    return fastaDict

def rename_headers(fastaDict, newComponent, joint=" "):
    newDict = {}
    for header in fastaDict.keys():
        newDict[">" + newComponent + joint + header.lstrip(">")] = fastaDict[header]
    return newDict

def get_seq_len(fastaDict):
    cumulativeLen = 0

    for header in fastaDict.keys():
        cumulativeLen += len(fastaDict[header])
    return cumulativeLen

def get_gc(seq):
    seqList = list(seq.upper())
    g = seqList.count('G')
    c = seqList.count('C')
    return g + c

def stats(fastaDict):
    statDict = {}

    lenList = list(map(len, list(fastaDict.values())))
    lenList.sort(reverse=True)

    statDict['TotalLength'] = sum(lenList)
    statDict['AverageLength'] = round(st.mean(lenList),2)
    statDict['Minimum'] = lenList[-1]
    statDict['Maximum'] = lenList[0]
    statDict['Median'] = st.median(lenList)
    statDict['StDev'] = round(st.pstdev(lenList),2)
    statDict['NumSeqs'] = len(lenList)
    cumulativeGC = 0
    for seq in fastaDict.values():
        cumulativeGC += get_gc(seq)
    statDict['GC'] = round(cumulativeGC/statDict['TotalLength'], 2)
    cumulativeSum = 0
    for length in lenList:
        cumulativeSum += length
        if cumulativeSum >= (statDict['TotalLength']/2.0):
            statDict['N50'] = length
            break
    return statDict

def rem_short(fastaDict, minLen):
    headers = fastaDict.keys()
    retDict = dict(fastaDict)

    for head in headers:
        if len(fastaDict[head]) < minLen:
            del retDict[head]
    return retDict

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Add two fasta files together")

    parser.add_argument("-f", "--fastas", nargs="+", help="Fasta files to process")
    parser.add_argument("-w", "--width", default=0, help="Fix line length at 'width' characters")
    parser.add_argument("-o", "--outfile", default=None, help="Output file")
    parser.add_argument("-n", "--name", type=str, default=None, help="Name to add to headers")
    parser.add_argument("-l", "--length", action='store_true', help="Print total number of nucleotides in fasta")
    parser.add_argument("-s", "--stats", action='store_true', help="Print some basic stats on the fasta file(s)")
    parser.add_argument("-t", "--trim", type=int, default=0, help="Remove seqs under given size")

    args = parser.parse_args()

    if args.fastas is None:
        raise NoFasta("You must input at least one fasta")

    if len(args.fastas) > 1:
        fastaDict = combine_fastas(args.fastas)
    else:
        fastaDict = fol.one_line_d(args.fastas[0])

    if args.trim > 0:
        fastaDict = rem_short(fastaDict, args.trim)

    if args.name is not None:
        fastaDict = rename_headers(fastaDict, args.name)

    if args.length:
        seqLen = get_seq_len(fastaDict)
        print("Cumulative Sequence Length: {}".format(seqLen))

    if args.stats:
        statDict = stats(fastaDict)
        for key in statDict.keys():
            print("{}: {}".format(key, statDict[key]))

    if args.outfile is not None:
        outfile = open(args.outfile, 'w')
        fixedWrite(fastaDict, outfile, args.width)
        outfile.close()

