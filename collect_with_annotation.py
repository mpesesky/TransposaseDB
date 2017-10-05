import Fasta_one_line as fol
import argparse

def gather_seqs(seqDict, targetAnnotation):
    retDict = {}

    for header in seqDict.keys():
        if targetAnnotation.lower() in header.lower():
            retDict[seqDict[header]] = header
    invDict = {v: k for k, v in retDict.items()}

    return invDict

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset a fasta file based on header keywords")

    parser.add_argument("Fasta", help="Input fasta file")
    parser.add_argument("Keyword", type=str, help="Word or phrase used to identify sequences to be collected.")
    parser.add_argument("-o", "--output", default="subset.fna", type=argparse.FileType('w'), help="Output file name")

    args = parser.parse_args()

    inDict = fol.one_line_d(args.Fasta)
    subDict = gather_seqs(inDict, args.Keyword)

    for header in subDict.keys():
        args.output.write("{}\n{}\n".format(header, subDict[header]))
    args.output.close()
