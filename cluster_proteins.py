import Fasta_one_line as fol

infile = open("medium_clusters.txt", 'r')
fasta = fol.blast_dict("/work/mpesesky/Plasmids/NCBI_Plasmids/taxonAnalysis/transposases.faa")

used_families = []

for line in infile:
    nodes = line.rstrip().split("\t")
    clusterName = nodes[0]
    if clusterName in used_families:
        print(clusterName)
        exit()
    outfileName = clusterName + ".faa"
    outfile = open(outfileName, 'w')

    for node in nodes[1:]:
        if "|" not in node:
            continue
        if node.startswith("ref"):
            nodeName = node
        else:
            nodeName = "ref|{}|".format(node.split("|")[1])

        try:
            seq = fasta[nodeName]
        except KeyError:
            print(nodeName)
            exit()
        outfile.write(">{}\n{}\n".format(nodeName, seq))
    outfile.close()
    used_families.append(clusterName)
infile.close()
