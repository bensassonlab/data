import sys, getopt
import random

def main(argv):
    infile = ''
    gendist = ''
    arg_list = sys.argv[1:]
    short_options = 'i:g:h:'
    long_options = ['infile=', 'gendist=', 'help']
    usage = 'removeSimilarStrains.py\n Usage:\n -h, --help \n -i, --infile, genetic distance matrix [required] \n -g, --gendist, genetic distance threshold [required] \n'

    try:
        options, args = getopt.getopt(arg_list,short_options,long_options)
    except getopt.GetoptError:
        print usage
        sys.exit(2)

    if len(arg_list)<1 or not options:
        print usage
        sys.exit(2)

    print '\n', 'Removing similar strains based on genetic distance threshold. Output will be written to outfile.txt', '\n'

    for arg, val in options:
        if arg in ('-i', '--infile'):
            print 'Using genetic distance matrix in ', val
            infile = val
        elif arg in ('-g', '--gendist'):
            print 'Using genetic distance threshold:', val
            gendist = val
        elif arg in ('-h', '--help'):
            print 'This program removes similar strains.\n'
            print usage
            sys.exit()
        else:
            print usage
            sys.exit()

    with open(infile, 'r') as f:
            distmat = f.readlines()

    # create list of strains with same indices as distance matrix
    strainlist = [distmat[i].split()[0] for i in range(len(distmat))]

    # set genetic distance threshold based on use input
    threshold = float(gendist)

    # initiate empty lists and break flag
    indexlist = [] # temp list of matrix indices
    similar = [] # final list of matrix indices
    newlist = [] # converted list from indices to strain IDs
    tokeep = [] # list of strains to to keep
    isbreak = 0 # loop break flag

    # loop through distance matrix and get strains with genetic distance below threshold
    for i in range(len(distmat)):
        isbreak = 0
        for item in similar:
            if(i in item):
                isbreak = 1
        if(isbreak == 1):
             continue
        c_mat = [float(item) for item in distmat[i].split()[1:]]
        for j in range(len(c_mat)):
            if(c_mat[j] <= threshold and j > i):
                indexlist.append(j)
        if indexlist:
            indexlist.insert(0,i)
            #avggd = Average([c_mat[i] for i in indexlist])
            #indexlist = [strainlist[id] for id in indexlist]
            #indexlist.append(avggd)
            similar.append(indexlist)
        indexlist = []

    # convert matrix indices to strain IDs
    for i in similar:
        newlist.append([strainlist[j] for j in i])

    # randomly select one strain per group of similar strains and add to a list called tokeep
    random.seed(12345)
    for i in newlist:
        tokeep.append(random.choice(i))

    # convert list of lists with similar strains to flat list to facilitate list operations
    flat_list = [item for sublist in newlist for item in sublist]

    # get strains to remove from list of similar strains based on strains to keep
    to_remove = list(set(flat_list) - set(tokeep))

    # get selected strains from full strain list
    selected_strains = list(set(strainlist) - set(to_remove))
    selected_strains.sort()

    outprefix = "outfile-" + str('{:f}'.format(threshold)) + ".txt"
    with open(outprefix, 'w') as f:
        for item in selected_strains:
            print >> f, item

    print('\n')
    print 'There are ', len(newlist), ' groups of similar strains (pairwise genetic distance <', '{:f}'.format(threshold), '):\n'
    print('\n'.join(' '.join(map(str,sl)) for sl in newlist))

    print '\n', len(to_remove), 'strains were removed, and the output file contains', len(selected_strains), 'strains.'


if __name__ == "__main__":
    main(sys.argv[1:])
