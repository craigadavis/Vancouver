import sys
from graphmaker import GraphMaker
from Bio import Phylo
import networkx as nx
from datetime import date, datetime
from csv import DictReader

def main():
    # unpack the command line arguments and process them
    """
    if len(sys.argv) < 4:
        print('Usage:')
        print('python graphmaker.py [<input> Newick tree file] [<input> cutoff] [<output> DOT file]')
    treefile = sys.argv[1]
    cutoff = sys.argv[2]
    dotfile = sys.argv[3]
    try:
        cutoff = float(cutoff)
    except:
        print('ERROR: expecting a float value for [cutoff]')
        sys.exit()
    """
    # point to the data files
    treefile = 'Vancourver_Bref_1302_aligned-out.tree'
    cutoff = 0.02
    dotfile = 'demo.dot'
    
    # load the tree
    try:
        tree = Phylo.read(treefile, 'newick')
    except:
        print('ERROR: tree must be in Newick format')
        raise
    
    # make an instance of GraphMaker object and call it GM
    GM = GraphMaker(tree)  
    
    # generates a dictionary of tips that are similar
    edges = GM.cluster(cutoff)
    
    # generate networkx object
    g = nx.Graph()
    for tip1, tip2, dist in edges:
        g.add_edge(tip1, tip2, dist=dist)
        
    # this code was written for when cluster() returned a dictionary
    #for tip1, neighbours in edges.iteritems():
    #    for tip2, dist in neighbours.iteritems():
    #        # tip2 is a Phylo.Clade object, not a string
    #        g.add_edge(tip1, tip2.name, dist=dist)
    
    # each cluster can be examined using networkx functions; for example,
    clusters = list(nx.connected_components(g))
    intermed = [(len(cluster), cluster) for cluster in clusters]
    intermed.sort(reverse=True)
    
    # this a list of node names (sequence labels)
    cluster0 = [tip if type(tip) is str else tip.name for tip in intermed[0][1]]
    
    # demonstrate how to parse out collection dates 
    # from the StudyID (sequence label (or patid variable))
    for patid in cluster0:
        ddmmyy = patid[-6:]  # last 6 characters of patid
        #day, month, year = [int(ddmmyy[i:(i+2)]) for i in range(0, 6, 2)]
        #coldate = date(year, month, day)  # this does not work :-/
        
        # the %X tokens tell Python how to interpret this string as a date
        # %d means two-digit days
        # %m means two-digit month (02 = February)
        # %y means two-digit year (15 = 2015)
        coldate = datetime.strptime(ddmmyy, '%d%m%y').date()
        
        # number of days since 1980, January 1
        origin = date(1980, 1, 1)
        #print (coldate - origin).days
        #print coldate.year
    
    # read in CSV file from spreadsheet - this is one way of doing it
    csvdata = {}
    handle = open('PHA_phylo_data_May15.csv', 'rU')
    #handle = open('PHA_phylo_data_May15/Sheet1-Table 1.csv', 'rU')
    header = handle.next()  # skip first line
    for line in handle:
        items = line.strip('\n').split(',')  # parse a row from CSV file
        #print items[0]
    handle.close()
    
    # this is a better way to read a CSV file (options for showing single
    # or multiple columns linked via StudyID/sequence name/patid fields )
    csvdata = {}
    reader = DictReader(open('PHA_phylo_data_May15.csv', 'rU'))
    #reader = DictReader(open('PHA_phylo_data_May15/Sheet1   -Table 1.csv', 'rU'))
    exposures = []  # collect all types of exposure values
    for row in reader:
        # csvdata.update({row['StudyID']: {'SumExposure': row['SumExposure']}})
        # csvdata.update({row['StudyID']: {'Expose': row['Expose']}})
        # csvdata.update({row['StudyID']: {'Subtype': row['Subtype']}})
        csvdata.update({row['StudyID']: {
            'Sex': row['Gender'], # MALE/FEMALE/blank
            'Subtype': row['Subtype'], 
            'YROnset': row['YROnset'], # year of onset
            'Expose': row['Expose'],
            'AIDSatDiag': row['LATE?'] # Y or blank
        }})
        #if row['Expose'] not in exposures:
        #    exposures.append(row['Expose'])

    handle.close()
    
    
    # export as GraphViz file using networkx
    #nx.write_dot(g, dotfile)
    
    # now let's write our own DOT file instead of using networkx
    handle = open(dotfile, 'w')  # prepare the file for writing
    handle.write('graph clusters\n{\n')  # opening statement for DOT
    handle.write('outputorder=edgesfirst;\n')  # global option
    handle.write('\tnode [width=0.25,height=0.25,style="filled",fontname="Helvetica",fontsize="9pt"];\n')  # set global options for nodes
    handle.write('\tedge [color="#77777730"];\n')  # global options for edges
    
    # now loop through all graph edges to make connections between nodes
    # note that we are not using the networkx.Graph object!
    """
    drawn_edges = {}  # track what edges are already in the DOT file
    for tip1, tip2, dist in edges:
        # this is how we check whether the reverse pair is already here
        dyad = [tip1, tip2]
        dyad.sort()
        if tuple(dyad) in drawn_edges:
            continue
        drawn_edges.update({tuple(dyad): None})
        
        handle.write('\t"%s"--"%s" [dist=%f];\n' % (
            tip1.replace('~', ''), 
            tip2.replace('~', ''), 
            dist))
    """    
    # now let's use the networkx.Graph object which can split our network
    # up into separate clusters
    nodes_in_clusters = []
    for cluster in clusters:
        if len(cluster) < 5:
            # skip this group!  (includes singletons)
            continue
        # each item is a list of tip names in the same cluster
        sg = g.subgraph(cluster)
        drawn_edges = {}  # see above
        for tip1, tip2, data in sg.edges_iter(data=True):
            dyad = [tip1, tip2]
            dyad.sort()
            if tuple(dyad) in drawn_edges:
                continue
                
            # track nodes in clusters
            if tip1 not in nodes_in_clusters:
                nodes_in_clusters.append(tip1)
            if tip2 not in nodes_in_clusters:
                nodes_in_clusters.append(tip2)
            
            drawn_edges.update({tuple(dyad): None})
            handle.write('\t"%s"--"%s" [len=%f];\n' % (
                tip1.replace('~', ''), 
                tip2.replace('~', ''), 
                data['dist']/0.02+0.5)  # between 0 and 1-ish
            )

    # set up conditional node attributes (how we want them to look)
    # handle.write('\tnode [shape="circle",style="filled",fontname="Helvetica"];\n')
    
    # make a dictionary that specifies node labels for risk factors
    exposures = {
        'OriginHP': 'O', 
        'HeteroWithOriginHP': 'HO', 
        'Unk/Miss': '', 
        'MSM': 'M', 
        'Bi': 'B', 
        'Hetero': 'H', 
        'Hetero+IDU': 'HI',
        'MotherRisk': 'W', 
        'IDU': 'I', 
        'MSM+IDU': 'IM', 
        'Bi+IDU': 'BI', 
        'Other': '', 
        'HeteroWithBi': 'BH', 
        'Hetro': 'H', 
        'HetroWithBi': 'HB', 
        '': ''
    }
    
    # an RGBA Hex code indicates the amount of red, green, blue and alpha
    # with hexadecimal numbers (F=15,FF=255)
    # alpha channel controls transparency level
    # #0000FF30 means 0 red, 0 green, 255 blue, 48 alpha
    gender = {'MALE': '#0000FF30', 'FEMALE': '#FF000030'}
    
    # loop through all subjects to write node statements in DOT
    # write ONLY the nodes that are in clusters
    for patid, data in csvdata.iteritems():  # key/value pairs
        year_onset = data['YROnset']
        if year_onset == 'Y' or year_onset == '':
            # for missing onset year, set to arbitrary value
            year_onset = 2009  # about the median
        else:
            year_onset = int(year_onset)
            if year_onset < 1970:
                year_onset = 2009  # handle as missing value (1 case)
           
        
        patid2 = patid.replace('~', '')  # eliminate tildes
        if patid2 not in nodes_in_clusters:
            continue
            
        handle.write('\t"%s"[label="%s",xlabel="%s",fillcolor="%s",shape="%s",width=%f];\n' % (
            # identifies the node internally (includes date)
            patid2,  
            
            # returns node label indicating exposure
            exposures[data['Expose']],  
            
            # external node label (shortened patid)
            patid2[:-6],  
            
            # colour node based on gender - default (missing) is white
            gender.get(data['Sex'], '#FFFFFF30'),
            
            # indicate if patient was chronic at diagnosis
            'hexagon' if (data['AIDSatDiag'] == 'Y') else 'circle',
            
            # scale node size to year of onset
            max(0., 0.5 * float(year_onset-1977) / (2013-1977))
            )
        )

    handle.write('}\n')  # closing statement
    handle.close()

# ===================================
#  This exists so that when you call this script from the command line it
#  executes main() with your arguments (in sys.argv)
if __name__ == "__main__":
    main()
