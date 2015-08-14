import sys
import math
from datetime import date
from Bio import Phylo
import networkx as nx

class GraphMaker:
    """
    Usage:
    import graphmaker as gm
    from Bio import Phylo
    t = Phylo.read('/Users/apoon//wip/drttree/data/foo.tree', 'newick')
    g = gm.GraphMaker(t)
    outfile = open('foo.dot', 'w')
    g.init_dotfile(outfile)
    g.draw_edges(outfile, 0.02)
    g.draw_nodes(outfile)
    outfile.close()
    """
    def __init__(self, tree, tip_labels = ['PATID', 'COLDATE'],
                    origin=date(1990,1,1), scaling_factor=75., 
                    colours = {True: 'firebrick', False: 'white'},
                    shapes = {}):
        self.tree = tree  # BioPython Tree object
        self.tip_labels = tip_labels
        
        self.origin = origin
        
        # settings for rendering GraphViz (graphing program)
        self.scaling_factor = scaling_factor
        self.colours = colours
        self.shapes = shapes
        
        # generate dictionary of child->parent node links
        self.parents = {}
        for clade in self.tree.find_clades(order='level'):
            for child in clade:
                self.parents[child] = clade
        
        # gather tips by patid
        self.nodelist = {}
    
    
    def walk_up (self, tips, curnode, pathlen, cutoff):
        """
        Recursive function for traversing up a tree.
        """
        pathlen += curnode.branch_length
        if pathlen < cutoff:
            if curnode.is_terminal():
                tips.append((curnode, pathlen))
            else:
                for c in curnode.clades:
                    tips = self.walk_up(tips, c, pathlen, cutoff)
        return tips
    
        
    def walk_trunk (self, curnode, cutoff):
        """
        Find all tips in the tree that are within a threshold distance
        of a reference tip.
        """
        # first go down to parent and up other branch
        tips = []
        pathlen = curnode.branch_length # 0.0184788
        p = self.parents[curnode]
        
        for c in p.clades:
            if c == curnode: continue
            if c.is_terminal():
                if pathlen + c.branch_length < cutoff:
                    tips.append((c, pathlen + c.branch_length))
            else:
                tips.extend(self.walk_up([], c, pathlen, cutoff))
        
        # next walk down trunk until path length exceeds cutoff or hit root
        while p in self.parents:
            curnode = p
            pathlen += p.branch_length # + 0.0104047
            p = self.parents[curnode]
            if pathlen >= cutoff: break
            for c in p.clades:
                if c == curnode:
                    continue
                if c.is_terminal():
                    if pathlen + c.branch_length < cutoff: # + 0.0503079
                        tips.append((c, pathlen + c.branch_length))
                else:
                    tips.extend(self.walk_up([], c, pathlen, cutoff))
        return tips
    
    def init_dotfile (self, outfile):
        """
        Initialize Graphviz DOT file
        """
        outfile.write('graph foo\n{\n')
        outfile.write('\tnode [shape="circle" style="filled", fontname="Helvetica", penwidth=0.5, regular="true"];\n')
        outfile.write('\trankdir=LR;\n')
        outfile.write('\tedge [color=\"#00000070\"];\n')
        outfile.write('\toutputorder=edgesfirst;\n')
    
    def cluster(self, cutoff):
        """
        Generate clusters based on tip-to-tip (patristic) distance cutoff
        This assumes that each individual (host) is represented by one tip
        in the tree only.  See find_short_edges() for another implementation
        with multiple sequences per patient.
        Returns a list of tuples (tip1, tip2, distance)
        
        cutoff: maximum distance for clustering
        """
        res = []  # a List to be returned
        tips = self.tree.get_terminals()  # returns all tips of the tree object
        
        # iterate over every tip in the tree
        for tip1 in tips:
            # find other tips within cutoff
            #res.update({tip1.name: {}})
            for tip2, dist in self.walk_trunk(tip1, cutoff):
                #res[tip1.name].update({tipname: dist})
                res.append((tip1.name, tipname.name, dist))
    
        return res

    

    
    def find_short_edges(self, cutoff, keep_ties=True, minimize=False):
        """
        Find the shortest edge from the earliest sequence of a patient to a 
        any sequence from any other patient.
        
        cutoff = tip-to-tip distance threshold for defining clusters
        minimize = keep only edge from earliest seq to the closest other seq
        keep_ties = [to be used in conjunction with minimize]
                    report all edges with the same minimum distance
        """
        res = []
        
        for patid in self.patid_to_tips.iterkeys():
            # get earliest sequence for this subject
            tips = self.patid_to_tips[patid]
            tipnames = [t.name for t in tips]
            
            intermed = []
            for tip in tips:
                tokens = tip.name.split('_')
                coldate = tokens[self.tip_labels.index('COLDATE')]
                intermed.append((map(int, coldate.split('-')), tip))
            
            #intermed = [(map(int, tip.name.split('_')[1].split('-')), tip) for tip in tips]
            intermed.sort()
            tip1 = intermed[0][1]
            
            # find the shortest distance in sequences that "cluster" with this one
            min_dist = 99999.
            tip2 = []
            for tipname, dist in self.walk_trunk (tip1, cutoff):
                if tipname in tipnames:
                    # omit sequences from same patient
                    continue
                if minimize and dist < min_dist:
                    min_dist = dist
                    tip2 = [[tipname, dist]]
                else:
                    tip2.append([tipname, dist])
            
            if tip2:
                if keep_ties:
                    for t2, dist in tip2:
                        res.append((tip1.name, t2, dist, True if len(tip2)>1 else False))
                else:
                    # take the first match only
                    if (tip1.name, tip2[0], min_dist, False) not in res:
                        res.append((tip1.name, tip2[0], min_dist, False))
        
        return res
    
    
    def draw_edges(self, outfile, cutoff):
        """
        Output edges to DOT file
        """
        self.nodelist = {} # reset dict
        edgestr = ''
        
        for i in range(len(self.patids)):
            patid1 = self.patids[i]
            tips1 = self.patid_to_tips[patid1]
            
            # all sequences in other subjects that cluster with ANY sequence
            # in the current subject
            cluster = {}
            for t1 in tips1:
                for partner, dist in self.walk_trunk (t1, cutoff):
                    cluster.update ( {partner: (t1.name, dist)} )
            
            # loop over other subjects
            links = []
            for j in range(i+1, len(self.patids)):
                patid2 = self.patids[j]
                tips2 = self.patid_to_tips[patid2]
                
                # which edge is the shortest between these two subjects?
                overlap = set(cluster.keys()).intersection(set([tip.name for tip in tips2]))
                if overlap:
                    min_dist = 99999.
                    for t2 in overlap:
                        t1, dist = cluster[t2]
                        if dist < min_dist:
                            min_dist = dist
                    
                    links.append( (min_dist, patid2) )
            
            if len(links) == 0: continue
            
            # else this patid is linked
            if not self.nodelist.has_key(patid1):
                self.nodelist.update( {patid1:0} )
            
            links.sort()
            for dist, patid2 in links:
                if not self.nodelist.has_key(patid2):
                    self.nodelist.update( {patid2:0} )
                edgestr += '\t"%s"--"%s" [len=%f];\n' % (patid1, patid2, (dist/cutoff+0.1)*4)
        
        outfile.write(edgestr)
    
    
    def draw_nodes (self, outfile):
        """
        Use edge list (string) object returned from make_edgelist() to generate
        a DOT-formatted file at the given file handle.  Node shape and color determined
        by contents of tip labels.
        """
        # only output patient sequences that have one or more links to another patient
        for patid in self.nodelist.iterkeys():
            nodes = self.patid_to_tips[patid]
            within = []
            
            # extract the median sample date and overall Virco result
            any_resistance = False
            risk_factor = None
            all_dates = []
            
            for node in nodes:
                items = node.name.split('_')
                anonid, coldate, vload = items[:3]
                try:
                    year, month, day = map(int, coldate.split('-'))
                    coldate = date(year, month, day)
                    all_dates.append (coldate)
                except:
                    # missing collection date
                    pass
                
                # if ANY samples from the patient has resistance
                if 'R' in items[5:26]: any_resistance = True
                
                # parse risk factors
                if items[-6] == '1' and items[-5] == '1':
                    risk_factor = 'both'
                elif items[-6] == '1':
                    risk_factor = 'idu'
                elif items[-5] == '1':
                    risk_factor = 'msm'
            
            # calculate median collection date
            if len(all_dates) > 0:
                all_dates.sort()
                num_dates = len(all_dates)
                if num_dates%2 == 0:
                    median_date = date.fromordinal((all_dates[num_dates/2-1].toordinal() + all_dates[num_dates/2].toordinal()) / 2)
                else:
                    median_date = all_dates[num_dates/2]
            else:
                median_date = date(2001,1,1) # arbitrary value
            
            outfile.write('\t"%s" [label="", shape="%s", fillcolor="%s", width=%f];\n' % (patid, 
                        self.shapes.get(risk_factor, 'circle'),
                        self.colours.get(any_resistance, 'white'), 
                        math.sqrt(median_date.toordinal()-self.origin.toordinal()) / self.scaling_factor))
        
        outfile.write('}\n')




