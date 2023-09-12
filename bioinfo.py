#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  4 15:46:08 2023

@author: mcatolos, stg, hohuqu, perryk
"""

from mygraph_final import fs_graph
import os
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

#%%#
class SeqRecord:
    '''
    An object that encapsulates both the sequence and features of a given
    genome.  Identified by an accession number for the organism (e.g. U00096.3).
    A SeqRecord object contains one Sequence object and a list of Feature
    objects for that sequence.
    '''
    
    '''Attribute for id count (for features w/o explicit id)'''
    id_count = 0
    
    def __init__(self,genome_directory,sequence='',features=[]): #remove: genome_directory,feature_list=[],sequence,features,
        '''Creates FileNode and FileSystem for genome directory'''
        self.rootdir = fs_graph.FileNode(genome_directory)
        
        '''Gets the directory for the .gff & .fa file and other attributes'''
        
        # Get list of all files in directory
        file_list = os.listdir(genome_directory)
        
        # Filter list for .gff files
        gff_files = [file for file in file_list if file.endswith('.gff')]
        
        
        self.accession = self.rootdir.basename()
        self.filedir_gff = os.path.join(str(self.rootdir) + "/" + str(gff_files[0]))
        #self.filedir_gff = os.path.join(str(self.rootdir),str(self.accession) + '.gff')
        #self.filedir_fa = os.path.join(str(self.rootdir),str(self.accession) + '.fa')
        
        '''Sequence and Feature attributes'''
        #self.sequence = 0
        #self.sequence = Sequence(self.filedir_fa)
        self.features = features
        
        '''Attribute for filter runs'''
        self.filter_runs = 0
    
    def __iter__(self): #iterable support
        self.count = -1
        return self
    
    def __next__(self): #iterable support
        while self.count <= len(self.features):
            self.count += 1
            return self.features[self.count].read
        raise(StopIteration)
    
    def parse_gff(self):
        '''
        Takes a path to the corresponding gff file and returns a list of
        feature objects (can include derived feature object)
        '''
        self.opt_attributes = []
        self.attributes = ['seqname','source','type','start','end','score','strand','frame']
        
        gff_file = open(self.filedir_gff)
        print("Success! Now parsing file...")
        with gff_file as f:
            for line in f:
                feature = Feature(line)
                subtypes = Subtypes(line)
                self.features.append(feature)
                self.features.append(subtypes)
                
                SeqRecord.id_count += 1
     
                '''finds optional attributes'''
                s = line.rstrip('\n') #Removes new line 
                s = s.split('\t') #Removes tabs
                s_opt = s[8].replace(';','~') #s[8] refers to optional attributes
                s_opt = s_opt.replace('=','~')
                s_opt = s_opt.split('~')
                self.opt_attributes.append(s_opt)
        print("Parsed " + str(int(len(self.features)/2)) + " features")
        print("\n")
        return self.features
    
    def filter(self,search):
        '''
        - Filter out the features for only those that match the arguments provided
        - Create and return a new SeqRecord object with the same sequence but
        only the features that pass the filter
        - How you pass arguments to filter() is up to you, but it should mimic
        the results of your gff_filter.py script.  In fact, you should consider
        pulling code from gff_filter.py.  And it should make it easy to
        implement seqrec_filter.py as described below.
        '''
        if search == '' or search == 'q':
            print("Thank you using our services! Please come again soon!")
            return ''
        elif search == 'back' or search == 'fasta' or search == 'full':
            attribute = search
        
        strings = ['back','fasta','full']
        search = search.split()
        if len(search) == 3:
            attribute = search[0]
        elif len(search) == 2:
            attribute = search[0]
            value = search[1]
            
        results_list = []
        j = 0
        for i in range(int(len(self.features)/2)):
            #combines the main and optional attributes of features
            l = self.features[i+j].read + "\t" + self.features[i+j+1].read
            l = l.split("\t")

            if len(search) == 3: #start/end pos1 pos2
                if attribute == 'start':
                    start = int(l[3])
                    if (start >= int(search[1]) and start <= int(search[2])):
                        results_list.append(self.features[i+j])
                        results_list.append(self.features[i+j+1])
                elif attribute == 'end':
                    end = int(l[4])
                    if (end >= int(search[1]) and end <= int(search[2])):
                        results_list.append(self.features[i])
            elif len(search) == 2: #attribute string
                #Main Attributes
                if self.attributes.count(attribute) != 0:
                    attribute_index = self.attributes.index(attribute)
                    if l[attribute_index] == value:
                        results_list.append(self.features[i+j])
                        results_list.append(self.features[i+j+1])
                else: #Extra Attributes
                    n = 0
                    for k in range(int(len(self.opt_attributes[i])/2)):
                        if self.opt_attributes[i][n+k] == attribute and l[-1].find(value) != -1:
                            #results_list.append(self.features[i])
                            results_list.append(self.features[i+j])
                            results_list.append(self.features[i+j+1])
                        n += 1
            elif attribute == 'back': #back
                if self.filter_runs == 0:
                    results_list = self.features.copy()
                elif self.filter_runs%2 == 0:
                    results_list = self.previous_results_list_even.copy()
                elif self.filter_runs%2 == 1:
                    results_list = self.previous_results_list_odd.copy()    
            elif attribute == 'fasta' or attribute == 'full':
                if self.filter_runs%2 == 1:
                    results_list = self.previous_results_list_even.copy()
                elif self.filter_runs%2 == 0:
                    results_list = self.previous_results_list_odd.copy()
                else:
                    results_list = self.features.copy()
            j += 1
        
        '''Sets copies'''
        if self.filter_runs%2 == 0:
            self.previous_results_list_even = results_list.copy()
        else:
            self.previous_results_list_odd = results_list.copy()
        self.filter_runs += 1
        
        '''Prints. has different functions for fasta, full, back'''
        h = 0
        if attribute == 'fasta':
            print("Found " + str(int(len(self.filtered.features)/2)) + " features!")
            for l in range(int(len(self.filtered.features)/2)):
                #Feature.fasta(results_list[l+h].id, self.sequence.seq,self.filtered.features[l+h])
                Feature.fasta(results_list[l+h].id,self.filtered.features[l+h])
                h += 1
        elif attribute == 'full':
            print("Found " + str(int(len(results_list)/2)) + " features!")
            for l in range(int(len(results_list)/2)):
                Feature.full(results_list[l+h].id,results_list[l+h],results_list[l+h+1])
                h += 1
        else: 
            print("Found " + str(int(len(results_list)/2)) + " features!")
            for l in range(int(len(results_list)/2)):
                Feature.short(results_list[l+h].id,results_list[l+h])
                h += 1

        #self.filtered = SeqRecord(str(self.rootdir),self.sequence,features=results_list)
        self.filtered = SeqRecord(str(self.rootdir),features=results_list)
        return self.filtered

"""
class Sequence:
    '''
    An object that encapsulates a DNA sequence.
    Identified by an accession number (e.g. U00096) for that sequence fasta file.
    '''
    def __init__(self,seq_dir):
        '''Creates FileNode and FileSystem for genome directory'''
        self.rootdir = fs_graph.FileNode(seq_dir)
        basename = self.rootdir.basename()

        '''Attributes for accession & file directory'''
        self.accession = basename[-len(basename):-5]
        
        '''Attributes for DNA sequence'''
        fa_file = open(seq_dir)
        fasta = ''
        with fa_file as f:
            lines = f.readlines()
            header = lines.pop(0)
            for l in lines:
                s = l.rstrip('\n')
                fasta = fasta + s
        fa_file.close()
        self.seq = fasta
"""

class Feature:
    '''
    An object that encapsulates a region of a genome sequence and the attributes
    associated with that region. Identified by a unique ID attribute.
    And associated with the attributes that must be provided for all features
    in a GFF file.
    '''

    def __init__(self,gff_line):
        
        '''Gets id'''
        if (gff_line.count("CDS") >= 1 or gff_line.count("gene") >= 1) and gff_line.count("ID=") >= 1:
            self.id = gff_line[gff_line.find("ID=")+3:gff_line.find(";")]
        elif gff_line.count("Protein") >= 1 and gff_line.count("protein_id"):
            self.id = gff_line[gff_line.find("protein_id=")+11:len(gff_line)-1]
        else:
            self.id = "other_id_" + str(SeqRecord.id_count)
        
        '''Gets type'''
        if gff_line.count("\tCDS\t") == 1:
            self.type = "CDS"
        elif gff_line.count("\tgene\t") == 1:
            self.type = "gene"
        elif gff_line.count("\tProtein\t") == 1:
            self.type = "Protein"
        elif gff_line.count("\tRegion\t") == 1:
            self.type = "Region"
        elif gff_line.count("\ttRNA\t") == 1:
            self.type = "tRNA"
        elif gff_line.count("\tncRNA\t") == 1:
            self.type = "ncRNA"
        elif gff_line.count("\trRNA\t") == 1:
            self.type = "rRNA"
        elif gff_line.count("\trepeat_region\t") == 1:
            self.type = "repeat_region"   
            
        self.read = self.remove_attributes(gff_line)
        
    def remove_attributes(self,gff_line):
        '''Removes \n and additional attributes'''
        gff_line = gff_line.replace('\n','')
        gff_line = gff_line.split("\t")
        gff_line = gff_line[0:8]
        gff_line = '\t'.join(gff_line)
        return gff_line
    
    def short(the_id,main):
        main_list = main.read.split("\t")
        short_print = the_id + "\t" + "[" + main_list[0] + "]" + "\t"+ \
            "\t" + main_list[2] + "\t" + main_list[3] + "\t" + \
            main_list[4] + "\t" + main_list[6]
        print(short_print + "\n")
    
    def fasta(the_id,sequence,main):
        main_list = main.read.split("\t")
        fasta_print = "> " + the_id + "\t" + "[" + main_list[0] + "]" + "\t"+ \
            "\t" + main_list[2] + "\t" + main_list[3] + "\t" + \
            main_list[4] + "\t" + main_list[6]
        print(fasta_print)
        print(sequence[int(main_list[3])-1:int(main_list[4])] + "\n")
    
    def full(the_id,main,optional):
        main_list = main.read.split("\t")
        opt_list = optional.read.split(";")
        print("ID: ",the_id)
        print("[Name: " + str(main_list[0]) + "]")
        print("Type: ",main_list[2])
        print("Start: ",main_list[3])
        print("Stop: ",main_list[4])
        print("Strand: ",main_list[6])
        for i in range(len(opt_list)):
            opt_list_2 = opt_list[i].split("=")
            print(opt_list_2[0] + ": " + opt_list_2[1])
        print("\n")
        
    def draw(self,ax,ypos,height=1,color='tab:blue',label=None):
        '''Draw the graph by for feature objects
        '''
        if self.read.split()[6] == '-':
            #draws arrow
            f_list = self.read.split()
            x=int(f_list[3])
            y=ypos
            dx=int(f_list[4]) - int(f_list[3])
            dy=0
            
            x = x+dx
            dx = -1*dx
            
            arrow = mpatches.Arrow(x,y,dx,dy,width=2*height,color=color)
            ax.add_patch(arrow)   
            
            #draws label
            if label is None:
                label=str(self.id)
            ax.text(x/2+(x+dx)/2,y,label,horizontalalignment='center',verticalalignment='center',color='white')
        elif self.read.split()[6] == '+':
            #draws arrow
            f_list = self.read.split()
            x=int(f_list[3])
            y=ypos
            dx=int(f_list[4]) - int(f_list[3])
            dy=0
            
            arrow = mpatches.Arrow(x,y,dx,dy,width=2*height,color=color)
            ax.add_patch(arrow)   
            
            #draws label
            if label is None:
                label=str(self.id)
            ax.text(x/2+(x+dx)/2,y,label,horizontalalignment='center',verticalalignment='center',color='white')

class Subtypes(Feature):
    '''
    Objects derived from the base Feature class that are specific to each unique
    feature type (e.g. GeneFeature,etc…).
    Associated with additional attributes specific to that subtype.
    '''
    def __init__(self,gff_line):   
        super().__init__(gff_line)
        self.read = self.remove_attributes(gff_line)
        
    def remove_attributes(self,gff_line):
        '''Removes main attributes and \n'''
        gff_line = gff_line.replace('\n','')
        gff_line = gff_line.split("\t")
        gff_line = gff_line[8]
        gff_line = ''.join(gff_line)
        return gff_line

#%%
print('\n')
if __name__ == '__main__':
    #----SeqRec Searching----#
    genome_directory = '/home/mcatolos/be500/datafiles/U00096.3'
    seqrec = SeqRecord(genome_directory)
    seqrec.parse_gff()
    filtered = seqrec.filter('type repeat_region')
    #filtered = seqrec.filter('fasta')
    #filtered = seqrec.filter('full')
    
    #----Plotting----#
    window_size = 10000
    gene = filtered.features[0].read.split()
    win_start = int(gene[3])-window_size
    win_stop = int(gene[4])+window_size
    
    search = 'type gene'
    filtered = seqrec.filter(search)
    
    fig,ax=plt.subplots(layout='constrained')
    #gene = filtered.features[0].read.split()
    plt.xlim([win_start,win_stop])
    plt.ylim([-.3,.3])

    for gene_feature in filtered.features:
        gene = gene_feature.read.split()
        if len(gene) > 1:
            gene_start = int(gene[3])
            gene_stop = int(gene[4])
            if gene_start >= win_start and gene_stop <= win_stop:
                
                y = 0
                gene_feature.draw(ax,y)
    plt.show()