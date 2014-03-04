#!/usr/bin/env python

# Author : Jean Cury
"""Script representing genes of .lst files in region from line_deb included to line_fin included. All genome represented with line_deb = 1 and line_fin = -1. See complete options with -h option.
"""

#######
# For further versions :
# - tree is midpoint rooted by default : allow out-group rooting.
# - set the PATH_TO_LST without having to modify this script.
# - Choose which annotation you want (between gene name (default), or full annotation)
# - <write here if you want to add something>
#######

############################## Default variable ##############################

PATH_TO_LST = "/net/abigfour/gembases/Prokaryotes_1113a/LSTINFO/"

tree_mode = False
annotations = False
circ= False

############################## Import ##############################
import sys
import re 

try:
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	mpl.rcParams['pdf.fonttype'] = 42 # Allows text to be read as text in Illustrator
except ImportError:
	print ">>>>>> Need Matplotlib 1.3.0"
	sys.exit(0)
try:
	import numpy as np
except ImportError:
	print ">>>>>> Need Numpy 1.7.1"
	sys.exit(0)


############################## Functions ##############################
def usage():
	print u"""        >>>> Basic use <<<<
	
	python plot_LST_vXX.py id f l

	where 
	f is int, first line in .lst file
	l is int, last line in .lst file
	id being the name of a replicon in the gembase format : ggee.000.c00.00

	If an id.config file exist, it will color the genes as many colors as there are classes in the .config file which is formatted as follow :
	GENE001c01 classe1
	GENE002c01 classe2
	GENE003c01 classe1

	>>>> with annotations <<<<

	python plot_LST_vXX.py id f l -a

	>>>> circular representation <<<<

	python plot_LST_vXX.py id f l -c
		
	>>>> On a tree <<<<

	python plot_LST_vXX.py id_tree -t
		
	with -t option, an id_tree.treeconfig file is needed, and is formatted as follow :

	id1 f1 l1
	id2 f2 l2 
	id3 f3 l3
	
	where idx are names of leaves in the necessary i_tree.tree file, in Newick format, like :
	((id1,id2),id3)
	id_tree is the name of the tree.
	To colorize genes, one .config file is needed where all desired genes from all classes are concatenated in it.

		"""

def uniq(seq):
    # Not order preserving
    return {}.fromkeys(seq).keys()

def open_file(file, split_line = "\n", tab = True, sep = "\t"):
	"""Opens a file a returns the content in a list of lines
	If tab=True, we consider the file being a table, so lines are split"""
	
	f_in = open(file, 'r')
	
	data = f_in.read().split(split_line)[:-1] # must exist a last blank line
	if tab :
		isfloat = re.compile("^[\+\-]*[\d]+\.?[\d]*$")
		data = [[float(i) if isfloat.match(i.strip()) else i for i in f.split(sep)] for f in data]
	f_in.close()
	return data	

def get_data(ID, line_deb, line_fin, config, PATH_TO_LST=PATH_TO_LST):
	"""Function which get data from lst files for gene between line_deb and line_fin.
	config is updated with genes not in config files."""

	if int(line_deb) * int(line_fin) == 0:
		print "line number cannot be 0 (zero). Line number is 1-based (from 1 to N)\n"
		sys.exit(0)
	
	line_deb, line_fin = int(line_deb)-1, int(line_fin)
	try:
		All_data = open_file(PATH_TO_LST + ID + ".lst", tab=1,sep="|")

	except IOError:
		print "No such file :", ID
		sys.exit(0)
		
	N_genes = len(All_data)-1
	All_data = All_data[line_deb:line_fin]

	size = [i[1] for i in open_file(PATH_TO_LST+ID+".inf", tab=1,sep="size:") if len(i)==2][0]	

	try:
		size = int(size.split()[0]) # for draft with many sizes
	except AttributeError:
		pass

	digit = re.compile("^[\+\-]*[\d]+$")
	dic_data = {}
	for j, alldat in enumerate(All_data):
		dic_data[alldat[0].split()[4]] = [
		int(d) if digit.match(d) else -1 if d=="C" else 1 if d=="D" else d 
		for d in alldat[0].split()[0:4]+alldat[0].split()[5:8]+alldat[1:]
										]
		if alldat[0].split()[4] not in config.keys():
			config[alldat[0].split()[4]] = 0

	return dic_data, config, N_genes, size
	
def plot_locus(dic_data, N_genes, size, treedata=[], tree=False, circ=False, annotations=False):
	"""Function which plot the genes of 1 locus. Informations are in dic_data. N_genes and size are used for circular representation. treedata is list containing where should be plot this locus on a tree graph. tree, circ and annotations cannot be used together. """
	g = []
	ax.set_axis_off()
	tmp = [dd.split("_") for dd in dic_data.keys()]    
	tmp.sort(key = lambda x:int(x[1]))
	keys_data = ["_".join(t) for t in tmp]
	deb = dic_data[keys_data[0]][0]
	fin = dic_data[keys_data[-1]][1]
	
	if circ:
		ax.set_theta_offset(np.pi/2+2*np.pi)
		ax.set_axis_off()
		ax.set_theta_direction(-1)
		ponderation = float(len(dic_data))/float(N_genes)
		if ponderation < 0.9:
			N_pi = 2 * ponderation * (1+10/len(dic_data))
		else:
			N_pi = 2
		for k in keys_data:
			g.append(ax.bar(dic_data[k][0]*N_pi*np.pi/(ponderation*size),dic_data[k][2], width=N_pi*np.pi*float(dic_data[k][1]-dic_data[k][0])/(ponderation*size), bottom = 3, color=col_classes[config[k]]))	
		
	elif tree:
		x_pos, y_pos,ratio = treedata[0],treedata[1], treedata[2]

		#/!\/!\ >> you can switch a locus according to the following condition : 		
		if "cas1" in [config[k] for k in keys_data[-3:]] : # If cas1 gene, in the 3 last genes of the locus : switch the locus. You can change that condition.
		
			deb = dic_data[keys_data[-1]][1]
			fin = dic_data[keys_data[0]][0]
			for k in keys_data[::-1]:
				g.append( ax.bar(x_pos+ratio*(deb-dic_data[k][1]),-0.25*dic_data[k][2],ratio*(dic_data[k][1]-dic_data[k][0]),color=col_classes[config[k]],bottom = y_pos))
		
		else: # else default drawing
			for k in keys_data:
				
				g.append( ax.bar(x_pos+ratio*(dic_data[k][0]-deb),0.25*dic_data[k][2],ratio*(dic_data[k][1]-dic_data[k][0]),color=col_classes[config[k]],bottom = y_pos))
		
			
	
	else:	
		if not annotations:
			for k in keys_data:
				g.append( ax.bar(dic_data[k][0],dic_data[k][2],dic_data[k][1]-dic_data[k][0],color=col_classes[config[k]]))	

		else:
			for k in keys_data:
				g.append(ax.bar(dic_data[k][0],dic_data[k][2],dic_data[k][1]-dic_data[k][0],color=col_classes[config[k]]))	
				 #dict(arrowstyle="->") xytext=(10,10), textcoords='offset points')

				if g[-1].patches[0].get_width()/float(len(dic_data[k][5])) < len(dic_data)*4.9-len(dic_data)*4.5*10/100.  : # 195 : seuil empirique pr que le text tienne dans la boite.
					ax.annotate(dic_data[k][5],xy=(dic_data[k][0]+(dic_data[k][1]-dic_data[k][0])/2.,dic_data[k][2]),arrowprops=dict(arrowstyle="->"),xytext=(dic_data[k][0],dic_data[k][2]*np.random.choice([3,3.5,4,4.5,5])), rotation = 90) # dic_data[k][5] is gene name, dic_data[k][8] is annotation.
					#annot[-1].set_visible(False)
				else:
					ax.annotate(dic_data[k][5],xy=(dic_data[k][0],dic_data[k][2]*0.5),arrowprops=None,xytext=(5,0), textcoords='offset points')

		ax.tick_params(top = 'off', left="off", right = 'off')
		ax.set_yticklabels(ax.get_yticklabels(), visible = False)
# 		ax.set_ybound( -5, 5)
# 		ax.set_xbound(deb-10, fin+10)
		
	if legend:
		n=[]
		for i in classes:
			n.append(ax.bar(0,0,zorder=0,color=col_classes[i]))
		plt.legend(n,classes,loc=(1.015,0.2))
	if not circ and not tree:
		ax.set_ybound( -5, 5)
	# 	ax.set_xbound(deb-10, fin+10)
	if not tree_mode:
		plt.show()
	
	return g

def no_tree_mode(dic_data, ID, N_genes, size, annotations = annotations, circ = circ):
	"""Function for basic use (see usage() function)"""
	
	global config, ax
	fig = plt.figure(figsize=(15,10))
	ax = fig.add_subplot(111, polar=circ)
	
	fig.suptitle(ID, fontsize = 14)
	
	
	plot_locus(dic_data, N_genes, size, annotations = annotations, circ = circ)
	
def tree_opt(id_tree, indiv=None, annotations = annotations, circ=circ):
	""" Function for tree usage. One might want to draw single locus from the .treeconfig file. To do so, one needs to call this function with indiv set to the line number (0-based, int), to have this locus represented out of a tree graph. You can represent it with full annotation or in a circular form. This pgm must have ran before in a python environment,"""
	global config, ax
	locus_instance = []
	dic_data, id_leaves = [], []
	N_genes, size = [], [] #size=genome; taille=locus
	data = open_file(id_tree + ".treeconfig", tab=1, sep=" ")
	
	if indiv != None:
		dic_data, config, N_genes, size = get_data( data[indiv][0], data[indiv][1], data[indiv][2], config )
		tmp = [dd.split("_") for dd in dic_data.keys()]    
		tmp.sort(key = lambda x:int(x[1]))
		keys_data = ["_".join(t) for t in tmp]
		deb = dic_data[keys_data[0]][0]
		fin = dic_data[keys_data[-1]][1]
		taille = fin - deb
		nbr_locus = len(dic_data)
		no_tree_mode(dic_data,data[indiv][0],nbr_locus,taille,annotations = annotations, circ=circ)
		plt.show()
	

	else:
		taille = []
		for id_taxon, line_deb, line_fin in data:
			gd = get_data(id_taxon, line_deb, line_fin, config)
			dic_data.append(gd[0])
			config = gd[1]
			id_leaves.append(id_taxon)
			N_genes.append(gd[2])
			size.append(gd[3])
			tmp = [dd.split("_") for dd in dic_data[-1].keys()]    
			tmp.sort(key = lambda x:int(x[1]))
			keys_data = ["_".join(t) for t in tmp]
			deb = dic_data[-1][keys_data[0]][0]
			fin = dic_data[-1][keys_data[-1]][1]
			taille.append(fin - deb)
			
		try:
			tree_correspondance = dict(open_file(ID+".treetab", tab=1, sep=" "))
		except IOError:
			tree_correspondance = {}
			for il in id_leaves:
				tree_correspondance[il] = il
						
		tree = Phylo.read(id_tree+".tree","newick")
		tree.root_at_midpoint()
		Phylo.draw(tree,do_show=False,show_confidence=False)
		fig = plt.gcf()
		ax = fig.axes[0]
# 		fig.set_figwidth(20) # useless whith plt.show(), dunno why. But might be good if 
# 		fig.set_figheight(10) # fig.save() is used.
		xmax = np.max([a.get_position()[0] for a in ax.texts])
		
		biggest_locus = np.max(taille)
		ratio = 2*xmax /float(biggest_locus)

		for a in ax.texts:
			a.set_x(xmax+10*xmax/100.)
	
		for i,j in enumerate(id_leaves):
			locus_instance.append(plot_locus(dic_data[i], N_genes[i], size[i], treedata = [ 2*xmax, [l.get_position()[1] for l in ax.texts if l.get_text().strip()==tree_correspondance[j]][0], ratio ] , tree = True))
	
		ax.set_xlim(ax.get_xlim()[0], ax.get_xlim()[1] + biggest_locus*ratio + (ax.get_xlim()[1]+biggest_locus*ratio)*25/100.)   
	
		plt.show()
		return dic_data,data,locus_instance
		
def classes_config(ID):
	"""Get info from config file and set color for each type"""
	try:
		config = dict(open_file(ID + ".config", tab=1,sep=" "))
		classes = uniq(config.values())
		legend = True
	except IOError:
		print ID + ".config file not found" 
		config = {}
		classes = []
		legend = False

	## Color classes

	cm = plt.get_cmap("jet") # You can change color her by changing the colormap (google colormap matplotlib)
	cgen = (cm(1.15*i/(len(classes))) for i in range(len(classes))) # 1.15 can be change to get other color. It's quite unpredictable though.
	
	col_classes = {}
	for s in classes:
		col_classes[s] = next(cgen)
	col_classes[0] = "lightgrey" # color for genes not in the config files.
	
	return config, col_classes, legend, classes

############################## Read options ##############################

if __name__ == "__main__" :

	if len(sys.argv) <= 2:
		if sys.argv[-1][:2] == "-h":
			usage()
			sys.exit(0)
		elif len(sys.argv) == 1:
			usage()
			sys.exit(0)
		else:
			print ">>>>>>> Wrong arguments \n\n\n", usage()
			sys.exit(0)
	elif len(sys.argv) == 3:
		# draw the locus on a tree.
		if sys.argv[-1][:2] == "-t":
			tree_mode = True
			ID = sys.argv[-2]
			try:
				from Bio import Phylo
			except ImportError:
				print ">>>>>> Need Biopython 1.63 with -t option"
				sys.exit(0)	

		else:
			print ">>>>>>>Wrong arguments \n\n", usage()
			sys.exit(0)
		
	elif len(sys.argv) == 4:
		ID = sys.argv[-3]
		line_deb = int(sys.argv[-2])
		line_fin = int(sys.argv[-1])
	elif len(sys.argv)>4:
		corr = len(sys.argv)-4
		ID = sys.argv[-3-corr]
		try:
			line_deb = int(sys.argv[-2-corr])
			line_fin = int(sys.argv[-1-corr])
		except ValueError:
			print ">>>>>>> 2nd and 3rd args need to be int\n\n", usage()
			sys.exit(0)
		while corr:
			if sys.argv[-corr][:2] == "-c" : 
				# Draw circular
				circ = True
			elif sys.argv[-corr][:2] == "-a" :
				# Annotate genes. Only without the -c option
				annotations = True
			else:
				print ">>>>>>> Unknown option\n\n", usage()
				sys.exit(0)
			corr-=1
		

	
	else:
		print ">>>>>>> Wrong args\n\n", usage()
		sys.exit(0)


############################## Main ##############################

	config, col_classes, legend, classes = classes_config(ID)
	
	if not tree_mode:
		
		dic_data, config, N_genes, size = get_data(ID, line_deb, line_fin, config)
		no_tree_mode(dic_data, ID, N_genes, size, annotations, circ)

	else:
		dic_data_tree, data_tree,locus_instance = tree_opt(ID)





