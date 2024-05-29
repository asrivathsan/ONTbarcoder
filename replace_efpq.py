import sys

with open(sys.argv[1]) as infile:
	with open(sys.argv[2],'w') as o:
		l=infile.readlines()
		for i,j in enumerate(l):
			if ">" in j:
			#	lens=[int(j.strip().split('_')[-3]),int(j.strip().split('_')[-2])]
			#	print lens
			#	if max([int(j.strip().split('_')[-3]),int(j.strip().split('_')[-2])])>1:
				seq=l[i+1].replace("E","A").replace("F","G").replace("Q","C").replace("P","T")
				o.write(j+seq)
