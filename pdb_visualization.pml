

load chrX_1Mb_coordinates.pdb
#load 3D_sklearn.pdb

util.performance(0)

bg_color white
show_as ribbon
spectrum

ray
png chrX_1Mb.png, dpi=300
#png 3D_sklearn.png, dpi=300