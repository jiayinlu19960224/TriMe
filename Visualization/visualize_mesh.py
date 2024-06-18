"""
Python visualization code for TriMe++
"""

#Import relevant Python libraries
import numpy as np
import matplotlib.pyplot as plt

#Import edges end point ID's: ..._tria_bar_ids.txt
tria_bar_pt_ids=np.genfromtxt(r"/path/to/file_tria_bar_ids.txt")
#Import point coordinates: ..._xy_id.txt
pt_coords=np.genfromtxt(r"/path/to/file_xy_id.txt")


#Visualize mesh
plt.figure(figsize = (30, 30))
#Loop through edges
for i in range(0,len(tria_bar_pt_ids)):
    #Edge end point IDs
    id0=int(tria_bar_pt_ids[i,0])
    id1=int(tria_bar_pt_ids[i,1])
    #Edge end point coordinates
    xl=[pt_coords[id0,0],pt_coords[id1,0]]
    yl=[pt_coords[id0,1],pt_coords[id1,1]]
    #Plot edge
    plt.plot(xl,yl, "-", color="dodgerblue", alpha=1,linewidth=4) #lime
plt.axis('off')
plt.scatter(pt_coords[:,0],pt_coords[:,1],s=600,c="fuchsia")
plt.show() 




