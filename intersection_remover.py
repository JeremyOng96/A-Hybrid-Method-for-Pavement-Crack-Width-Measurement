#!/usr/bin/env python
# coding: utf-8

# # Intersection remover

# In[1]:


import cv2, glob, os, sknw, math
import numpy as np
from skimage import morphology,feature
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from rasterizing import *
from DSE_prune import *


# In[28]:


def remove_intersection(binary_image,skl,boundary_coords,graph, verbose=False):
    """
    This method removes the intersection in a crack binary map after the pruning step.
    
    (1) Obtain the nodes of the graph
    (2) Identify the degree corresponding to each node
    (3) Use midpoint circle algorithm to identify appropriate radius. The appropriate radius can be found when the
        number of intersections exceeds two times the number of degree.
    (4) Apply a circle dilation and remove the intersection given the radius.
    
    Argument:
    binary_img - The binary image of a crack
    skl - Skeleton image
    boundary_coords - Boundary of the crack
    graph - A graph after the pruning process
    
    Return:
    skeleton - An updated skeleton
    binary_img - An updated binary image
    """
    binary_img = binary_image.copy()*1
    skeleton_img = skl.copy()*1
    nodes = graph.nodes()
    index = list(nodes)
    degrees = [graph.degree(node) for node in nodes if graph.degree(node) > 2]
    points = [nodes[index[i]]['o'].squeeze() for i in range(len(index)) if graph.degree(index[i]) > 2]
    
    
    if verbose:
        print(f"The nodes are {nodes}")
        print(f"The degree corresponding to each node is {degrees}")
        print(f"The points that have deg more than 2 are {points}")
    if len(points) > 0:
        for i,point in enumerate(points):
            r,c = point.astype(np.int)
            radius = 0
            while True:
                circle_points = set(midpoint_circle(r,c,radius))
                intersections = list(circle_points & set(boundary_coords)) # Gets the intersection
                if len(intersections) >= 2*degrees[i]:
                    # morphology dilation to remove intersection with radius at center (r,c)
                    radius += 1 # This is done to clear some surroundings only, can remove if not wanted
                    binary_img[r-radius:r+radius+1,c-radius:c+radius+1] -= morphology.disk(radius)
                    skeleton_img[r-radius:r+radius+1,c-radius:c+radius+1] -= morphology.disk(radius) * 1
                    break
                else:
                    radius += 1
            if verbose:
                print(f"The radius used is:{radius}")
    return (binary_img == 255).astype(np.uint8) , (skeleton_img > 0).astype(np.uint8)
    


# # Test case

# In[29]:


if __name__ == "__main__":
    folder_path = "./data_for_crack_width_measurement/field_data_2"
    mask_paths = sorted(glob.glob(os.path.join(folder_path,'*png')))
    mask_files = []
    padded_mask = []
    for path in mask_paths:
        mask = cv2.imread(path,cv2.IMREAD_UNCHANGED)
        mask_files.append(mask)

    for mask in mask_files:
        pad_mask = np.pad(mask,((3,3),(3,3)),"constant")
        padded_mask.append(pad_mask)

    test_img = padded_mask[1]
    thin = morphology.thin(test_img)
    skel,dist = morphology.medial_axis(test_img,return_distance=True)
    skeleton_img, graph = skel_pruning_DSE(skel, dist, 70,return_graph=True)
    canny = feature.canny(test_img)
    
    fs = 100
    fig,ax = plt.subplots(nrows = 2, ncols = 2,figsize = (70,70), squeeze = True)
    plt.subplots_adjust(top = 0.99, bottom=0.01, hspace=-0.5, wspace=0.2)

    ax[0][0].imshow(test_img,cmap="Greys_r")
    ax[0][0].set_title(f"Binary Image",fontsize=fs)

    ax[0][1].imshow(skel,cmap="Greys_r")
    ax[0][1].set_title(f"Medial axis",fontsize=fs)

    ax[1][0].imshow(thin,cmap="Greys_r")
    ax[1][0].set_title(f"Thinning",fontsize=fs)

    ax[1][1].imshow(skeleton_img,cmap="Greys_r")
    ax[1][1].set_title(f"DSE",fontsize=fs)
    
    r,c = np.nonzero(canny)
    B = list(zip(r,c))
    b,s = remove_intersection(test_img,skeleton_img,B,graph,verbose=True)
    plt.figure()
    fig,ax = plt.subplots(nrows = 1, ncols = 2,figsize = (70,70), squeeze = True)
    ax[0].imshow(b,cmap="Greys_r")
    ax[0].set_title(f"Binary_map")
    
    ax[1].imshow(s,cmap="Greys_r")
    ax[1].set_title(f"Updated Skeleton")


# In[30]:


def visualize(img,center,r):
    im = np.reshape(img,(*img.shape,1)).astype(np.int)
    im = np.concatenate([im,im,im],axis=-1) * 255
    plt.figure()
    plt.imshow(im)
    plt.axis(False)
    color = (0,255,0)
    thickness = 1
    
    for i in range(len(center)):
        c = (int(center[i][1]),int(center[i][0]))
        im = cv2.circle(im, c, r[i], color, thickness)
        
    return im


# In[ ]:




