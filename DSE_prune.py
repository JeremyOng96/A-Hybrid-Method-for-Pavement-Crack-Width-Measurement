#!/usr/bin/env python
# coding: utf-8

# In[1]:


from skimage import morphology,feature
import numpy as np
import cv2, glob, os, sknw, math
import matplotlib.pyplot as plt
from skimage.draw import line


# In[2]:


def flatten(l):
    return [item for sublist in l for item in sublist]


# In[3]:


def recnstrc_by_disk(pts,dist_tr,reconstruct):
    """
    Attempt to reconstruct the binary image from the skeleton using distance transform
    
    Arguments:
    pts (list) - branch points
    dist_tr (nd) - Distance transform matrix
    reconstruct (nd) - reconstructed image only on branch
    
    Return:
    reconstruct (nd) - reconstructed image
    """
    reconstruct = reconstruct * 1
    for pt in pts:
        r, c = pt
        radius = math.ceil(dist_tr[r,c])
        stel = morphology.disk(radius)
        reconstruct[r-radius:r+radius+1,c-radius:c+radius+1] += stel
        
    return reconstruct


# In[4]:


def get_weight(recn,branch):
    """
    The weights are calculated using reconstruction pixel loss.
    i.e. the number of pixels loss when comparing reconstruction and original image
    """
    w = 0
    reconstruction = np.reshape(recn,-1)
    branch_image = np.reshape(branch,-1)
    
    for i in range(recn.shape[0]*recn.shape[1]):
        # Given a pixel in the reconstruction, if the branch gives the reconstruction pixel, then add its weight.
        # If reconstruction pixel > branch image then it means other branches have an effect on the pixel.
        # Therefore, the branch that provides to the reconstruction pixel obtains a higher weight.
        w += (reconstruction[i] > 0) ^ ((reconstruction[i] - branch_image[i]) > 0)
    
    return w


# In[5]:


def _remove_branch_by_DSE(G, recn, dist, max_px_weight, checked_terminal=set()):
    deg = dict(G.degree())
    terminal_points = [i for i, d in deg.items() if d == 1]
    edges = list(G.edges())
    for s, e in edges:
        branch_recn = np.zeros_like(recn, dtype=np.int32)
        if s == e:
            G.remove_edge(s, e)
            continue
        
        pts = flatten([[v] for v in G[s][e].values()])[0]['pts']
        branch_recn = recnstrc_by_disk(pts,dist,branch_recn)
        weight = get_weight(recn, branch_recn)

        if s in terminal_points:
            checked_terminal.add(s)
            if weight < max_px_weight:
                G.remove_node(s)
                recn = recn - branch_recn
        if e in terminal_points:
            checked_terminal.add(e)
            if weight < max_px_weight:
                G.remove_node(e)
                recn = recn - branch_recn

    return G, recn


# In[6]:


def _remove_mid_node(G):
    start_index = 0
    while True:
        nodes = [x for x in G.nodes() if G.degree(x) == 2]
        if len(nodes) == start_index:
            break
        i = nodes[start_index]
        connected_nodes = list(G[i])
        # assert len(nbs)==2, 'degree not match'
        if len(connected_nodes) != 2:
            start_index = start_index + 1
            continue

        edge1 = G[i][connected_nodes[0]][0]
        edge2 = G[i][connected_nodes[1]][0]

        s1, e1 = edge1['pts'][0], edge1['pts'][-1]
        s2, e2 = edge2['pts'][0], edge2['pts'][-1]
        dist = np.array(list(map(np.linalg.norm, [s1-s2, e1-e2, s1-e2, s2-e1])))
        if dist.argmin() == 0:
            line = np.concatenate([edge1['pts'][::-1], [G.nodes[i]['o'].astype(np.int32)], edge2['pts']], axis=0)
        elif dist.argmin() == 1:
            line = np.concatenate([edge1['pts'], [G.nodes[i]['o'].astype(np.int32)], edge2['pts'][::-1]], axis=0)
        elif dist.argmin() == 2:
            line = np.concatenate([edge2['pts'], [G.nodes[i]['o'].astype(np.int32)], edge1['pts']], axis=0)
        elif dist.argmin() == 3:
            line = np.concatenate([edge1['pts'], [G.nodes[i]['o'].astype(np.int32)], edge2['pts']], axis=0)

        G.add_edge(connected_nodes[0], connected_nodes[1], weight=edge1['weight']+edge2['weight'], pts=line)
        G.remove_node(i)
    return G


# In[7]:


def skel_pruning_DSE(skel, dist, min_area_px=100, return_graph=False):
    """Skeleton pruning using dse
    
    Arguments:
        skel {ndarray} -- skeleton obtained from skeletonization algorithm
        dist {ndarray} -- distance transfrom map
    
    Keyword Arguments:
        min_area_px {int} -- branch reconstruction weights, measured by pixel area. Branch reconstruction weights smaller than this threshold will be pruned. (default: {100})
        return_graph {bool} -- return graph

    Returns:
        ndarray -- pruned skeleton map
    """
    graph = sknw.build_sknw(skel, multi=True)
    dist = dist.astype(np.int32)
    graph = _remove_mid_node(graph)
    edges = list(set(graph.edges()))
    pts = []
    for s, e in edges:
        temp_pts = flatten([[v] for v in graph[s][e].values()])[0]['pts']
        pts.extend(temp_pts)

    recnstrc = np.zeros_like(dist, dtype=np.int32)
    recnstrc = recnstrc_by_disk(np.array(pts, dtype=np.int32), dist, recnstrc) 
    num_nodes = len(graph.nodes())
    checked_terminal = set()
    while True:
        # cannot combine with other pruning method because the reconstruction map is not updated in other approach
        graph, recnstrc = _remove_branch_by_DSE(graph, recnstrc, dist, min_area_px, checked_terminal=checked_terminal)
        if len(graph.nodes()) == num_nodes:
            break
        graph = _remove_mid_node(graph)
        num_nodes = len(graph.nodes())
    if return_graph:
        return graph2im(graph, skel.shape), graph
    else:
        return graph2im(graph, skel.shape)      
     


# In[8]:


def graph2im(graph, shape):
    mask = np.zeros(shape, dtype=np.bool)
    for s,e in graph.edges():
        vals = flatten([[v] for v in graph[s][e].values()])
        for val in vals:
            coords = val.get('pts')
            coords_1 = np.roll(coords, -1, axis=0)
            for i in range(len(coords)-1):
                rr, cc = line(*coords[i], *coords_1[i])
                mask[rr, cc] = True
            mask[tuple(graph.nodes[s]['pts'].T.tolist())] = True
            mask[tuple(graph.nodes[e]['pts'].T.tolist())] = True
    return mask


# # Test case

# In[ ]:


if __name__ == "__main__":
    folder_path = "../Crack500/validation"
    mask_paths = sorted(glob.glob(os.path.join(folder_path,'mask','*png')))
    mask_files = []
    for path in mask_paths:
        mask = cv2.imread(path,cv2.IMREAD_UNCHANGED)
        mask_files.append(mask)


    # padding 
    padded_mask = []
    for mask in mask_files:
        pad_mask = np.pad(mask,((2,2),(2,2)),"constant")
        padded_mask.append(pad_mask)

    test_img = padded_mask[80]
    thin = morphology.thin(test_img)
    skel,dist = morphology.medial_axis(test_img,return_distance=True)
    skeleton_img = skel_pruning_DSE(skel, dist, 100)
    canny = feature.canny(test_img)
    
    fs = 100
    fig,ax = plt.subplots(nrows = 2, ncols = 2,figsize = (70,70), squeeze = True)
    plt.subplots_adjust(top = 0.99, bottom=0.01, hspace=-0.7, wspace=0.2)

    ax[0][0].imshow(test_img,cmap="Greys_r")
    ax[0][0].set_title(f"Binary Image",fontsize=fs)

    ax[0][1].imshow(skel,cmap="Greys_r")
    ax[0][1].set_title(f"Medial axis",fontsize=fs)
    ax[0][1].axis("off")

    ax[1][0].imshow(thin,cmap="Greys_r")
    ax[1][0].set_title(f"Thinning",fontsize=fs)
    ax[1][0].axis("off")

    ax[1][1].imshow(skeleton_img,cmap="Greys_r")
    ax[1][1].set_title(f"DSE",fontsize=fs)

