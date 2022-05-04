#!/usr/bin/env python
# coding: utf-8

# # Midpoint circle algorithm

# In[1]:


import numpy as np


# In[21]:


def midpoint_circle(r_c,c_c,r):
    """
    The midpoint circle algorithm.
    This algorithm calculates all points based on point (0,0).
    """
    x = r
    y = 0
    P = 1-r
    points = [(x,y)]
    
    while x > y:
        y += 1
        
        if P<=0:
            P = P + 2*y + 1
        else:
            x -= 1
            P = P + 2*y - 2*x + 1 
            
        points.append((x,y))
            
    # reflect on y = x 
    for point in points:
        temp_x = point[1]
        temp_y = point[0]
        temp_point = (temp_x,temp_y)
        if temp_point not in points:
            points.append(temp_point)
        
    # reflect on y axis
    for point in points:
        temp_x = -point[0]
        temp_y = point[1]
        temp_point = (temp_x,temp_y)
        if temp_point not in points:
            points.append(temp_point)
            
    # reflect on x axis
    for point in points:
        temp_x = point[0]
        temp_y = -point[1]
        temp_point = (temp_x,temp_y)
        if temp_point not in points:
            points.append(temp_point)
             
    points = [(j,i) for (i,j) in points] # invert points
    points = np.asarray(points) + np.asarray([r_c,c_c]) # translate back to origin
    points = [(point[0], point[1]) for point in points] # convert to set 
    return points


# In[24]:


def bresenham(pt1,pt2):
    y1,x1 = pt1
    y2,x2 = pt2
    m_new = 2 * (y2 - y1)
    slope_error_new = m_new - (x2 - x1)
    pts = []
    y=y1
    if y1 == y2:
        for x in range(min(x1,x2),max(x1,x2)+1):
            pts.append((y,x))
    elif x1 == x2:
        for y_t in range(min(y1,y2),max(y1,y2)+1):
            pts.append((y_t,x1))
    else:
        for x in range(min(x1,x2),max(x1,x2)+1):
            pts.append((y,x))
            slope_error_new =slope_error_new + m_new

            if (slope_error_new >= 0):
                y=y+1
                slope_error_new =slope_error_new - 2 * (x2 - x1)
            
    return pts

# # Test case

# In[22]:


if __name__ == "__main__":
    x = midpoint_circle(2,3,5)
    print(x)


# In[ ]:




