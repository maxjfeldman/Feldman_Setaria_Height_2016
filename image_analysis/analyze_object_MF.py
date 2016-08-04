### Analyzes an object and outputs numeric properties

import cv2
import numpy as np
import plantcv as pcv
from . import print_image
from . import fatal_error

def analyze_object_MF(img,imgname,obj,mask,line_position,device,debug=False,filename=False):
  # Outputs numeric properties for an input object (contour or grouped contours)
  # Also color classification?
  # img = image object (most likely the original), color(RGB)
  # imgname= name of image
  # obj = single or grouped contour object
  # line_position = boundary line
  # device = device number. Used to count steps in the pipeline
  # debug= True/False. If True, print image
  # filename= False or image name. If defined print image
  device += 1
  ori_img=np.copy(img)
  if len(np.shape(img))==3:
    ix,iy,iz=np.shape(img)
  else:
    ix,iy=np.shape(img)
  # Change line postion to coordinate location on image
  line_position=int(ix)-int(line_position)
  # size is black and white image
  # size1 is dimensions of the image
  size = ix,iy,3
  size1 = ix,iy
  background = np.zeros(size, dtype=np.uint8)
  background1 = np.zeros(size1, dtype=np.uint8)
  background2 = np.zeros(size1, dtype=np.uint8)
  
  # Check is object is touching image boundaries (QC)
  frame_background = np.zeros(size1, dtype=np.uint8)
  frame=frame_background+1
  frame_contour,frame_heirarchy=cv2.findContours(frame,cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)
  ptest=[]
  vobj=np.vstack(obj)
  for i,c in enumerate(vobj):
      xy=tuple(c)
      pptest=cv2.pointPolygonTest(frame_contour[0],xy, measureDist=False)
      ptest.append(pptest)
  in_bounds=all(c==1 for c in ptest)
    
  # Convex Hull
  hull = cv2.convexHull(obj)
  hull_vertices = len(hull)
  # Moments
  #  m = cv2.moments(obj)
  m = cv2.moments(mask, binaryImage=True)
  ## Properties
  # Area
  area = m['m00']
  
  if area:
    # Convex Hull area
    hull_area = cv2.contourArea(hull)
    # Solidity
    solidity = 1
    if int(hull_area) != 0:
      solidity = area / hull_area
    # Perimeter
    perimeter = cv2.arcLength(obj, closed=True)
    # x and y position (bottom left?) and extent x (width) and extent y (height)
    x,y,width,height = cv2.boundingRect(obj)
    # Centroid (center of mass x, center of mass y)
    cmx,cmy = (m['m10']/m['m00'], m['m01']/m['m00'])
    # Ellipse
    center, axes, angle = cv2.fitEllipse(obj)
    major_axis = np.argmax(axes)
    minor_axis = 1 - major_axis
    major_axis_length = axes[major_axis]
    minor_axis_length = axes[minor_axis]
    eccentricity = np.sqrt(1 - (axes[minor_axis]/axes[major_axis]) ** 2)
    
    #Longest Axis: line through center of mass and point on the convex hull that is furthest away
    cv2.circle(background, (int(cmx),int(cmy)), 4, (255,255,255),-1)
    center_p = cv2.cvtColor(background, cv2.COLOR_BGR2GRAY)
    ret,centerp_binary = cv2.threshold(center_p, 0, 255, cv2.THRESH_BINARY)
    centerpoint,cpoint_h = cv2.findContours(centerp_binary,cv2.RETR_TREE,cv2.CHAIN_APPROX_NONE)

    
    dist=[]
    vhull=np.vstack(hull)
    
    for i,c in enumerate(vhull):
      xy=tuple(c)
      pptest=cv2.pointPolygonTest(centerpoint[0],xy, measureDist=True)
      dist.append(pptest)
    
    abs_dist=np.absolute(dist)
    max_i=np.argmax(abs_dist)
    
    caliper_max_x, caliper_max_y=list(tuple(vhull[max_i]))
    caliper_mid_x, caliper_mid_y=[int(cmx),int(cmy)]

    xdiff = float(caliper_max_x-caliper_mid_x)
    ydiff= float(caliper_max_y-caliper_mid_y)
    
    if xdiff!=0: 
      slope=(float(ydiff/xdiff))
    if xdiff==0:
      slope=1
    b_line=caliper_mid_y-(slope*caliper_mid_x)
    
    if slope==0:
      xintercept=0
      xintercept1=0
      yintercept='none'
      yintercept1='none'
      cv2.line(background1,(iy,caliper_mid_y),(0,caliper_mid_y),(255),1)
    else:
      xintercept=int(-b_line/slope)
      xintercept1=int((ix-b_line)/slope)
      yintercept='none'
      yintercept1='none'
      if 0<=xintercept<=iy and 0<=xintercept1<=iy:
        cv2.line(background1,(xintercept1,ix),(xintercept,0),(255),1)
      elif xintercept<0 or xintercept>iy or xintercept1<0 or xintercept1>iy:
        if xintercept<0 and 0<=xintercept1<=iy:
          yintercept=int(b_line)
          cv2.line(background1,(0,yintercept),(xintercept1,ix),(255),1)
        elif xintercept>iy and 0<=xintercept1<=iy:
          yintercept1=int((slope*iy)+b_line)
          cv2.line(background1,(iy,yintercept1),(xintercept1,ix),(255),1)          
        elif 0<=xintercept<=iy and xintercept1<0:          
          yintercept=int(b_line)
          cv2.line(background1,(0,yintercept),(xintercept,0),(255),1)          
        elif 0<=xintercept<=iy and xintercept1>iy:
          yintercept1=int((slope*iy)+b_line)
          cv2.line(background1,(iy,yintercept1),(xintercept,0),(255),1)          
        else:  
          yintercept=int(b_line)
          yintercept1=int((slope*iy)+b_line)
          cv2.line(background1,(0,yintercept),(iy,yintercept1),(255),1)
    
    ret1,line_binary = cv2.threshold(background1, 0, 255, cv2.THRESH_BINARY)
    #print_image(line_binary,(str(device)+'_caliperfit.png'))

    cv2.drawContours(background2, [hull], -1, (255), -1)
    ret2,hullp_binary = cv2.threshold(background2, 0, 255, cv2.THRESH_BINARY)
    #print_image(hullp_binary,(str(device)+'_hull.png'))
    
    caliper=cv2.multiply(line_binary,hullp_binary)    
    #print_image(caliper,(str(device)+'_caliperlength.png'))
    
    caliper_y,caliper_x=np.array(caliper.nonzero())
    caliper_matrix=np.vstack((caliper_x,caliper_y))
    caliper_transpose=np.transpose(caliper_matrix)
    caliper_length=len(caliper_transpose)

    caliper_transpose1 = np.lexsort((caliper_y, caliper_x))
    caliper_transpose2 = [(caliper_x[i],caliper_y[i]) for i in caliper_transpose1]
    caliper_transpose=np.array(caliper_transpose2)
    
    
    ##### Measure Canopy Height
    height_ab = line_position - y

    # If height is greater than 20 pixels make 20 increments (5% intervals)
    if height_ab >= 20:
      inc = height_ab /20

      # Define variable for max points and min points
      pts_max = []
      pts_min = []
      # Get max and min points for each of the intervals
      for i in range(1,20):
        if (i == 1):
          pt_max = y
          pt_min = y + (inc * i)
        else:
          pt_max = y + (inc * (i-1))
          pt_min = y + (inc * i)
        # Put these in an array
        pts_max.append(pt_max)
        pts_min.append(pt_min)

      # Combine max and min into a set of tuples
      point_range = list(zip(pts_max,pts_min))

      # define some list variables to fill
      row_median=[]
      row_ave=[]
      max_width=[]

      # For each of the 20 intervals
      for pt in point_range:
        # Get the lower and upper bounds  (lower and higher in terms of value; low point is actually towards top of photo, higher is lower of photo)
        low_point, high_point = pt
        # Get all rows within these two points
        rows=[]
        # Get a continuous list of the values between the top and the bottom of the interval save as vals
        vals = list(range(low_point, high_point))
        # For each row... get all coordinates from object contour that match row
        for v in vals:
          # Value is all entries that match the row
          value = obj[v == obj[:,0,1]]
          if len(value) > 0:
            # Could potentially be more than two points in all contour in each pixel row
            # Grab largest x coordinate (column)
            largest = value[:,0,0].max()
            # Grab smallest x coordinate (column)
            smallest = value[:,0,0].min()
            # Take the difference between the two (this is how far across the object is on this plane)
            row_width = largest - smallest
            # Append this value to a list
            rows.append(row_width)
          if len(value) == 0:
            row_width = 1
            rows.append(row_width)
        # For each of the points find the median and average width
        row_median.append(np.median(np.array(rows)))
        row_ave.append(np.mean(np.array(rows)))
        max_width.append(np.max(np.array(rows)))


      # Get the indicie of the largest median/average x-axis value (if there is a tie it takes largest index)
      indice_median = row_median.index(max(row_median))
      indice_ave = row_ave.index(max(row_ave))
      median_value = row_median[indice_median]
      ave_value = row_ave[indice_ave]
      max_value = max_width[indice_ave]

      # Canopy height as the height at which the average pixel width across a scoring window is maximized
      indice_reported = point_range[indice_ave]

      # Now you can get indice of point_range and make plots (lower and higher in terms of value; low point is actually towards top of photo, higher is lower of photo)

      lp, hp = indice_reported
      # Report Canopy height as the average of the upper and lower values of the scoring window (lower and higher in terms of value; low point is actually towards top of photo, higher is lower of photo)
      canopy_height = (lp + hp)/2
      # Define canopy width as mean value across scoring window
      canopy_width = ave_value
      # Find the center of the window
      w_center = (x+(x+int(max_value)))/2

      each_side = ave_value/2
      lside = w_center - each_side
      rside = w_center + each_side

      # Make rectangle and draw line at canopy height
  
      img_copy=np.copy(img)
      ix,iy,iz=np.shape(img)

      size = ix,iy,3
      size1 = ix,iy
      background = np.zeros(size, dtype=np.uint8)
      background_ch = np.zeros(size1, dtype=np.uint8)

      # Make a gray scale image with color area masked out
      gray = cv2.cvtColor(img,cv2.COLOR_RGB2GRAY)
      gray_rgb = cv2.cvtColor(gray,cv2.COLOR_GRAY2RGB)
      # Fill the rectable to be totally black (-1)
      cv2.rectangle(gray_rgb, (x, lp), (x+int(max_value), hp), (0,0,0), -1)


      cv2.rectangle(background_ch, (x, lp), (x+int(max_value), hp), (255,255,255), -1)
      device, masked_img = pcv.apply_mask(img_copy, background_ch, 'black', device, debug)

      # LOGICAL OR statement to combine background (gray_rgb) and the scoring window of canopy height (masked_img)
      device, example = pcv.logical_or(gray_rgb, masked_img, device, debug)
      # Draw scoring window rectangle
      cv2.rectangle(example, (x, lp), (x+int(max_value), hp), (255,0,0), 5)
      # Draw lines for height and max width in rectangle
      cv2.line(example, (int(lside), int(canopy_height)), (int(rside), int(canopy_height)), (0,0,255), 3)
      cv2.line(example, (int(w_center), int(canopy_height)), (int(w_center), int(y+height)), (0,0,255), 3)

      # Print image
      #cv2.imwrite('example_img.png', example)
  
    # If height is < 20 pixels Get widest point and report  
    if height_ab < 20:
      #rows=[]
      # Get a continuous list of the values between the top and the bottom of the interval save as vals
      #vals = list(range(y, y+height_ab))
      # For each row... get all coordinates from object contour that match row
      #for v in vals:
        # Value is all entries that match the row
        #value = obj[v == obj[:,0,1]]
        # Could potentially be more than two points in all contour in each pixel row
        # Grab largest x coordinate (column)
        #largest = value[:,0,0].max()
        # Grab smallest x coordinate (column)
        #smallest = value[:,0,0].min()
        # Take the difference between the two (this is how far across the object is on this plane)
        #row_width = largest - smallest
        # Append this value to a list
        #rows.append(row_width)
    
      # For each of the points find the median and average width  (if there is a tie it takes largest index)
      #max_width = np.max(np.array(rows))
      #canopy_height = rows.index(max(rows))
      #canopy_width = max_width
      max_width = width
      canopy_height = height
      #else:
      #  hull_area, solidity, perimeter, width, height, cmx, cmy = 'ND', 'ND', 'ND', 'ND', 'ND', 'ND', 'ND'
      
    # Change the values to reflect actual measurments not just point coordinates
    canopy_height = line_position - canopy_height
    centroid_y = line_position - cmy
    ellipse_y = line_position - center[1]
        
  #Store Shape Data
  shape_header=(
    'HEADER_SHAPES',
    'area',
    'hull-area',
    'solidity',
    'perimeter',
    'width',
    'height',
    'longest_axis',
    'center-of-mass-x',
    'center-of-mass-y',
    'hull_vertices',
    'in_bounds',
    'ellipse_center_x',
    'ellipse_center_y',
    'ellipse_major_axis',
    'ellipse_minor_axis',
    'ellipse_angle',
    'ellipse_eccentricity',
    'canopy_height',
    'canopy_width'
    )

  shape_data = (
    'SHAPES_DATA',
    area,
    hull_area,
    solidity,
    perimeter,
    width,
    height,
    caliper_length,
    cmx,
    centroid_y,
    hull_vertices,
    in_bounds,
    center[0],
    ellipse_y,
    major_axis_length,
    minor_axis_length,
    angle,
    eccentricity,
    canopy_height,
    canopy_width
    )

  analysis_images = []
      
   #Draw properties
  if area and filename:
    cv2.drawContours(ori_img, obj, -1, (255,0,0), 1)
    cv2.drawContours(ori_img, [hull], -1, (0,0,255), 1)
    cv2.line(ori_img, (x,y), (x+width,y), (0,0,255), 1)
    cv2.line(ori_img, (int(cmx),y), (int(cmx),y+height), (0,0,255), 1)
    cv2.line(ori_img,(tuple(caliper_transpose[caliper_length-1])),(tuple(caliper_transpose[0])),(0,0,255),1)
    cv2.circle(ori_img, (int(cmx),int(cmy)), 10, (0,0,255), 1)
    # Output images with convex hull, extent x and y
    extention = filename.split('.')[-1]
    #out_file = str(filename[0:-4]) + '_shapes.' + extention
    out_file = str(filename[0:-4]) + '_shapes.jpg'
    out_file1 = str(filename[0:-4]) + '_mask.jpg'
    
    print_image(ori_img, out_file)
    analysis_images.append(['IMAGE', 'shapes', out_file])
    
    print_image(mask,out_file1)
    analysis_images.append(['IMAGE','mask',out_file1])
    
  else:
    pass
  
  if debug:
    cv2.drawContours(ori_img, obj, -1, (255,0,0), 1)
    cv2.drawContours(ori_img, [hull], -1, (0,0,255), 1)
    cv2.line(ori_img, (x,y), (x+width,y), (0,0,255), 1)
    cv2.line(ori_img, (int(cmx),y), (int(cmx),y+height), (0,0,255), 1)
    cv2.circle(ori_img, (int(cmx),int(cmy)), 10, (0,0,255), 1)
    cv2.line(ori_img,(tuple(caliper_transpose[caliper_length-1])),(tuple(caliper_transpose[0])),(0,0,255),1)
    print_image(ori_img,(str(device)+'_shapes.jpg'))
 
  return device, shape_header, shape_data, analysis_images