from qgis.core import QgsProject, QgsVectorLayer, QgsRasterLayer, QgsDistanceArea, QgsCoordinateReferenceSystem, QgsPointXY, QgsPoint, QgsCircle,QgsGeometry
from math import pi,acos,cos,sqrt,tan,hypot
# from Geometry3D import *
from osgeo import gdal as gd
import numpy as np
from datetime import datetime

distance = QgsDistanceArea()
# distance.setSourceCrs(crs)
distance.setEllipsoid('EPSG:3857') # To change depending on the project CRS
driver = gd.GetDriverByName('GTiff')

def getAllPointLayers():
    #Fetch the current loaded layer
    layers = [l for l in QgsProject.instance().mapLayers().values() if l.type() == QgsVectorLayer.VectorLayer]
    for layer in layers:
        feats = list(layer.getFeatures())
        if feats != []:
            geom = (feats)[0].geometry()
            if not(geom.wkbType() in [1,4]):
                layers.remove(layer)
        else:
            layers.remove(layer)
    return layers

def getAllRasterLayers():
    #Fetch the current loaded layer
    layers = [l for l in QgsProject.instance().mapLayers().values() if l.type() == QgsRasterLayer.RasterLayer]
    return layers
    
def calculateBeamAngle(lumen, candela):
    beamAngle_rad = acos(1-((lumen/candela)/(2*pi)))
    return(beamAngle_rad * (180/pi))

def calculateLumen(beamAngle, candela):
    beamAngle_rad = beamAngle / (180/pi)
    return(2*pi*candela*(1-cos(beamAngle_rad)))
    
def calculateCandela(lumen, beamAngle):
    beamAngle_rad = beamAngle / (180/pi)
    return(lumen/(2*pi*(1-cos(beamAngle_rad))))

def getDistance(point_1,point_2,height_1, height_2):
    dist = distance.measureLine(point_1,point_2)
    return sqrt(((height_1-height_2)**2) +(dist**2))

# def getLuxViewshedValue(target, src, coef, lay,src_height):
#     height_target = lay.dataProvider().sample(target,1)
#     if (height_target[1]):
#         dist = getDistance(target,src,height_target[0],src_height)
#     else:
#         dist = getDistance(target,src,0,src_height)
#     return (coef/(dist**2))

# def getLuxTargetedValue(target, src, point, coef, lay,src_height, vis_lay):
#     if vis_lay.sample(point,1)[0] == 0:
#         return 0
#     else:
#         return getLuxViewshedValue(target, src, coef, lay, src_height)
    
def getLuxViewshedValue(target,src,coef, lumen,beamSolidAngle,layer, src_height,px_size, x_step,y_step):
    # target_circle = QgsCircle(QgsPoint(target),px_size/2)
    distance = sqrt((target[0]-src[0])**2 + (target[1]-src[1])**2)
    # print("src:"+str(src[0]))
    vector_ts = ((1/distance)*(target[0]-src[0]),(1/distance)*(target[1]-src[1])) if (distance != 0) else (1,1)
    first_dif = QgsPointXY(target[0] - (px_size/2) * vector_ts[0],target[1] - (px_size/2) * vector_ts[1])
    second_dif = QgsPointXY(target[0] + (px_size/2) * vector_ts[0],target[1] + (px_size/2) * vector_ts[1])
    first_dif_height = layer.sample(first_dif,1)[0]
    second_dif_height = layer.sample(second_dif,1)[0]
    area_cov = x_step * hypot(abs(first_dif_height-second_dif_height),x_step)
    # print(target,src)
    # print(first_dif,second_dif)
    first_deriv = np.array([target[0] - (px_size/2)*vector_ts[0],target[1] - (px_size/2)*vector_ts[1],coef*first_dif_height + (1-coef)*second_dif_height])
    second_deriv = np.array([target[0] + (px_size/2)*vector_ts[0],target[1] + (px_size/2)*vector_ts[1],coef*second_dif_height + (1-coef)*first_dif_height])
    src_deriv = np.array([src[0],src[1],src_height])
    ba = first_deriv - src_deriv
    bc = second_deriv - src_deriv
    ba = ba/np.linalg.norm(ba)
    bc = bc/np.linalg.norm(bc)
    area_angle = np.arccos(np.dot(ba, bc))
    solid_area_angle = 2*pi*(1-cos(area_angle/2))
    # return (ba,bc)
    avg_lumen = lumen * solid_area_angle / beamSolidAngle
    return avg_lumen/area_cov

def getLuxTargetedValue(target, src, point, beamAngle, coef, lumen, beamSolidAngle, lay,src_height, px_size, x_step,y_step):
    distance = sqrt((target[0]-src[0])**2 + (target[1]-src[1])**2)
    vector_ts = ((1/distance)*(target[0]-src[0]),(1/distance)*(target[1]-src[1])) if (distance != 0) else (1,1)
    vector_ts = ((1/distance)*(target[0]-src[0]),(1/distance)*(target[1]-src[1])) if (distance != 0) else (1,1)
    first_dif = QgsPointXY(point[0] - (px_size/2) * vector_ts[0],point[1] - (px_size/2) * vector_ts[1])
    second_dif = QgsPointXY(point[0] + (px_size/2) * vector_ts[0],point[1] + (px_size/2) * vector_ts[1])
    first_dif_height = lay.sample(first_dif,1)[0]
    second_dif_height = lay.sample(second_dif,1)[0]
    t_deriv = np.array([target[0],target[1], 0.5*first_dif_height + 0.5*second_dif_height])
    p_deriv = np.array([point[0],point[1], 0.5*first_dif_height + 0.5*second_dif_height])
    src_deriv = np.array([src[0],src[1],src_height])
    st = t_deriv - src_deriv
    sp = p_deriv - src_deriv
    st = st/np.linalg.norm(st)
    sp = sp/np.linalg.norm(sp)
    area_angle = np.arccos(np.dot(st, sp))
    if area_angle > beamAngle:
        return 0
    else:
        area_cov = x_step * hypot(abs(first_dif_height-second_dif_height),x_step)
        first_deriv = np.array([target[0] - (px_size/2)*vector_ts[0],target[1] - (px_size/2)*vector_ts[1],coef*first_dif_height + (1-coef)*second_dif_height])
        second_deriv = np.array([target[0] + (px_size/2)*vector_ts[0],target[1] + (px_size/2)*vector_ts[1],coef*second_dif_height + (1-coef)*first_dif_height])
        ba = first_deriv - src_deriv
        bc = second_deriv - src_deriv
        ba = ba/np.linalg.norm(ba)
        bc = bc/np.linalg.norm(bc)
        area_angle = np.arccos(np.dot(ba, bc))
        solid_area_angle = 2*pi*(1-cos(area_angle/2))
        avg_lumen = lumen * solid_area_angle / beamSolidAngle
        return avg_lumen/area_cov


#Add color to the graph and send to Paul
# def getVisibilityValue(current_point, src_coord, src_height, features, relief_lay,avg_height):
#     line = QgsGeometry.fromPolyline([QgsPoint(current_point[0],current_point[1],relief_lay.sample(QgsPointXY(current_point[0],current_point[1]),1)[0]),QgsPoint(src_coord[0],src_coord[1],src_height)])
#     for feat in features:
#         if (line.intersects(feat.geometry())):
#             print(line.intersection(feat.geometry()))
#             print(feat.geometry())
#             return 0
#     return 1

def getVisibilityValue(current_point, src_coord, src_height, features, relief_lay,avg_height):
    # line = QgsGeometry.fromPolyline([QgsPoint(current_point[0],current_point[1],relief_lay.sample(QgsPointXY(current_point[0],current_point[1]),1)[0]),QgsPoint(src_coord[0],src_coord[1],src_height)])
    # line =Line(Point(src_coord[0],src_coord[1],src_height),Point(current_point[0],current_point[1],relief_lay.sample(QgsPointXY(current_point[0],current_point[1]),1)[0]))
    l,m,n=current_point[0] - src_coord[0],current_point[1] - src_coord[1],relief_lay.sample(QgsPointXY(current_point[0],current_point[1]),1)[0] - src_height
    for feat in features:
        centroid_height = relief_lay.sample(((feat.geometry()).centroid().asPoint()),1)[0] + avg_height
        x,y = ((l/n)*(centroid_height-src_height))+src_coord[0],((m/n)*(centroid_height-src_height))+src_coord[1]
        if feat.geometry().contains(QgsPointXY(x,y)):
        # if (line.intersects(feat.geometry())):
        #     print(line.intersection(feat.geometry()))
        #     print(feat.geometry())
            return 0
    return 1
    
def createViewshed(lay, vis_lay, px_size, src_coord, lumen, beamAngle, src_height, size_q_render, size_m_render):
    #Get current time to avoid overwriting file
    date_time = datetime.now().strftime('%m_%d_%Y_%H_%M_%S')
    fn="/home/rayan/Documents/project/tifRaster/Viewshed_raster_" + date_time + ".tif"
    # coef_calc = lumen/(pi*(tan(beamAngle/(180/pi))**2))
    total_solid_angle = 2*pi*(1-cos(beamAngle/(180/pi)))
    # geot=[src_coord.x()-0.02,0.0004,0, src_coord.y()+0.02,0,-0.0004]
    geot=[src_coord[0].x()-(size_m_render / 2.0),size_m_render / size_q_render,0, src_coord[0].y()+(size_m_render / 2.0),0,-(size_m_render / size_q_render)]
    ds = driver.Create(fn, xsize=size_q_render, ysize=size_q_render, bands=1, eType=gd.GDT_Float32)
    rasterBand = np.zeros((size_q_render, size_q_render))
    coef = geot[1]/px_size
    for x in range(0,size_q_render):
        for y in range(1,size_q_render + 1): 
        # for y in range(0,size_q_render):
            # current_point = QgsPointXY(geot[0]+ x*geot[1], geot[3] + y*geot[5])
            current_point = (geot[0]+ x*geot[1], geot[3] + y*geot[5])
            # rasterBand[x][y] = sum(map(lambda a: (getLuxViewshedValue(current_point,a,coef_calc,lay,src_height) if (vis_lay.sample(current_point,1)[0] != 0) else 0),src_coord))
            rasterBand[x][size_q_render - y] = sum(map(lambda a: (getLuxViewshedValue(current_point, a, coef, lumen, total_solid_angle, lay, src_height, px_size, geot[1], geot[5]) if (vis_lay.sample(QgsPointXY(current_point[0],current_point[1]),1)[0] != 0) else 0),src_coord))
            # print("(" + str(x) + ", " + str(y) + "): " + str(rasterBand[x][y]))
        print("(" + str(x) + "): " + str(rasterBand[x][size_q_render//2]))
    ds.GetRasterBand(1).WriteArray(rasterBand)
    ds.SetGeoTransform(geot)
    # newLayer=QgsRasterLayer(fn,)
    return fn

def createTargeted(relief_lay, vis_lay, px_size, slope_lay, target_coord, src_coord, lumen, beamAngle, src_height, size_q_render, size_m_render):
    #Get current time to avoid overwriting file
    date_time = datetime.now().strftime('%m_%d_%Y-%H:%M:%S')
    fn="/home/rayan/Documents/project/tifRaster/Targeted_raster_" + date_time + ".tif"
    # coef_calc = lumen/(pi*(tan(beamAngle/(180/pi))**2))
    beamAngle_rad = beamAngle/(180/pi)
    total_solid_angle = 2*pi*(1-cos(beamAngle_rad))
    # ds = driver.Create(fn, xsize=size_q_render, ysize=size_q_render, bands=len(target_coord), eType=gd.GDT_Float32)
    ds = driver.Create(fn, xsize=size_q_render, ysize=size_q_render, bands=1, eType=gd.GDT_Float32)
    # for i in range(0,len(target_coord)):
        # geot=[target_coord[i].x()-(size_m_render / 2.0),size_m_render / size_q_render,0, target_coord[i].y()+(size_m_render / 2.0),0,-(size_m_render / size_q_render)]
    geot=[target_coord[0].x()-(size_m_render / 2.0),size_m_render / size_q_render,0, target_coord[0].y()+(size_m_render / 2.0),0,-(size_m_render / size_q_render)]
    coef = geot[1]/px_size
    rasterBand = np.zeros((size_q_render, size_q_render))
    for x in range(0,size_q_render):
        for y in range(1,size_q_render + 1):
        # for y in range(0,size_q_render): 
            current_point = (geot[0]+ x*geot[1], geot[3] + y*geot[5])
            rasterBand[x][size_q_render - y] = getLuxTargetedValue(target_coord[0],src_coord[0], current_point,beamAngle_rad, coef, lumen, total_solid_angle, relief_lay, src_height, px_size, geot[1], geot[5]) if (vis_lay.sample(QgsPointXY(current_point[0],current_point[1]),1)[0] != 0) else 0
            # rasterBand[x][y] = getLuxTargetedValue(target_coord[0],src_coord[0], current_point,beamAngle_rad, coef, lumen, total_solid_angle, relief_lay, src_height, px_size, geot[1], geot[5]) if (vis_lay.sample(QgsPointXY(current_point[0],current_point[1]),1)[0] != 0) else 0
        print("(" + str(x) + "): " + str(rasterBand[x][size_q_render//2]))
    # ds.GetRasterBand(i+1).WriteArray(rasterBand)
    ds.GetRasterBand(1).WriteArray(rasterBand)
    ds.SetGeoTransform(geot)

def createVisibility(relief_lay, build_lay, vis_layer, src_coord, src_height, px_size, size_q_render, size_m_render, avg_height=20):
    #Get current time to avoid overwriting file
    date_time = datetime.now().strftime('%m_%d_%Y-%H:%M:%S')
    fn="/home/rayan/Documents/project/tifRaster/Targeted_raster_" + date_time + ".tif"
    geot=[src_coord.x()-(size_m_render / 2.0),size_m_render / size_q_render,0, src_coord.y()+(size_m_render / 2.0),0,-(size_m_render / size_q_render)]
    ds = driver.Create(fn, xsize=size_q_render, ysize=size_q_render, bands=1, eType=gd.GDT_Float32)
    rasterBand = np.zeros((size_q_render, size_q_render))
    features=list(build_lay.getFeatures())
    coef = geot[1]/px_size
    for x in range(1,size_q_render+1):
    # for x in range(0,size_q_render):
        # for y in range(1,size_q_render + 1): 
        for y in range(0,size_q_render):
            current_point = (geot[0]+ x*geot[1], geot[3] + y*geot[5])
            rasterBand[size_q_render - x][y] = getVisibilityValue(current_point, src_coord, src_height, features, relief_lay,avg_height) if (vis_layer.sample(QgsPointXY(current_point[0],current_point[1]),1)[0] != 0) else 0
            # rasterBand[x][y] = getVisibilityValue(current_point, src_coord, src_height, features, relief_lay,avg_height) if (vis_layer.sample(QgsPointXY(current_point[0],current_point[1]),1)[0] != 0) else 0
            # rasterBand[x][size_q_render - y] = getVisibilityValue(current_point, src_coord, src_height, features, relief_lay,avg_height) if (vis_layer.sample(QgsPointXY(current_point[0],current_point[1]),1)[0] != 0) else 0
        # print("(" + str(x) + "): " + str(rasterBand[x][size_q_render//2]))
        print("(" + str(x) + "): " + str(rasterBand[size_q_render - x][size_q_render//2]))
    ds.GetRasterBand(1).WriteArray(rasterBand)
    ds.SetGeoTransform(geot)