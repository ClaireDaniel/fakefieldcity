#%%
from igraph import *
import pandas as pd
import numpy as np
from random import randint
import itertools

import matplotlib.pyplot as plt
import matplotlib.patches as patch



from shapely.geometry import Polygon
from shapely.ops import cascaded_union
from descartes import PolygonPatch
import matplotlib.collections as plt_c


#%%
class field:
  
    def __init__(self, x_holdings = 2, y_holdings = 2, x_blocks = 8, y_blocks = 8, block_length = 16, road_width = 15, lot_width = 15, lot_length = 30):
        
        #%%
        '''
        x_holdings = 2
        y_holdings = 2
        x_blocks = 8
        y_blocks = 8
        block_length = 8
        road_width = 15 
        lot_width = 15
        lot_length = 30
        '''
        #%%
        
        if (x_holdings%2 != 0) | (y_holdings%2 != 0):
            raise Exception('Land holdings input for each axis must be an even number')
    
        if(x_blocks/x_holdings%2 != 0) | (y_blocks/y_holdings%2 != 0):
            raise Exception('Block number input along each axis must be an even multiple of land holdings')
        
        if block_length not in [8,16]:
            raise Exception('Block length must be either 8 or 16')

        self.x_holdings = x_holdings
        self.y_holdings = y_holdings
        self.x_blocks = x_blocks
        self.y_blocks = y_blocks
        self.block_length = block_length
        self.road_width = road_width
        self.lot_width = lot_width
        self.lot_length = lot_length

        #%%
    #def create_field(self):
    
        # Number of roads along the x and the y axis
        x_roads = x_blocks + 1
        y_roads = y_blocks + 1
        
        # Full length of the x and the y axis
        x_length = (x_blocks*2 + x_roads)
        y_length = (y_blocks*block_length + y_roads)
        
        # Full number of nodes in the graph
        v_number = x_length * y_length
        
        #Initial trails y and y coordinates(public road reserve land) around perimeter of each initial land holding
        init_trail_y = list(range(0, y_length))[0::math.floor(y_length/y_holdings)]
        init_trail_x = list(range(0, x_length))[0::math.floor(x_length/x_holdings)]
        
        # x and y positions of all lots in holdings
        holding_xpos = list(set(list(range(0, x_length))) - set(init_trail_x))
        holding_ypos = list(set(list(range(0, y_length))) - set(init_trail_y))
        
        # Number of possible subdivisions of holdings along the x and y axis
        y_holding_length = int((y_length - len(init_trail_y))/y_holdings)
        x_holding_length = int((x_length - len(init_trail_x))/x_holdings)
        
        # List of y axis coordinate lists for lots in each land holding
        holding_ypos_chunks = [holding_ypos[y:y+y_holding_length] for y in range(0, len(holding_ypos), y_holding_length)]
        holding_xpos_chunks = [holding_xpos[x:x+x_holding_length] for x in range(0, len(holding_xpos), x_holding_length)]

        # Coordinates for the road placement along the x and y axes
        road_x = list(range(0, x_length))[0::3]
        road_y = list(range(0, y_length))[0::block_length+1]
        
        # X coordinates for blocks
        blocks_xpos = list(set(list(range(0, x_length))) - set(road_x))
        #Tuple for X axis neightbour relationship
        blocks_xpos_tuple = list(zip(*[iter(blocks_xpos)]*2))
        
        # Y coordinates for blocks
        blocks_ypos = list(set(list(range(0, y_length))) - set(road_y))
        # List of y axis coordinate lists for lots in each block   
        blocks_ypos_chunks = [blocks_ypos[x:x+block_length] for x in range(0, len(blocks_ypos), block_length)]
        # Tuple for Y axis lot neightbour relationship
        blocks_ypos_tuple = []
        
        for i in range(0, len(blocks_ypos_chunks)):
            for j in blocks_ypos_chunks[i][:-1]:
                tup = (j, j+1)
                blocks_ypos_tuple.append(tup)
        
        # Tuple list for relationships along the x axis - which lots ajoin a road
        block_road_tuple = []
        for i in range(0, x_length):
        
            if (i not in blocks_xpos) & (i != 0):
                block_road_tuple.append((i-1, i))
                
            if (i in road_x) & (i != x_length-1):
                block_road_tuple.append((i, i+1))
 
        # List of y axis coordinate lists for lots in each group of 8
        g8_ypos_chunks = [blocks_ypos[y:y+4] for y in range(0, len(blocks_ypos), 4)]
        xpos_chunks = [blocks_xpos[x:x+2] for x in range(0, len(blocks_xpos),2)]
        
        # List of y axis coordinate lists for lots in each group of 4
        y_g4_no = int(len(blocks_ypos)/2)
        g4_ypos_chunks = [blocks_ypos[y:y+2] for y in range(0, len(blocks_ypos), 2)]
        #%%

        # ---------------------------------------------------------------------------------------------------------------------
        # STEP 5 - Create graph
        # ---------------------------------------------------------------------------------------------------------------------
        
        g = Graph()
        g.add_vertices(v_number)
            
        # ---------------------------------------------------------------------------------------------------------------------
        # STEP 6 - Add vertex attributes
        # ---------------------------------------------------------------------------------------------------------------------
        
        #ID
        g.vs["id"] = list(range(0, v_number))
        
        # X and Y coordinates
        g.vs["x"] = list(range(0, x_length)) * y_length
        g.vs["y"] = np.repeat(list(range(0, y_length)), x_length)
        
        #Road starting attributes
        g.vs.select(x_in = road_x)["type"] = "road"
        g.vs.select(y_in = road_y)["type"] = "road"
        g.vs.select(x_in = road_x)["landuse"] = "UD"
        g.vs.select(y_in = road_y)["landuse"] = "UD"
        g.vs.select(x_in = init_trail_x)["landuse"] = "RR"
        g.vs.select(y_in = init_trail_y)["landuse"] = "RR"
        g.vs.select(x_in = init_trail_x)["lot_number"] = "RR"
        g.vs.select(y_in = init_trail_y)["lot_number"] = "RR"
        
        #Add attribute for development 
        for x in blocks_xpos:
            for y in blocks_ypos:
                g.vs.select(y=y, x = x)["type"] = "lot"
                g.vs.select(y=y, x = x)["landuse"] = "UD"
                
        #Assign holding number to individual lots
        holding_count = 0
        
        for y in range(0, y_holdings):
            ypos = holding_ypos_chunks[y]
                
            for x in range(0, x_holdings):  
                xpos = holding_xpos_chunks[x]  
        
                g.vs.select(y_in = ypos, x_in = xpos)["holding_number"] = holding_count
                g.vs.select(y_in = ypos, x_in = xpos)["lot_number"] = holding_count
                
                self.latest_lot_number = holding_count
                holding_count += 1
                
        #Assign block number to individual lots
        block_count = 0
        
        for y in range(0, y_blocks):
            ypos = blocks_ypos_chunks[y]
                
            for x in range(0, x_blocks):    
                xpos = list(blocks_xpos_tuple[x])  
                
                g.vs.select(y_in = ypos, x_in = xpos)["block_number"] = block_count
            
                block_count += 1
                
        
        #Assign group of 8 number to individual lots
        g8_count = 0
        
        for y in g8_ypos_chunks:
                
            for x in xpos_chunks:  
        
                g.vs.select(y_in = y, x_in = x)["g8_number"] = g8_count
        
                g8_count += 1
                
        #Assign group of 4 number to individual lots
        g4_count = 0
        
        for y in g4_ypos_chunks:       
            for x in xpos_chunks:  
        
                g.vs.select(y_in = y, x_in = x)["g4_number"] = g4_count
        
                g4_count += 1
        
        #Assign group of 2 number to individual lots
        g2_count = 0
        
        for y in blocks_ypos:       
            for x in xpos_chunks:  
        
                g.vs.select(y = y, x_in = x)["g2_number"] = g2_count
        
                g2_count += 1  
                
        #Assign coordinates to each to make polygons
    
        x_select = g.vs.select(y = 1)
        y_select = g.vs.select(x = 1)
        
        x_coords = [0]
        for x in x_select:
            if x['type'] == 'road':
                x_coords.append(x_coords[-1] + road_width)
            elif x['type'] == 'lot':
                x_coords.append(x_coords[-1] + lot_length)
                
        y_coords = [0]
        for y in y_select:
            if y['type'] == 'road':
                y_coords.append(y_coords[-1] - road_width)
            elif y['type'] == 'lot':
                y_coords.append(y_coords[-1] - lot_width)
                
        for v in g.vs():
        
            x = v['x']
            y = v['y']
            
            coords = [(x_coords[x], y_coords[y]),
                     (x_coords[x + 1], y_coords[y]),
                     (x_coords[x + 1], y_coords[y + 1]),
                     (x_coords[x], y_coords[y + 1])]
            
            v['coords'] = coords
        
        #Assign a lot number
        
        
        

        # ---------------------------------------------------------------------------------------------------------------------
        # STEP 7 - Assign vertex size for visualisation
        # ---------------------------------------------------------------------------------------------------------------------
        
        #Different size for blocks
        g.vs.select(type = "lot")["size"] = 10
        g.vs.select(type = "road")["size"] = 5
 
        # ---------------------------------------------------------------------------------------------------------------------
        # STEP 8 - Add edges to the network to define relationships between different components
        # ---------------------------------------------------------------------------------------------------------------------
        
        #Rear boundary neighbours
        for x in blocks_xpos_tuple:
            for y in blocks_ypos:
                v1 = g.vs.select(y=y, x = x[0])[0].index
                v2 = g.vs.select(y=y, x = x[1])[0].index        
                g.add_edge(v1, v2, type = 'rear neighbour')
        
        #Side boundary neighbours
        for y in blocks_ypos_tuple:
            for x in blocks_xpos:
                v1 = g.vs.select(y=y[0], x = x)[0].index
                v2 = g.vs.select(y=y[1], x = x)[0].index        
                g.add_edge(v1, v2, type = 'side neighbour')
        
        #Block to road links
        for x in block_road_tuple:
            for y in blocks_ypos:
                
                v1 = g.vs.select(y=y, x = x[0])[0].index
                v2 = g.vs.select(y=y, x = x[1])[0].index        
                g.add_edge(v1, v2, type = 'road to lot')
                #print(x,y)
                
        #Across road neighbours
        for i in range(0,len(blocks_xpos_tuple)-1):
            for y in blocks_ypos:
                v1 = g.vs.select(y=y, x = blocks_xpos_tuple[i][1])[0].index
                v2 = g.vs.select(y=y, x = blocks_xpos_tuple[i+1][0])[0].index
                g.add_edge(v1, v2, type = 'accross road neighbour', curved = True)
                
        #Road to road links
        for x in road_x:
            for y in range(0,y_length-1):
                v1 = g.vs.select(y=y, x = x)[0].index
                v2 = g.vs.select(y=y+1, x = x)[0].index
                g.add_edge(v1, v2, type = 'road to road', status = 'undev')
                
        for y in road_y:
            for x in range(0,x_length-1):
                v1 = g.vs.select(y=y, x = x)[0].index
                v2 = g.vs.select(y=y, x = x+1)[0].index
                g.add_edge(v1, v2, type = 'road to road', status = 'undev')
                
                
        #Road to block and accross neighbours on y axis
        for x in blocks_xpos:
            for y in road_y:
                
                if y != road_y[len(road_y)-1]:
                    v1 = g.vs.select(y= y, x = x)[0].index
                    v2 = g.vs.select(y= y+1, x = x)[0].index
                
                    g.add_edge(v1, v2, type = 'road to lot')
    
                if y != road_y[0]:
                    v1 = g.vs.select(y= y, x = x)[0].index
                    v2 = g.vs.select(y= y-1, x = x)[0].index
                    
                    g.add_edge(v1, v2, type = 'road to lot')
    
                if (y != road_y[len(road_y)-1]) & (y != road_y[0]):
                    
                    v1 = g.vs.select(y= y-1, x = x)[0].index
                    v2 = g.vs.select(y= y+1, x = x)[0].index
                    
                    g.add_edge(v1, v2, type = 'accross road neighbour')
                    
        #Road to block and accross neighbours on intersections
        
        for x in road_x:
    
            for y in road_y:
                
                if (y != road_y[-1]) & (x != road_x[-1]):
                
                    v1 = g.vs.select(y= y, x = x)[0].index
                    v2 = g.vs.select(y= y+1, x = x+1)[0].index
                        
                    g.add_edge(v1, v2, type = 'road to lot')
                    
                if (y != road_y[0]) & (x != road_x[0]):
                    
                    v1 = g.vs.select(y= y, x = x)[0].index
                    v2 = g.vs.select(y= y-1, x = x-1)[0].index
                        
                    g.add_edge(v1, v2, type = 'road to lot')
                    
                if (y != road_y[0]) &  (x != road_x[-1]):
                    
                    v1 = g.vs.select(y= y, x = x)[0].index
                    v2 = g.vs.select(y= y-1, x = x+1)[0].index
                        
                    g.add_edge(v1, v2, type = 'road to lot')
                    
                if (y != road_y[-1]) &  (x != road_x[0]):
                    
                    v1 = g.vs.select(y= y, x = x)[0].index
                    v2 = g.vs.select(y= y+1, x = x-1)[0].index
                        
                    g.add_edge(v1, v2, type = 'road to lot')
 
        # ---------------------------------------------------------------------------------------------------------------------
        # STEP 9 - Add color to edges
        # ---------------------------------------------------------------------------------------------------------------------
        
        col_dic = {'rear neighbour':'red', 'side neighbour':'red', 'road to lot':'grey', 'road to road':'black', 
                   'accross road neighbour':'pink'}
        
        col = []
        for i in g.es["type"]:
            col.append(col_dic[i])
        g.es['color'] = col  
        
        g.es.select(status = 'undev')['color'] = 'grey'
    
        # ---------------------------------------------------------------------------------------------------------------------
        # STEP 10 - Set up for automata model
        # ---------------------------------------------------------------------------------------------------------------------
        self.landuse_lst = ['R1', 'R2', 'C', 'I', 'P', 'UD', 'US', 'RR']
        self.dev_landuse_lst = ['R1', 'R2', 'C', 'I', 'P']
        
        self.landuse_name_dic = {'R1':'High Density Residential', 'R2':'Low Density Residential', 'C':'Commerical', 
                                 'I':'Industry', 'P':'Park', 'UD':'Undeveloped', 'US':'Unserviced', 'RR':'Road Reserve',}
        self.landuse_proportion_dic =  {'R1':0.07, 'R2':0.53, 'C':0.15, 'I':0.15, 'P':0.10}
        self.landuse_min_lots_dic= {'R1':4, 'R2':2, 'C':4, 'I':8, 'P':4}
        self.landuse_color_dic = {'R1':'red', 'R2':'pink', 'C':'blue', 'I':'purple', 'P':'green', 'UD':'grey', 'US':'brown', 'RR':'black'}
        
        
        # Add colour
        g.vs['frame_color'] = 'white'
        g.vs.select()['color'] = [self.landuse_color_dic[lu] for lu in g.vs.select()['landuse']]

        
        # ---------------------------------------------------------------------------------------------------------------------
        # STEP 12 - Create copy and remove neighbour edges in order to measure distance
        # ---------------------------------------------------------------------------------------------------------------------
        g_distance = g.copy()
        e_selection = g_distance.es.select(type_in = ['accross road neighbour', 'side neighbour', 'rear neighbour'])
        e_remove = [e.index for e in e_selection]
        g_distance.delete_edges(e_remove)
                
        
        # Set roads around the seed to developed
        
        #Select all undeveloped road noes directly attached to developed land - change to developed
        
        #Select all undeveloped road edges between each road node - change to developed
        

        self.y_length = y_length
        self.x_length = x_length
        self.x_coords = x_coords
        self.y_coords = y_coords
        self.g = g
        self.g_distance = g_distance


    def set_seed(self, block_1 = [('I', 2), ('P', 1), ('R2', 6)], block_2 = [('R2', 8), ('P', 1), ('R1', 1), ('C', 2)]):

        #Getting y coordinates for seed
        y_mid = int(self.y_length/2)+1
        ypos_finish = y_mid-1
        
        if self.block_length == 16:
            ypos_start = ypos_finish - self.block_length
            ypos_seed = list(range(ypos_start, ypos_finish))
        
        elif self.block_length == 8:
            ypos_start = ypos_finish - 2 * self.block_length - 1
            ypos_break = ypos_start + self.block_length
            ypos_restart = ypos_break + 1
        
            ypos_seed = list(range(ypos_start, ypos_break)) + list(range(ypos_restart, ypos_finish))
            
        #Getting x coordinates for seed
        x_mid = int(self.x_length/2)+1
        xpos_seed = [x_mid - 3, x_mid - 2, x_mid, x_mid + 1]
        
        #Develop the roads around the blocks first
        v_selection = self.g.vs.select(x_in = xpos_seed, y_in = ypos_seed)
        
        for v1 in v_selection:
            v_neighbours = v1.neighbors()
            
            for v2 in v_neighbours:
                if v2['type'] == 'road':
                    v2['landuse'] = 'RR'
                    v2['lot_number'] = 'RR'
        
        for block in [block_1, block_2]:
            
            if block == block_1:
                x = 0
            else:
                x = 2
                    
            shift = 0  
            
            for item in block:
                
                lu = item[0]

                for app_no in range(0,item[1]):
        
                    if lu == 'I':
                        g_number = self.g.vs.select(x = xpos_seed[x], y = ypos_seed[shift])['g8_number']
                        self.g.vs.select(g8_number = g_number[0])["landuse"] = lu
                        self.g.vs.select(g8_number = g_number[0])["lot_number"] = self.latest_lot_number + 1
        
                        self.latest_lot_number = self.latest_lot_number + 1
        
                        shift = shift + 4
        
                    if lu in ['P', 'C', 'R1']:
        
                        g_number = self.g.vs.select(x = xpos_seed[x], y = ypos_seed[shift])['g4_number']
                        self.g.vs.select(g4_number = g_number[0])["landuse"] = lu
                        self.g.vs.select(g4_number = g_number[0])["lot_number"] = self.latest_lot_number + 1
        
                        self.latest_lot_number = self.latest_lot_number + 1
        
                        shift = shift + 2
        
                    if lu == 'R2':
                        g_number = self.g.vs.select(x = xpos_seed[x], y = ypos_seed[shift])['g2_number']
                        self.g.vs.select(g2_number = g_number[0])["landuse"] = lu
                        self.g.vs.select(g2_number = g_number[0])["lot_number"] = [self.latest_lot_number + 1, self.latest_lot_number + 2]
        
                        self.latest_lot_number = self.latest_lot_number + 2
        
                        shift = shift + 1
            
        # Add metro station
        self.g.vs.select(x = x_mid, y = y_mid-2)["add_infras"] = 'Metro Station'
            
        self.g.vs.select()['color'] = [self.landuse_color_dic[lu] for lu in self.g.vs.select()['landuse']]
        self.g.vs.select(add_infras = 'Metro Station')['frame_color'] = 'orange'
        self.g.vs.select(add_infras = 'Metro Station')["frame_width"] = 3
          
    def required_development(self):

        summary_lu = self.g.vs.select(landuse_in = self.dev_landuse_lst)['landuse']
        lu_dic = {lu:summary_lu.count(lu) for lu in summary_lu}
        lu_df = pd.DataFrame.from_dict(lu_dic, orient = 'index', columns = ['current_lots'])
    
        #Calculate each landuse as proportion of total
        total = sum(lu_df.current_lots)
        lu_df['landuse_prop'] = lu_df.current_lots/total
    
        #Calculate the difference between the current proportion and the idea balance
        lu_df['assigned_prop'] = [self.landuse_proportion_dic[lu] for lu in lu_df.index]
        lu_df['prop_diff'] = lu_df.assigned_prop - lu_df.landuse_prop
    
        #Find the landuse which is most over supplied and undersupplied
        max_lu = lu_df[lu_df.prop_diff == max(lu_df.prop_diff)]
        min_lu = lu_df[lu_df.prop_diff == max(lu_df.prop_diff)]
    
        # Scale the other land uses accordingly
        lu_df['required_lots'] = lu_df.prop_diff * int(max_lu.current_lots[0]/max_lu.assigned_prop[0])
        lu_df['required_lots'] = [int(i) for i in lu_df.required_lots]
    
        #Adjust so each total can be met by the minimum lots requirement
        lu_df['min_lots'] = [self.landuse_min_lots_dic[lu] for lu in lu_df.index]
        lu_df['required_lots'] = lu_df.required_lots + lu_df.min_lots-(lu_df.required_lots%lu_df.min_lots)
        lu_df['app_no'] = [int(i) for i in lu_df.required_lots/lu_df.min_lots]
    
        #Number of blocks to develop
        blocks_required = sum(lu_df.required_lots)/(self.block_length*2)
        print(lu_df)
        print("Bocks required: ", blocks_required)
        
        self.required_dev_df = lu_df
    
    def make_application(self,lu):

        vs_metro = self.g.vs.select(add_infras = 'Metro Station')['id'] #Check this later
        dic_attributes = ['id', 'x', 'y', 'type', 'block_number', 'g8_number', 'g4_number', 'g2_number', 'holding_number','landuse', 'add_infras']
        vs_lu = self.g.vs.select(landuse = lu)
    
        app_df = pd.DataFrame(columns = ['vertex_from'] + dic_attributes + ['metro_distance'])
    
        for vf in vs_lu:
            v_from = vf.index
            v_to = vf.neighbors()
    
    
            for vt in v_to:
    
                if vt['type'] != 'road':
    
                    metro_distance = self.g_distance.shortest_paths_dijkstra(source = vt.index, target = vs_metro)
                    dic = vt.attributes()
                    row = [vf.index] + [dic[a] for a in dic_attributes] + min(metro_distance)
                    app_df.loc[len(app_df)] = row
    
        # Criteria - Undeveloped land closed to metro
        
        if len(app_df[app_df['landuse'] == 'UD']) > 0:
        
            sl_df = app_df[app_df['landuse'] == 'UD']
            min_distance = min(sl_df['metro_distance'])
            sl_df = sl_df[sl_df['metro_distance'] == min_distance]
            choice = sl_df.sample(n=1)
        
        #Otherwise take over some other land
        else:
            sl_df = app_df[app_df['landuse'] == 'R2']
            min_distance = min(sl_df['metro_distance'])
            sl_df = sl_df[sl_df['metro_distance'] == min_distance]
            choice = sl_df.sample(n=1)
    
        
        #Change all 
        if lu == 'I':
            vs = self.g.vs.select(g8_number = int(choice['g8_number']))
        elif lu in ['R1', 'C', 'P']:
            vs = self.g.vs.select(g4_number = int(choice['g4_number']))
        else:
            vs = self.g.vs.select(g2_number = int(choice['g2_number']))
            
            
        #Develop the roads around the blocks first
        block_number = vs[0]['block_number']
        v_selection = self.g.vs.select(block_number = block_number)
        
        for v1 in v_selection:
            v_neighbours = v1.neighbors()
            
            for v2 in v_neighbours:
                if v2['type'] == 'road':
                    v2['landuse'] = 'RR'
                    v2['lot_number'] = 'RR'
                    v2['color'] = 'black'
        
        vs["landuse"] = lu
        self.latest_lot_number = self.latest_lot_number + 1
        vs["lot_number"] = self.latest_lot_number
        
        if lu == 'R2':
            vs["lot_number"] = self.latest_lot_number, self.latest_lot_number + 1
            self.latest_lot_number = self.latest_lot_number + 1
        
        vs["color"] = self.landuse_color_dic[lu]
        
    def plot_field(self):
            
        lot_polygons = []
        landuse = []
        colours = []
        
        for i in list(set(self.g.vs.select()['lot_number'])):

            lots = self.g.vs.select(lot_number = i)
            landuse = lots[0]['landuse']
            colours.append(self.landuse_color_dic[landuse])
            
            points = lots['coords']
            polygons = []
        
            for lot in points:
        
                polygon = Polygon(lot)
                polygons.append(polygon)
        
            merge = cascaded_union(polygons)
            merge_patch = PolygonPatch(merge)
            lot_polygons.append(merge_patch)
            
            
        pc = plt_c.PatchCollection(lot_polygons, edgecolor='black', facecolor = colours)
        
        fig, ax = plt.subplots(figsize=(30,50))
        ax.set_aspect("equal")
        
        ax.add_collection(pc)
        ax.set_title('Fakefield City')
        xrange = [-1, max(self.x_coords)]
        yrange = [min(self.y_coords), -1]
        ax.set_xlim(xrange)
        
        ax.set_ylim(yrange)
        
        ax.set_aspect(1)
        
        