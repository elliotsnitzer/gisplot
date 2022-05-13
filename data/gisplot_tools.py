#!/usr/bin/env python

def update_progress(progress):
    #generates loading animation throughout plotting process
    title = 'Plotting Uploaded Data:'
    bar_length = 20
    block = int(20.0*progress)
    text = title+" [{0}] {1:.1f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    output_widget.clear_output(wait = True)
    with output_widget:
        print(text)


def load_icesheet_data(file_names):
    #update_progress function is used to show the loading animation while running plotting functions
    update_progress(0)
    #compiling file names into a list
    counter = 0
    global file_vars
    file_vars = []
    for file in file_names:
        file_vars.append(file)
        counter = counter+1
    
    global mapping_var
    if 'AIS' in file_vars[0]:
        mapping_var = 'mapping'
    elif 'GIS' in file_vars[0]:
        mapping_var = 'Polar_Stereographic'
        
    #converting netCDF files into Xarrays and storing them in a list
    counter = 0
    global model_vars
    model_vars = []  
    for f in file_vars:
        #decoded_times = False is needed to process these files
        #As far as I can tell it is only needed for AIS data
        model_vars.append(xr.open_dataset(f,engine='netcdf4',decode_times=False)) 
        counter = counter+1
    
    update_progress(0.1)  
    #collect titles for each of the subplots using the file names
    global ctvs, titles
    ctvs = []
    titles = []
    for f in file_names:
        fbase = os.path.basename(f)
        t = fbase.replace('.nc','')
        lt = t.split('_')
        #try catch to catch error in case the 2D variable name does not exist in the file
        try:
            #ctvs lilst stores the 2D variable names      ##might no longer be needed after adding references
            ctvs.append(lt[0])
        except:
            with output_widget:
                print('Data Variable missing from netCDF file')
            #"return" returns None, which simply exits the function
            return
        title = lt[2]+'_'+lt[3]
        titles.append(title)
    
    #function grabs and stores polar stereographic data from each Xarray 
    counter = 0
    global ctv_vars, ctv_proj_vars
    ctv_vars = []
    ctv_proj_vars = []
    try:
        for m in model_vars:
            ctv_vars.append(m[ctvs[counter]])
            ctv_proj_vars.append(m[mapping_var])
            counter = counter+1
    except:
        with output_widget:
            print('')
        
    update_progress(0.2) 
    
    #sets the stereographic projection based on the user's choice between a model based and standard projection
    global polar_stereographic
    if stereograph_choice.value == 'standard':
        if mapping_var == 'Polar_Stereographic':
            #Set Standard Polar Sterographic Projection definition
            polar_stereographic = ccrs.Stereographic(
                central_latitude=90.0,
                central_longitude=-45.0,
                false_easting=0.0,
                false_northing=0.0,
                true_scale_latitude=70.0,
                globe=ccrs.Globe('WGS84'))
        elif mapping_var == 'mapping':
            polar_stereographic = ccrs.Stereographic(
                central_latitude=-90.0,
                central_longitude=0.0,
                false_easting=0.0,
                false_northing=0.0,
                true_scale_latitude=-71.0,
                globe=ccrs.Globe('WGS84'))
    
    else:
        polar_stereographic = ccrs.Stereographic(
            central_latitude=ctv_proj_vars[0].latitude_of_projection_origin,
            central_longitude=ctv_proj_vars[0].straight_vertical_longitude_from_pole,
            false_easting=ctv_proj_vars[0].false_easting,
            false_northing=ctv_proj_vars[0].false_northing,
            true_scale_latitude=ctv_proj_vars[0].standard_parallel,
            globe=ccrs.Globe('WGS84') )
            
    update_progress(0.3)
    
    #calls the next function in the plotting process
    try:
        check = transform_projection()
        if check=='failed':
            upload_button.clear_output()
            return check
    except:
        with output_widget:
            print('ERROR: Projection Transformation Failed')
        return
    
    
def transform_projection(): 
    ######################
    # Transform projection
    #setting cartopy map values based on WGS84
    geodetic = ccrs.Geodetic(globe=ccrs.Globe('WGS84'))
    
    yv, xv = np.meshgrid(model_vars[0].y.data, model_vars[0].x.data)
    ll = geodetic.transform_points(src_crs=polar_stereographic, x=xv.flatten(), y=yv.flatten())
    global lons,lats
    lons = ll[:,0]
    lats = ll[:,1]
    
    update_progress(0.4)
    
    counter = 0
    global ctv_mean_vars
    ctv_mean_vars = []
    for l in ctv_vars:
        ctv_mean_vars.append(l.mean(dim='time').data)
        ctv_mean_vars[counter] = ctv_mean_vars[counter].transpose()
        ctv_mean_vars[counter] = ctv_mean_vars[counter].flatten()
        counter = counter+1
        
    update_progress(0.5)
    
    try:
        check = plot_icesheet_data()
        if check == 'failed':
            return check
    except:
        output_widget.clear_output()
        with output_widget:
            print('ERROR: Plotting Failed')
        time.sleep(3)
        return 'failed'
    
    #line not needed just used to profile the following function
#    cProfile.run('plot_icesheet_data()')



def plot_icesheet_data():    
    ####################
    #Plot Transformed Ice Sheet Data
    update_progress(0.6)
    
    #get dimensions for subplotting images
    num_frames = len(ctv_vars)
    if num_frames == 1:
        grid_cols = 1
        grid_rows = 1
        fig_dims = [10,10]
    else:
        grid_cols = 2
        grid_rows = math.ceil(math.sqrt(num_frames))
        fig_dims = [16,8*grid_rows]
    #creates a figure using the dimensions calculated
    plt.figure(figsize=(fig_dims[0],fig_dims[1]))
    
    update_progress(0.7)
    #using all previously gathered data to create the subplots for each dataset
    try:
        counter = 0
        ax_vars = []
        for l in ctv_mean_vars:
            frame = counter+1
            ax_vars.append(plt.subplot(grid_rows,grid_cols,frame, projection=polar_stereographic))
            if mapping_var == 'Polar_Stereographic':
                ax_vars[counter].set_extent([-65, -20, 57, 84]) #not needed for ais plots
#        elif mapping_var == 'mapping':
#            pass
            #ax_vars[counter].set_extent([-180,-160,183,30])#shortest = [-180,-160,180,-10]
            #if this doesn't end up speeding it up just delete
            #ax_vars[counter].set_extent([-65, -20, 57, 84]) #changing these values may speed up the plotting
            ax_vars[counter].coastlines(resolution='10m', zorder=7)
            ax_vars[counter].gridlines(zorder=8)
        #appropriately names the reference data subplot
            if counter==0:
                ax_vars[counter].set_title('Reference Data ('+titles[0]+')', fontsize=20)
            else:
                ax_vars[counter].set_title(titles[counter], fontsize=20)
#        if stereograph_choice.value == 'standard':
            plt.scatter(lons, lats, 1, c=ctv_mean_vars[counter], transform=ccrs.Geodetic(), zorder=0, cmap='viridis')
#        else:
#            plt.scatter(lons_vars[counter], lats_vars[counter], 1, c=ctv_mean_vars[counter], transform=ccrs.Geodetic(), zorder=0, cmap='viridis')
        #sets subplots axis name and units based on the data in the Xarray
            data_table = getattr(model_vars[counter],ctvs[counter])
            name = data_table.attrs['standard_name']
            units = data_table.attrs['units']
            c = plt.colorbar(fraction=0.046, pad=0.04)
            c.set_label('{0} ({1})'.format(name, units), size=16)
            counter = counter+1
        axes = plt.gca()
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    except:
        return 'failed'

    update_progress(0.8)
    #sets title for entire figure
    plt.suptitle((name+' ('+ctvs[0]+')'), fontsize=30) #subtitle name may need to be adjusted later
    plt.subplots_adjust(top=0.88)
    #saves figure for user to download later if they choose
    plt.savefig('GIS_Ice_Sheet_Model_Comparison.png', dpi='figure')
    
    update_progress(0.9)
    #print(ax_vars[0].get_extent())
    output_widget.clear_output(wait = True)
    with output_widget:
        plt.show() 
    with filename_output:
        print('Plotted Files:')
        for f in file_names:
            fbase = os.path.basename(f)
            print(fbase)