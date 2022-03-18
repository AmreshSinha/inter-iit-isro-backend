from flask import Flask, render_template, url_for, request, jsonify
from flask_cors import CORS, cross_origin
from werkzeug.utils import secure_filename

import sys
import os
import glob
import re
import csv
import json

from astropy.utils.data import get_pkg_data_filename
import magic
from astropy.table import Table
from astropy.io import fits
import numpy as np
from scipy.signal import find_peaks, peak_prominences
import pandas as pd
from statsmodels.tsa.seasonal import STL
from scipy.integrate import simps
from numpy import trapz

app=Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

def reduce_noise_by_stl_trend(file_data):
    window_width = 60
    
    rate = (pd.Series(file_data['Rate'].byteswap().newbyteorder()).rolling(window=window_width).mean().iloc[window_width-1:].values)
    time = (pd.Series(file_data['Time'].byteswap().newbyteorder()).rolling(window=window_width).mean().iloc[window_width-1:].values)
    
    
    return rate,time

def find_peak(rate,time):
    x = rate
    peaks, _ = find_peaks(x)
    prominences, _, _ = peak_prominences(x, peaks)
    selected = prominences > 0.5 * (np.min(prominences) + np.max(prominences))
    top = peaks[selected]
    topp = []
    for i in range(0,len(top)):
        if x[top][i]>600:
             topp.append(top[i])
    
    return topp

def get_start_point(top,rate, time):
    x = rate   
    i = top
    start = []
    start_index = []
    start_time = []
    peak = []
    peak_time = []
    for i in top:
        print("peak coordinates : ",time[i],x[i])
        peak.append(x[i])
        peak_time.append(time[i])
        while i>0:
            t = time[i]
            
            if((x[i]-x[i-1])<0.00001 and (x[i-1]-x[i-2])<0.00001 and (x[i-2]-x[i-3])<0.00001 and x[i]<0.5*(np.max(x)-np.min(x))):
                print("start coordinates : ",t,x[i])
                start.append(x[i])
                start_time.append(time[i])
                start_index.append(i)
                break

            if( x[i-1]>=x[i] and x[i-2]>=x[i-1] and x[i-3]>=x[i-2] and x[i-3]>=1.005*x[i]):
                
                print("start coordinates : ",t,x[i])
                start.append(x[i])
                start_index.append(i)
                start_time.append(time[i])
                break
                
            i = i-1
        
    return start, start_index, start_time, peak , peak_time  

def get_end_time(top,start,rate,time):
    x = rate   
    i = top
    end = []
    end_index = []
    end_time = []
    for i in top:
        print("peak coordinates : ",time[i],x[i])
        while i:
            t = time[i]
            
            if((x[i]-x[i+1])<0.00001 and (x[i+1]-x[i+2])<0.00001 and (x[i+2]-x[i+3])<0.00001 and x[i]<0.5*(np.max(x)-np.min(x))):
                print("end coordinates : ",t,x[i])
                end.append(x[i])
                end_index.append(i)
                end_time.append(time[i])
                break

            if( x[i+1]>=x[i] and x[i+2]>=x[i+1] and x[i+3]>=x[i+2] and x[i+3]>=1.005*x[i]):
                print("end coordinates : ",t,x[i])
                end.append(x[i])
                end_index.append(i)
                end_time.append(time[i])
                break
                
            i = i+1
            if(i>=len(time)):
                break
        
    return end, end_index, end_time



def area_under_curve(rate, start_index, end_index):
    area = []
    for i in range(len(start_index)):
        y = rate[start_index[i]:end_index[i]+1]
        area.append(trapz(y, dx=1))
    return area


def flux_curve(df):
    df.fillna(0, inplace=True)
    window_width = 80
    flux = (pd.Series(df[2]).rolling(window=window_width).mean().iloc[window_width-1:].values)
    time = (pd.Series(df.index).rolling(window=window_width).mean().iloc[window_width-1:].values)
    x = flux
    peaks, _ = find_peaks(x)
    prominences, _, _ = peak_prominences(x, peaks)
    selected = prominences > 0.5 * (np.min(prominences) + np.max(prominences))
    top = peaks[selected]
    print(x[top])


    x = flux  
    i = top
    start = []
    start_index = []
    for i in top:
        print("peak coordinates : ",x[i])
        while i>0: 
            if((x[i]-x[i-1])<0.00001 and (x[i-1]-x[i-2])<0.00001 and (x[i-2]-x[i-3])<0.00001 and x[i]<0.1*(np.max(x)-np.min(x))):
                print("start coordinates : ",time[i],x[i])
                start.append(x[i])
                break

            if( x[i-1]>=x[i] and x[i-2]>=x[i-1] and x[i-3]>=x[i-2] and x[i-3]>=1.005*x[i]):

                print("start coordinates : ",time[i], x[i])
                start.append(x[i])
                start_index.append(i)
                break

            i = i-1

    x = flux   
    j = top
    end = []
    end_index = []
    for j in top:
        print("peak coordinates : ",time[j],x[j])
        while j:
            t = time[j]

            if((x[j]-x[j+1])<0.00001 and (x[j+1]-x[j+2])<0.00001 and (x[j+2]-x[j+3])<0.00001 and x[j]<0.1*(np.max(x)-np.min(x))):


                print("end coordinates : ",t,x[j])
                end.append(x[j])
                break

            if( x[j+1]>=x[j] and x[j+2]>=x[j+1] and x[j+3]>=x[j+2] and x[j+3]>=1.005*x[j]):
                print("end coordinates : ",t,x[j])
                end.append(x[j])
                end_index.append(j)
                break

            j = j+1
    for s in range(len(start)):
        si = np.where(x == start[s])[0]
        ei = np.where(x == end[s])[0]
        x_f = np.delete(x, slice(si[0], ei[0]), 0)
    flux_bc = np.mean(x_f, axis=0)
    
    return time[top],x[top], flux_bc

def get_bc(start, end, rate):
    x = rate
    for s in range(len(start)):
        si = np.where(x == start[s])[0]
        ei = np.where(x == end[s])[0]
        x = np.delete(x, slice(si[0], ei[0]), 0)
    bc = np.mean(x, axis=0)
    return bc

def classification_by_area(area):
    area_class = []
    for i in range(0,len(area)):
        if(area[i]>=1e6):
            area_class.append("BRIGHT")
        elif(area[i]<1e6 and area[i]>=1e5):
            area_class.append("NORMAL")
        else:
            area_class.append("FAINT")
    return area_class

def classification_by_duration(start_time,end_time):
    duration_class = []
    for i in range(0,len(start_time)):
        duration = end_time[i]-start_time[i]
        if(duration<=3600):
            duration_class.append("SHORT DURATION OR IMPULSIVE EVENT")
        else:
            duration_class.append("LONG DURATION OR GRADUAL EVENT")
    return duration_class


def append_to_dataframe(df,name,start,start_time,end,end_time,peak,peak_time,area,bc,area_class,duration_class):
    burst_time = []
    rise_time = []
    decay_time = []
    for i in range(0,len(peak)):
        burst_time.append(end_time[i]-start_time[i])
        rise_time.append(peak_time[i]-start_time[i])
        decay_time.append(end_time[i]-peak_time[i])
    dict = {'file_name':name,'start coordinate (x)':start_time, 'start coordinate (y)':start, 'peak coordinate (x)':peak_time, 'peak coordinate (y)':peak, 'end coordinate (x)':end_time, 'end coordinate (y)':end, 'total burst time':burst_time, 'rise time':rise_time, 'decay time':decay_time, 'area under curve':area,'background count Rate vs Time':bc, 'classfication by area':area_class, 'classification by duration':duration_class}
    df2 = pd.DataFrame(dict)

  
    df3 = pd.concat([df, df2], ignore_index = True)
    df3.reset_index()

    return df3

def classification_by_flux_peak(flux_peak):
    flux_class = []
    for i in range(0,len(flux_peak)):
        if(flux_peak[i]<1e-7):
            flux_class.append('A')
        elif(flux_peak[i]>1e-7 and flux_peak[i]<1e-6):
            flux_class.append('B')
        elif(flux_peak[i]>1e-6 and flux_peak[i]<1e-5):
            flux_class.append('C')
        elif(flux_peak[i]>1e-5 and flux_peak[i]<1e-4):
            flux_class.append('M')
        else:
            flux_class.append('X')
            
    return flux_class

def classification_by_flux_peak_by_bc(flux_peak,flux_bc):
    flux_class_bc = []
    for i in range(0,len(flux_peak)):
        if((flux_peak[i]/flux_bc)<10):
            flux_class_bc.append('Type 1')
        elif((flux_peak[i]/flux_bc)>10 and (flux_peak[i]/flux_bc)<100):
            flux_class_bc.append('Type 2')
        else:
            flux_class_bc.append('Type 3')
            
    return flux_class_bc

def append_to_dataframe(df,name,start,start_time,end,end_time,peak,peak_time,area,bc,area_class,duration_class):
    burst_time = []
    rise_time = []
    decay_time = []
    for i in range(0,len(peak)):
        burst_time.append(end_time[i]-start_time[i])
        rise_time.append(peak_time[i]-start_time[i])
        decay_time.append(end_time[i]-peak_time[i])
    dict = {'file_name':name,'start coordinate (x)':start_time, 'start coordinate (y)':start, 'peak coordinate (x)':peak_time, 'peak coordinate (y)':peak, 'end coordinate (x)':end_time, 'end coordinate (y)':end, 'total burst time':burst_time, 'rise time':rise_time, 'decay time':decay_time, 'area under curve':area,'background count Rate vs Time':bc, 'classfication by area':area_class, 'classification by duration':duration_class}
    df2 = pd.DataFrame(dict)

  
    df3 = pd.concat([df, df2], ignore_index = True)
    df3.reset_index()

    return df3

def flux_dataframe(df1,flux_file,flux_peak_time, flux_peak,flux_bc,flux_class,flux_class_bc):
    dict = {'flux_file_name':flux_file,'Peak Flux (x)':flux_peak_time,'Peak Flux (y)':flux_peak,'background count Flux vs Time':flux_bc,'Classification by Flux Peak':flux_class,'Classification by Flux Peak By Background Count':flux_class_bc}
    df2 = pd.DataFrame(dict)

  
    df3 = pd.concat([df1, df2], ignore_index = True)
    df3.reset_index()

    return df3

def store_data(zipname):
    name = zipname.split("_")[-2]
    year = name[:4]
    month = name[4:6]
    day = name[6:8]

    os.system("unzip "+zipname+" -d temp")#unzfileip
    os.system("rsync -av temp/xsm/data/ data")#merge
    os.system("rm -r temp")#temp removal
    os.system("rm -r "+zipname)#zip
    return year,month,day

#data/year/month/day/calibrated/ch2_xsm_yeardaymonth_v1_level2.lc
#                               fluxc.txt


def path(year,month,day):
    lcpath = "data/"+year+"/"+month+"/"+day+"/calibrated/"+"ch2_xsm_"+year+month+day+"_v1_level2.lc"
    flux_path = "data/"+year+"/"+month+"/"+day+"/calibrated/fluxc.txt"
    return lcpath,flux_path


def generate_flux(year,month,day):        
    os.system("xsmgenspec l1file=data/"+year+"/"+month+"/"+day+"/raw/ch2_xsm_"+year+month+day+"_v1_level1.fits specfile=data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+month+day+"_v1_flux.txt spectype=time-resolved tstart=0 tstop=0 tbinsize=1 hkfile=data/"+year+"/"+month+"/"+day+"/raw/ch2_xsm_"+year+month+day+"_v1_level1.hk safile=data/"+year+"/"+month+"/"+day+"/raw/ch2_xsm_"+year+month+day+"_v1_level1.sa gtifile=data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+month+day+"_v1_level2.gti")
    os.system("xsmcomputeflux  data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+""+month+""+day+"_v1_flux.txt data/"+year+"/"+month+"/"+day+"/calibrated/fluxc.txt 1.5498 12.398")
    os.system("rm -r data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+""+month+""+day+"_v1_flux.txt")
    os.system("rm -r data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+""+month+""+day+"_v1_flux.arf") 

def make_json(csvFilePath, jsonFilePath):
     
    # create a dictionary
    data = {}
     
    # Open a csv reader called DictReader
    with open(csvFilePath, encoding='utf-8') as csvf:
        csvReader = csv.DictReader(csvf)
         
        # Convert each row into a dictionary
        # and add it to data
        for rows in csvReader:
             
            # Assuming a column named 'No' to
            # be the primary key
            key = rows['No']
            data[key] = rows
 
    # Open a json writer, and use the json.dumps()
    # function to dump data
    with open(jsonFilePath, 'w', encoding='utf-8') as jsonf:
        jsonf.write(json.dumps(data, indent=4))

@app.route('/', methods=['GET'])
@cross_origin()
def home():
    return "<h1>Hi from Backend!</h1>"

@app.route('/api/upload', methods=['POST'])
@cross_origin()
def upload():
        # Get the file from post request
        f = request.files['imgfile']
        f.save(secure_filename(f.filename))
        file_path=f.filename

        year,month,day=store_data(file_path)
        generate_flux(year,month,day)
        lcpath,flux_path=path(year,month,day)


        df = pd.DataFrame(columns = ['file_name','start coordinate (x)', 'start coordinate (y)', 'peak coordinate (x)', 'peak coordinate (y)', 'end coordinate (x)', 'end coordinate (y)', 'total burst time', 'rise time', 'decay time', 'area under curve','background count Rate vs Time', 'classfication by area', 'classification by duration'])
        flux_df = pd.DataFrame(columns = ['flux_file_name','Peak Flux (x)','Peak Flux (y)','background count Flux vs Time','Classification by Flux Peak','Classification by Flux Peak By Background Count'])
        df1 = pd.read_table(flux_path, delimiter=' ', header=None)
        filetype = magic.from_file(lcpath)
        if 'ASCII' in filetype:
            table = Table.read(lcpath, format='ascii')
            table.write(lcpath, format='fits')
        elif 'XLS' in filetype:
            file_xls = pd.read_excel(lcpath)
            file_xls.to_csv(lcpath)
            table2 = Table.read(lcpath+".csv", format='pandas.csv')
            table2.write(lcpath, format='fits')
        df_flux = pd.read_csv(flux_path, delimiter = ' ')
#         df_flux.to_csv(flux_path+'.csv', index = None)
        image_file = fits.open(lcpath)
        file_data = image_file[1].data
        rate,time = reduce_noise_by_stl_trend(file_data)
        rate_time_array = np.transpose(np.array([time,rate]))
        df_rate = pd.DataFrame(rate_time_array)
#         df_rate.to_csv(path+file_name+'.csv', index=None, header=False)
        top = find_peak(rate,time)
        start, start_index, start_time,peak,peak_time = get_start_point(top,rate,time)
        end, end_index,end_time = get_end_time(top,start,rate,time)
        flux_peak_time, flux_peak, flux_bc = flux_curve(df1)
        area = area_under_curve(rate, start_index, end_index)
        bc = get_bc(start, end, rate)
        area_class = classification_by_area(area)
        duration_class = classification_by_duration(start_time,end_time)
        flux_class = classification_by_flux_peak(flux_peak)
        flux_class_bc = classification_by_flux_peak_by_bc(flux_peak,flux_bc)
        df = append_to_dataframe(df,flux_path,start,start_time,end,end_time,peak,peak_time,area,bc,area_class,duration_class)
        flux_df = flux_dataframe(flux_df,lcpath,flux_peak_time, flux_peak, flux_bc,flux_class,flux_class_bc)


        try:
            lc_orig_df = pd.read_csv("CSV/lc.csv")
            flux_orig_df = pd.read_csv("CSV/flux.csv")
            all_lc_orig_df = pd.read_csv("CSV/all_lc.csv")
            all_flux_orig_df = pd.read_csv("CSV/all_flux.csv")
            pd.concat([lc_orig_df, df], ignore_index = True).to_csv("CSV/lc.csv", index=False)
            pd.concat([flux_orig_df, flux_df], ignore_index = True).to_csv("CSV/flux.csv", index=False)
            pd.concat([all_lc_orig_df, df_rate], ignore_index = True).to_csv("CSV/all_lc.csv", index=False)
            pd.concat([all_flux_orig_df, df_flux], ignore_index = True).to_csv("CSV/all_flux.csv", index=False)
            df.to_csv(f'./CSV/lc.csv')
            flux_df.to_csv(f'./CSV/flux.csv')
            df_rate.to_csv(f'./CSV/all_lc.csv')
            df_flux.to_csv(f'./CSV/all_flux.csv')
        except:
            df.to_csv(f'./CSV/lc.csv')
            flux_df.to_csv(f'./CSV/flux.csv')
            df_rate.to_csv(f'./CSV/all_lc.csv')
            df_flux.to_csv(f'./CSV/all_flux.csv')

        return jsonify({'status': 'ok'})

@app.route('/api/data/lc', methods=['GET'])
@cross_origin()
def lcData():
    try:
        lc_csv = pd.read_csv(r'CSV/lc.csv')
        lc_csv.columns = lc_csv.columns.str.replace(' ','_')
        lc_csv.columns = lc_csv.columns.str.replace('(','_')
        lc_csv.columns = lc_csv.columns.str.replace(')','_')
        lc_csv.to_json(r'JSON/lc.json')
        with open('JSON/lc.json', 'r') as file:
            lcJSON = file.read()
        return jsonify(lcJSON)
    except:
        return "No File Provided"

@app.route('/api/data/flux', methods=['GET'])
@cross_origin()
def fluxData():
    try:
        flux_csv = pd.read_csv(r'CSV/flux.csv')
        flux_csv.columns = flux_csv.columns.str.replace(' ','_')
        flux_csv.columns = flux_csv.columns.str.replace('(','_')
        flux_csv.columns = flux_csv.columns.str.replace(')','_')
        flux_csv.to_json(r'JSON/flux.json')
        with open('JSON/flux.json', 'r') as file:
            fluxJSON = file.read()
        return jsonify(fluxJSON)
    except:
        return "No File Provided"

if __name__ == "__main__":
    app.run(debug=False, host='0.0.0.0', port=8080)

# debug was initially True
