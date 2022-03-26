from flask import Flask, render_template, url_for, request, jsonify,send_file
from flask_cors import CORS, cross_origin
from werkzeug.utils import secure_filename
from scipy import stats
import sys
import os
import glob
import re
import csv
import json
import cdflib
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
from matplotlib import pyplot as plt

app=Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

def reduce_noise_by_stl_trend(rate, time):
    window_width = 60
    byteorder = rate.dtype.byteorder
    if byteorder=='=':
        byteorder=sys.byteorder
    print(byteorder)
    if byteorder == '<'or byteorder=='little':
        rate = (pd.Series(rate).rolling(window=window_width).mean().iloc[window_width-1:].values)
        time = (pd.Series(time).rolling(window=window_width).mean().iloc[window_width-1:].values)
    else:
        rate = (pd.Series(rate.byteswap().newbyteorder()).rolling(window=window_width).mean().iloc[window_width-1:].values)
        time = (pd.Series(time.byteswap().newbyteorder()).rolling(window=window_width).mean().iloc[window_width-1:].values)
#     fig = plt.figure(figsize=(18, 8))
    
    return rate,time
  

def find_peak(rate,time):
    x = rate
    peaks, _ = find_peaks(x)
    prominences, _, _ = peak_prominences(x, peaks)
    selected = prominences > 0.3 * (np.min(prominences) + np.max(prominences))
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
    window_width = 1
    flux = (pd.Series(df['flux']).rolling(window=window_width).mean().iloc[window_width-1:].values)
    time = (pd.Series(df['time']).rolling(window=window_width).mean().iloc[window_width-1:].values)
    x = flux
    peaks, _ = find_peaks(x)
    prominences, _, _ = peak_prominences(x, peaks)
    selected = prominences > 0.3 * (np.min(prominences) + np.max(prominences))
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
            try:
                if((x[j]-x[j+1])<0.00001 and (x[j+1]-x[j+2])<0.00001 and (x[j+2]-x[j+3])<0.00001 and x[j]<0.1*(np.max(x)-np.min(x))):


                    print("end coordinates : ",t,x[j])
                    end.append(x[j])
                    break
            except:
                    print("end coordinates : ",t,x[j-1])
                    end.append(x[j-1])
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
        try:
            x_f = np.delete(x, slice(si[0], ei[0]), 0)
        except:
            x_f=x
    flux_bc = np.mean(x_f, axis=0)
    
    return time[top],x[top], flux_bc

def get_bc(start, end, rate):
    x = rate
    for s in range(len(start)):
        si = np.where(x == start[s])[0]
        ei = np.where(x == end[s])[0]
        try:
            x = np.delete(x, slice(si[0], ei[0]), 0)
        except:
            continue
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
        flag=0
        os.system("xsmgenspec l1file=data/"+year+"/"+month+"/"+day+"/raw/ch2_xsm_"+year+month+day+"_v1_level1.fits specfile=data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+month+day+"_v1_flux.txt spectype=time-resolved tstart=0 tstop=0 tbinsize=1 hkfile=data/"+year+"/"+month+"/"+day+"/raw/ch2_xsm_"+year+month+day+"_v1_level1.hk safile=data/"+year+"/"+month+"/"+day+"/raw/ch2_xsm_"+year+month+day+"_v1_level1.sa gtifile=data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+month+day+"_v1_level2.gti")
        os.system("xsmcomputeflux  data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+""+month+""+day+"_v1_flux.txt data/"+year+"/"+month+"/"+day+"/calibrated/fluxc.txt 1.5498 12.398")
        os.system("rm -r data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+""+month+""+day+"_v1_flux.txt")
        flag = os.system("rm -r data/"+year+"/"+month+"/"+day+"/calibrated/ch2_xsm_"+year+""+month+""+day+"_v1_flux.arf") 
        if(flag!=0):
            os.system("touch data/"+year+"/"+month+"/"+day+"/calibrated/fluxc.txt")
            return 0
        else:
            return 1

def choose1(x,y,maxsize):
    binsize = int(len(x)/maxsize)
    xnew = np.array([])
    ynew = np.array([])
    count = 0
    while count < maxsize:
        count = int(count)
        xnew = np.append(xnew,x[count*binsize])
        ynew = np.append(ynew,y[count*binsize])
        count=count+1
    return xnew,ynew

def assign_status(df_rate,peak_time,start_time,end_time):
    df_rate['Status']='Normal'
    for i in range(0,len(df_rate)):    
        if(df_rate.iloc[i,1] in peak_time):
            df_rate.iloc[i,3] = 'Peak'
        elif(df_rate.iloc[i,1] in start_time):
            df_rate.iloc[i,3] = 'Start'
        elif(df_rate.iloc[i,1] in end_time):
            df_rate.iloc[i,3] = 'End'
        else:
            df_rate.iloc[i,3] = 'Normal'
    return df_rate



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
        file_path = f.filename
        if(file_path[-4:]!=".zip"):
            lcpath = f.filename
            is_flux=0
            os.system("touch tempflux.csv")
            flux_path = "tempflux.csv"
        else:
            year, month, day = store_data(file_path)
            is_flux = generate_flux(year, month, day)
            lcpath, flux_path = path(year, month, day)
        df = pd.DataFrame(columns=['file_name', 'start coordinate (x)', 'start coordinate (y)', 'peak coordinate (x)',
                                   'peak coordinate (y)', 'end coordinate (x)', 'end coordinate (y)', 'total burst time',
                                   'rise time', 'decay time', 'area under curve', 'background count Rate vs Time',
                                   'classfication by area', 'classification by duration'])
        flux_df = pd.DataFrame(columns=['flux_file_name', 'Peak Flux (x)', 'Peak Flux (y)', 'background count Flux vs Time',
                                        'Classification by Flux Peak', 'Classification by Flux Peak By Background Count'])
        # df1 = pd.read_table(flux_path, delimiter=' ', header=None)
        filetype = magic.from_file(lcpath)
        print("\n\n\n"+filetype+"\n\n\n")
        if 'ASCII' in filetype:
            df_x = pd.read_csv(lcpath, sep=" ", skipinitialspace=True)
            df_x.columns = ["time", "rate"]
            rate, time = reduce_noise_by_stl_trend(np.array(df_x["rate"], dtype=float), np.array(df_x["time"], dtype=float))
#             print(df_x)
            # df_x.to_excel("excel_try.xlsx",index=False)
        elif 'CSV' in filetype:
            df_x = pd.read_csv(lcpath)
            df_x.columns = ["time", "rate"]
            rate, time = reduce_noise_by_stl_trend(np.array(df_x["rate"], dtype=float), np.array(df_x["time"], dtype=float))
        elif 'Excel' in filetype:
            df_x = pd.read_excel(lcpath)
            df_x.columns = ["time", "rate"]
            #         plt.plot(df_x['time'], df_x['rate'])
#             print(df_x)
            rate = np.array(df_x['rate'])
            time = np.array(df_x['time'])
            # plt.plot(time,rate)
            rate, time = reduce_noise_by_stl_trend(rate, time)
        elif 'FITS' in filetype:
            image_file = fits.open(lcpath)
            file_data = image_file[1].data
            rate, time = reduce_noise_by_stl_trend(file_data["rate"], file_data["time"])
        elif 'FPT' in filetype or 'data' in filetype:
            cdf_file = cdflib.CDF(lcpath)
            arr = np.array((cdf_file.varget(variable='Sample Light Curve')[0]))
            rate = np.array([element[1] for element in arr])
            time = np.array([element[0] for element in arr])
            rate, time = reduce_noise_by_stl_trend(rate, time)
        #         df_flux.to_csv(flux_path+'.csv', index = None)
        # image_file = fits.open(lcpath)
        # file_data = image_file[1].data
        #rate, time = reduce_noise_by_stl_trend(file_data)
        #         rate_time_array = np.transpose(np.array([time,rate]))
        # time2,rate2 = time,rate
        if(len(time)>=1000):
            time,rate=choose1(time,rate,1000)
        df_rate = pd.DataFrame({ 'time':time, 'rate':rate}, index=None)
#         df_rate.to_csv(path+file_name+'.csv', index=None, header=False)
        top = find_peak(rate,time)
        start, start_index, start_time,peak,peak_time = get_start_point(top,rate,time)
        end, end_index,end_time = get_end_time(top,start,rate,time)
        
        print(start)
        print(end)
        print(rate)
        area = area_under_curve(rate, start_index, end_index)
        bc = get_bc(start, end, rate)
        area_class = classification_by_area(area)
        duration_class = classification_by_duration(start_time,end_time)
        

        

        df = append_to_dataframe(df,flux_path,start,start_time,end,end_time,peak,peak_time,area,bc,area_class,duration_class)
        df_rate['status'] ='Normal'
        #df_rate = assign_status(df_rate,peak_time,start_time,end_time)
        df_rate['status'] = df_rate['time'].apply(lambda x: 'Peak' if x in peak_time else('Start' if x in start_time else('End' if x in end_time else 'Normal')))
        #df_flux = assign_status(df_flux,flux_peak_time,[],[])
        if is_flux:
            df_flux = pd.read_csv(flux_path, delimiter = ' ',usecols = [2])
            df_flux.columns = ['flux']
            df_flux['time'] = df_flux.index
            df_flux = df_flux[['time', 'flux']]
            tm,rt = choose1(df_flux['time'],df_flux['flux'],1000)
            df_temp=pd.DataFrame()
            df_temp['time']=tm
            df_temp['flux'] =rt
            df_flux=df_temp
            flux_peak_time, flux_peak, flux_bc = flux_curve(df_flux)
            flux_class = classification_by_flux_peak(flux_peak)
            flux_class_bc = classification_by_flux_peak_by_bc(flux_peak,flux_bc)
            flux_df = flux_dataframe(flux_df,lcpath,flux_peak_time, flux_peak, flux_bc,flux_class,flux_class_bc)
            df_flux['status'] ='Normal'
            df_flux['status'] = df_flux['time'].apply(lambda x: 'Peak' if x in flux_peak_time else 'Normal')
        else:
            df_flux = pd.DataFrame(columns = ['flux_file_name','Peak Flux (x)','Peak Flux (y)','background count Flux vs Time','Classification by Flux Peak','Classification by Flux Peak By Background Count'])
        try:
            # lc_orig_df = pd.read_csv("CSV/lc.csv")
            # flux_orig_df = pd.read_csv("CSV/flux.csv")
            # all_lc_orig_df = pd.read_csv("CSV/all_lc.csv")
            # all_flux_orig_df = pd.read_csv("CSV/all_flux.csv")
            # pd.concat([lc_orig_df, df], ignore_index = True).to_csv("CSV/lc.csv", index=False)
            # pd.concat([flux_orig_df, flux_df], ignore_index = True).to_csv("CSV/flux.csv", index=False)
            # pd.concat([all_lc_orig_df, df_rate], ignore_index = True).to_csv("CSV/all_lc.csv", index=False)
            # pd.concat([all_flux_orig_df, df_flux], ignore_index = True).to_csv("CSV/all_flux.csv", index=False)
            df.to_csv(f'./CSV/lc.csv')
            df_rate.to_csv(f'./CSV/all_lc.csv')
            if is_flux:
                flux_df.to_csv(f'./CSV/flux.csv')
                df_flux.to_csv(f'./CSV/all_flux.csv')
        except:
            df.to_csv(f'./CSV/lc.csv')
            flux_df.to_csv(f'./CSV/flux.csv')
            df_rate.to_csv(f'./CSV/all_lc.csv')
            df_flux.to_csv(f'./CSV/all_flux.csv')

        return jsonify({'status': 'ok'})

@app.route('/api/data/lcfull', methods=['GET'])
@cross_origin()
def lcfulldata():
    try:
        lc_csv = pd.read_csv(r'CSV/all_lc.csv')
        lc_csv.to_json(r'JSON/all_lc.json')
        with open('JSON/all_lc.json', 'r') as file:
            lcJSON = file.read()
        return jsonify(lcJSON)
    except:
        return "No File Provided"

@app.route('/api/data/fluxfull', methods=['GET'])
@cross_origin()
def fluxfulldata():
    try:
        flux_csv = pd.read_csv(r'CSV/all_flux.csv')
        flux_csv.to_json(r'JSON/all_flux.json')
        with open('JSON/all_flux.json', 'r') as file:
            fluxJSON = file.read()
        return jsonify(fluxJSON)
    except:
        return "No File Provided"

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
