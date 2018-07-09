#import Tkinter,tkFileDialog

#root = Tkinter.Tk()
#f = tkFileDialog.askopenfilenames(parent=root,title='Choose a file')
#print(root.tk.splitlist(f))

import os
import pandas as pd
import numpy as np
from tkinter.filedialog import askopenfilenames
from scipy.stats import mstats

def ir_txt2csv():
    f = askopenfilenames(title='Select Files')

    header_end = 7
    line_csv = ''
    for ind in f:
       with open(ind[0:-3]+'csv', 'w') as outfile, open(ind, 'r', encoding='utf-8') as infile:
            fname = os.path.basename(ind)
            bird = fname[0:4]

            # need to edit the below for other experiments 
            loc = fname.find('O')
            eye = fname[loc:loc+2]
            loc2 = fname.find('P')
            day = fname[loc2+1:loc2+2]
            print(day)
            outfile.write('bird,eye,day,num,axis,pupil,D\n')
            for num, line in enumerate(infile,1):
                  if num > header_end:
                      outfile.write(bird + ',' + eye + ',' + day  + ',' + line.replace('\t',','))    
    if not f:
        return True
    else:
        return False

def ir_csv2avgaxis():
    # empty dict
    data = []
    # winsorizing limits
    w_limits = [0.05, 0.05]

    f = askopenfilenames(title='Select Files')
    if f:
        # get parent directory
        uppath = lambda _path, n: os.sep.join(_path.split(os.sep)[:-n])
        parent_dir = uppath(f[0],2)
        for ind in f:
            with open(ind, 'r', encoding='utf-8') as infile:
                df = pd.read_csv(ind)
                re_h_v = df['D'].groupby(df['axis']).apply(mstats.winsorize, limits = w_limits)
                re = np.mean([np.mean(re_h_v[0]), np.mean(re_h_v[1])])
                data.append([df.iloc[0]['bird'], df.iloc[0]['eye'], df.iloc[0]['day'], re])
        out_df = pd.DataFrame(data, columns = ['bird','eye', 'day', 'D'])
        # need average the two eyes and replace day variable to get mean refraction of the two eyes
        avg_df = out_df['D'].groupby(out_df['bird']).mean().reset_index()
        avg_df['day'] = out_df['day']
        avg_df.to_csv(parent_dir + '/DAY_' + str(df.iloc[0]['day']) + '.csv', index = False)
        return False
    else:
        # send True if we're empty and done
        return True

def ir_avgaxis2eye_day_combine():
    condition = input('Enter condition: ')
    f = askopenfilenames(title = "Is this even used?")
    for ind in f:
        with open(ind, 'r', encoding = 'utf-8') as infile:
            df = pd.read_csv(ind)
            re_eyes = df['D'].groupby(df['bird']).mean().unstack

# convert IR txt file to a csv
done = False
while done == False:
    done = ir_txt2csv()

# convert raw csv to a csv with axis averaged
done = False
while done == False:
    done = ir_csv2avgaxis()
# average OD and OS and combine across days
#done