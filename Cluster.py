#!/usr/bin/env python
'''
Cornell University ASTR 1102 Lab: Star Clusters and Galactic Nebulae
Original Java Applet by Terry Herter et al.
Written by Michael Lam
'''

import matplotlib
matplotlib.use('TkAgg')

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
#from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from matplotlib.figure import Figure
#from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.ticker import *#MultipleLocator, FormatStrFormatter, LogLocator

import sys
import time

if sys.version_info[0] < 3:
    import Tkinter as Tk
else:
    import tkinter as Tk


DIR = "./"
SEPARATOR_COLOR="#CCCCCC"
SLEEP_TIME=1

##==================================================
## Physics
##==================================================

#Period-Luminosity Relationship
#Using Matsunaga et al. 2006
#M_V=mu_V*(log(P)-1.2)+eta_V (eq 10)
#Observational values: mu_V = -1.64, eta_V=-1.92
#Or for students: M_V = -1.64*log(P) + 0.048
#Period of V42: 25.735+/-0.015
def PLR(P,mu=-1.64,eta=-1.92):
    return mu*(np.log10(P)-1.2) + eta
    

M_sun=4.83 #http://nssdc.gsfc.nasa.gov/planetary/factsheet/sunfact.html
m_sun=-26.74
d_sun=4.822e-6
def get_age(V,distance=None):
    if distance==None:
        L=10**((M_sun-V)/2.5) #in solar units
    else:
        L=(distance/d_sun)**2 * 10**((m_sun-V)/2.5) #in solar units
    mass=L**(1/3.5)
    age=10*(mass**-2.5) #in billion years
    return age

##==================================================
## Data Processing
##==================================================

MAX_PROBABILITY=95



def loadM5():
    retval_V=list()
    retval_I=list()
    ID_list=list()
    temp=[0,0]
    FILE=open('M5_1.dat','r')
    lines=FILE.readlines()[54:3127]
    for i in range(len(lines)):
        if len(lines[i])<2 or lines[i][0:2]!='15':
            continue
        splitline=lines[i].replace('       ',' @ ').strip().split()

        if len(splitline)<10:#8:
            continue

        ID=splitline[6]

        if ID not in ID_list:
            ID_list.append(ID)
        if ID in ID_list and temp[0]!=0 and temp[1]!=0:
            retval_V.append(temp[0])
            retval_I.append(temp[1])
            temp=[0,0]
            continue #Just ignore re-writing points
        probability=splitline[7]
        #print ID,len(splitline)
        if probability != '@' and float(probability)>=MAX_PROBABILITY*0.01:
            if len(splitline)==10:
                V=splitline[9]
                I='@'
            elif len(splitline)==11:
                V='@'
                I=splitline[10]
#            print "foo",V,I
            if V != '@' or I != '@':
                if V != '@' and temp[0]==0:
                    temp[0]=float(V)
                elif I != '@' and temp[1]==0:
                    temp[1]=float(I) 
                ID_list.append(ID)
    FILE.close()
    FILE=open('M5_2.dat','r')
    lines=FILE.readlines()[39:10447]
    for i in range(len(lines)):
        splitline=lines[i].strip().split()
        probability=float(splitline[2])
        if probability >=MAX_PROBABILITY*0.01:
            retval_V.append(float(splitline[0]))
            retval_I.append(float(splitline[1]))
    
    FILE.close()

    FILE=open('M5_3.dat','r') #relatively useless
    lines=FILE.readlines()[47:178]
    for i in range(len(lines)):
        splitline=lines[i].strip().split()
        if len(splitline)<9:
            continue
        probability=float(splitline[6])
        if probability >=MAX_PROBABILITY*0.01:
            V=float(splitline[7])
            I=V-float(splitline[8])
            retval_V.append(V)
            retval_I.append(I)

    FILE.close()
    return np.array(retval_V),np.array(retval_I)



def loadPleiades(): #trial
    FILE=open('pleiadesdata.csv')
    lines=FILE.readlines()[5:]
    V=np.zeros(len(lines))
    B=np.zeros(len(lines))
    for i in range(len(lines)):
        stripline=lines[i].strip().split(',')
        V[i]=float(stripline[1])
        B[i]=float(stripline[2])
    FILE.close()

    testplot(B,V)
    


def loadM45():
    retval_V=list()
    retval_I=list() 
    FILE=open('M45.dat','r')
    lines=FILE.readlines()[11:3289]
    for i in range(len(lines)):
        splitline=lines[i].strip().split()
        probability=splitline[2]
        if probability!='~' and int(probability)>=MAX_PROBABILITY:
            V=splitline[-2]
            I=splitline[-1]
            if V!='~' and I!='~':
                retval_V.append(float(V))
                retval_I.append(float(I))
    FILE.close()
    return np.array(retval_V),np.array(retval_I)




def loadM67():
    retval_V=list()
    retval_I=list() 
    FILE=open('M67.dat','r')
    lines=FILE.readlines()[47:2456]
    for i in range(len(lines)):
        splitline=lines[i].replace('       ',' @ ').strip().split()
        probability=int(splitline[-1])
        if probability>=MAX_PROBABILITY:
            V=splitline[7]
            I=splitline[8]
            if V!='@' and I!='@':
                retval_V.append(float(V))
                retval_I.append(float(I))
    FILE.close()
    return np.array(retval_V),np.array(retval_I)



def loadMS(age):
    data = np.loadtxt('outputiso.dat')
    data = np.transpose(data)
    inds = np.where(np.logical_and(np.abs(data[0]-age)<0.001,data[4]<10))[0] #M<10 Msun
    retval_V = data[9][inds]
    retval_I = data[11][inds]
    return retval_V,retval_I


def loadZAMS():
    return loadMS(age=7.0)

    retval_V=list()
    retval_I=list()
    '''
    data=np.loadtxt('output.dat')
    data=np.transpose(data)
    inds=np.where(data[1]==6.0)[0]
    retval_V=data[10][inds]
    retval_I=data[12][inds]
    return retval_V,retval_I
    '''
    data=np.loadtxt('outputiso.dat')
    data=np.transpose(data)
    inds=np.where(np.logical_and(data[0]==7.0,data[4]<10))[0]
    retval_V=data[9][inds]
    retval_I=data[11][inds]
    return retval_V,retval_I



    FILE=open('isocsummz0.dat','r')
    
    lines=FILE.readlines()
    for line in lines:
        if line[0]=='#':
            continue
        splitline=line.strip().split()
        if splitline[-1]=='TO':
            V=float(splitline[6]) #this is absolute magnitude
            I=V-float(splitline[9])
            retval_V.append(V)
            retval_I.append(I)
                  
        
        
    FILE.close()
    return np.array(retval_V),np.array(retval_I)



def loadCepheid():
    FILE=open("cepheid.txt",'r')
    tV=[]
    V=[]
    eV=[]
    tI=[]
    I=[]
    eI=[]
    for line in FILE.readlines()[46:600]:
        splitline=line.strip().split()
        if splitline[1]=="V42": #Can use V84 as well
            if splitline[0]=="V":
                tV.append(float(splitline[2]))
                V.append(float(splitline[3]))
                eV.append(float(splitline[4]))
            elif splitline[0]=="I":
                tI.append(float(splitline[2]))
                I.append(float(splitline[3]))
                eI.append(float(splitline[4]))
    # In order to get D=7.5 kpc, using the value of M_V from the PLR above, we need to shift the magnitude so that the mean value is 12.11. This corresponds to a dV=0.79. Not sure why off.
    V=V-np.mean(V)+12.11
    return tV,V,eV,tI,I,eI #now including errors




def testplots():
    import matplotlib.pyplot as plt
    fig=plt.figure()
    def testplot(ax,V,I,title):
        ax.plot(V-I,V,'k.')
        ax.invert_yaxis()
        ax.set_xlabel('V-I')
        ax.set_ylabel('V')
        ax.set_title(title)
    VZ,IZ=loadZAMS()
    V,I=loadM5()
    ax=fig.add_subplot(221)
    
    testplot(ax,V,I,'M5')
    ax.plot(VZ-IZ,VZ+14,'r.') #17
    V,I=loadM45()
    ax=fig.add_subplot(222)
    testplot(ax,V,I,'M45')
    ax.plot(VZ-IZ,VZ+5,'r.') #7.5
    V,I=loadM67()
    ax=fig.add_subplot(223)
    testplot(ax,V,I,'M65')
    ax.plot(VZ-IZ,VZ+9,'r.') #13
    ax=fig.add_subplot(224)

    testplot(ax,VZ,IZ,'MS')            
    for age in np.arange(7.0,10.5,0.5):
        V,I=loadMS(age)
        ax.plot(V-I,V,'.')

    ax.set_ylim(14,-7)
    plt.tight_layout()
    plt.show()
    
#testplots()
#exit()




##==================================================
## GUI 
##==================================================

    

root = Tk.Tk()
#root.geometry('+1400+100')
root.geometry('+100+100')
root.wm_title("Star Clusters and Galactic Nebulae")



## ----------
## Build primary GUI containers
## ----------

mainframe = Tk.Frame(root)
mainframe.grid(row=0)

figframe = Tk.Frame(mainframe)#, bd = 6, bg='red')
fig = Figure(figsize=(8,5), dpi=75)
canvas = FigureCanvasTkAgg(fig, figframe)
canvas.get_tk_widget().grid(row=0)#,side=Tk.TOP)#,fill='x')
canvas.show()

canvas._tkcanvas.grid(row=1)#, fill=Tk.BOTH, expand=1)

figframe.grid(row=0,column=0)

## ----------
## Tkinter Variables
## ----------
var_mode = Tk.StringVar() # To determine which mode of operation
var_mode.set("ZAMS")
var_ZAMS_on = Tk.IntVar()
var_ZAMS_mag = Tk.StringVar()
var_distance = Tk.StringVar()
var_age = Tk.StringVar()
var_m_V = Tk.StringVar()
var_M_V = Tk.StringVar()
var_message = Tk.StringVar()
var_period = Tk.StringVar()

image_distance_modulus = Tk.PhotoImage(file=DIR+"DistanceModulus.gif")
image_period_magnitude = Tk.PhotoImage(file=DIR+"PeriodMagnitude.gif")

## ----------
## Primary Color-Magnitude Diagram
## ----------

ax_CMD = fig.add_subplot(111)

## pre-load MSs, ZAMS file. Takes a few seconds to start, but will not lag in the middle
ages=np.arange(7.0,10.5,0.5) #do not hardwire ages?
V_MS=dict()
I_MS=dict()
for age in ages: 
    V_MS[age],I_MS[age]=loadMS(age)

V_ZAMS,I_ZAMS=loadZAMS()



## primary plot updater
def update_CMD(label=None):
    if label==None: #Either by Mode menu or buttons
        label = var_mode.get()
    else:
        var_mode.set(label)
    ax_CMD.cla()
    if label=="ZAMS":
        #busy()
        turnoffs=[]
        for age in ages:#np.arange(7.0,10.5,0.5): #do not hardwire ages?
            V,I=V_MS[age],I_MS[age]
            VmI=V-I
            ax_CMD.plot(VmI,V,'.')
            if age!=7.0: #No turnoff here
                turnoff_ind=np.where(np.diff(VmI)>1e-3)[0][0] #rough cutoff
                turnoffs.append((VmI[turnoff_ind],V[turnoff_ind]))
                ax_CMD.plot(VmI[turnoff_ind],V[turnoff_ind],'o',color='0.5',markersize=7)

        ax_CMD.set_xlim(-0.9,5.5) #do not hard wire!
        ax_CMD.set_ylim(13,-6) #pre-invert
        ax_CMD.set_xlabel('V-I')
        ax_CMD.set_ylabel('V')
        #Draw on ages

        Vlim=ax_CMD.get_ylim()

        ax_CMD.text(-0.7,5,"Age (Gyr)")
        for VmI,V in turnoffs:
            age=get_age(V)
            ax_CMD.text(-0.7,V,"%0.2f"%age)
        #notbusy()
    elif label=="Cepheid":
        tV,V,eV,tI,I,eI=loadCepheid()
        tV-=np.min(tV)
        period = var_period.get()

        if period!="":
            try:
                tV=tV%float(period)
            except ValueError:
                busy("Error: Bad Input, ignoring period",sleep=SLEEP_TIME)
                period=""

        ax_CMD.errorbar(tV,V,yerr=eV,fmt='k.')

        if period=="":
            #ax_CMD.plot(tV,V,'k')
            ax_CMD.set_xlabel('Days')
        else:
            ax_CMD.set_xlabel('Days (Wrapped by Period)')
        dt=np.ptp(tV)
        ax_CMD.set_xlim(np.min(tV)-0.1*dt,np.max(dt)+0.1*dt)
        ax_CMD.set_ylabel('V')
        
    else:
        ZAMS_on = var_ZAMS_on.get()
        #exec("V,I=load%s()"%label)
        if label=="M5":
            V,I=loadM5()
        elif label=="M45":
            V,I=loadM45()
        elif label=="M67":
            V,I=loadM67()
        VmI=V-I
        ax_CMD.plot(VmI,V,'k.')
        ax_CMD.set_xlim(min(VmI)-0.5,max(VmI)+0.5)
        ax_CMD.set_ylim(max(V)+0.5,min(V)-0.5) #pre-invert
                
        if ZAMS_on!=0:
            VZ,IZ=V_ZAMS,I_ZAMS
            ZAMS_mag = var_ZAMS_mag.get()
            if ZAMS_mag == "":
                ZAMS_mag = 0
            else:
                try:
                    ZAMS_mag = float(ZAMS_mag)
                except ValueError:
                    busy("Error: Bad Input, setting ZAMS to 0",sleep=SLEEP_TIME)
                    ZAMS_mag=0
            ax_CMD.plot(VZ-IZ,VZ+ZAMS_mag,'r',linewidth=3)


        ax_CMD.set_xlabel('V-I')
        ax_CMD.set_ylabel('V')
    
    ax_CMD.set_title(label)

    canvas.draw()


## ----------
## Distance calculator
## ----------

def distance_calculator():
    m = var_m_V.get()
    M = var_M_V.get()
    if m!="" or M!="":
        try:
            m = float(m)
            M = float(M)
            dist = 10**((m-M+5)/5.0)
            if len("%0.2f"%dist)>11:
                var_distance.set("Overflow")
            else:
                var_distance.set("%0.2f" % dist)
        except ValueError:
            busy("Error: Bad Input",sleep=SLEEP_TIME)
            var_distance.set("")
            
    else:
        var_distance.set("")

def age_calculator():
    M = var_M_V.get()
    m = var_m_V.get()
    if m!="" and M!="":
        try:
            m=float(m)
            M=float(M)
            dist = 10**((m-M+5)/5.0)
            age = get_age(m,distance=dist)
            if len("%0.2f"%age)>10: #could be 11, but need some space
                var_age.set("Overflow")
            else:
                if age<1:
                    var_age.set("%0.3f" % age)
                else:
                    var_age.set("%0.2f" % age)
        except ValueError:
            busy("Error: Bad Input",sleep=SLEEP_TIME)
            var_age.set("")
    else:
        var_age.set("")

#NEW MODIFICATION:
toolbar_frame = Tk.Frame(mainframe)
toolbar_frame.grid(row=1,sticky=Tk.W)
toolbar = NavigationToolbar2TkAgg(canvas, toolbar_frame)
toolbar.update()
#toolbar.grid(row=1,sticky=Tk.W)

separator = Tk.Frame(mainframe,width=600,height=2,bg=SEPARATOR_COLOR,bd=1, relief=Tk.SUNKEN).grid(row=2,pady=2)

frame_buttons = Tk.Frame(mainframe)
frame_buttons.grid(row=3,sticky=Tk.W)

frame_mode = Tk.Frame(frame_buttons)
frame_mode.grid(row=0,column=0,sticky=Tk.W)

label_mode = Tk.Label(frame_mode,text="Mode Select:").grid(row=0,column=0,columnspan=2,sticky=Tk.N+Tk.W)

#Radio button construction
modes=["M5","M45","M67","ZAMS","Cepheid"]
for c in range(0,2):
    for r in range(1,4):
        ind = 3*c+r-1
        if ind==5:
            break
        mode=modes[ind]
        radiobutton = Tk.Radiobutton(frame_mode,text=mode,variable=var_mode,value=mode,indicatoron=0,command=update_CMD)
        radiobutton.grid(row=r,column=c,sticky=Tk.W)
#lambda mode=mode:update_CMD(mode)).grid(row=r,column=0,sticky=Tk.W)
'''
ROWVAR=1
radiobutton_ZAMS = Tk.Radiobutton(frame_mode,text="ZAMS",variable=var_mode,value="ZAMS",indicatoron=0,command=update_CMD).grid(row=ROWVAR,column=1,sticky=Tk.W)
ROWVAR+=1
radiobutton_Cepheid = Tk.Radiobutton(frame_mode,text="Cepheid",variable=var_mode,value="Cepheid",indicatoron=0,command=update_CMD).grid(row=ROWVAR,column=1,sticky=Tk.W)
'''
separator = Tk.Frame(frame_buttons,width=2,height=100, bg=SEPARATOR_COLOR,bd=1, relief=Tk.SUNKEN).grid(row=0,column=1,padx=2)


frame_options = Tk.Frame(frame_buttons)
frame_options.grid(row=0,column=2)

checkbutton_ZAMS = Tk.Checkbutton(frame_options,text="Overplot ZAMS",variable=var_ZAMS_on,command=update_CMD)
checkbutton_ZAMS.grid(row=0,column=0)

#checkbutton_age = Tk.Checkbutton(frame_options,text="Overplot Ages",variable=var_ages_on,command=update_CMD)
#checkbutton_age.grid(row=0,column=1)

label_ZAMS = Tk.Label(frame_options,text="Shift magnitude:")
label_ZAMS.grid(row=1,column=0)
entry_ZAMS = Tk.Entry(frame_options,width=7,textvariable=var_ZAMS_mag)
entry_ZAMS.grid(row=1,column=1)

label_period = Tk.Label(frame_options,text="Cepheid Period Estimate")
label_period.grid(row=2,column=0)
entry_period = Tk.Entry(frame_options,width=7,textvariable=var_period)
entry_period.grid(row=2,column=1)

button_redraw = Tk.Button(frame_options,text="Redraw",command=update_CMD)
button_redraw.grid(row=3,column=1)
#radiobutton_Cepheid = Tk.Radiobutton(frame_mode,text="Cepheid",variable=radiobutton_var,value=5,indicatoron=0).grid(row=ROWVAR,column=1,sticky=Tk.W)


separator = Tk.Frame(frame_buttons,width=2,height=100, bg=SEPARATOR_COLOR,bd=1, relief=Tk.SUNKEN).grid(row=0,column=3,padx=2)

frame_calculator=Tk.Frame(frame_buttons)
frame_calculator.grid(row=0,column=4)

label_calculator = Tk.Label(frame_calculator,text="Distance/Lifetime Calculator:").grid(row=0,column=0,columnspan=3,sticky=Tk.N+Tk.W)

label_m=Tk.Text(frame_calculator,width=2,height=1.6,borderwidth=0,background=frame_calculator.cget("background"))
label_m.tag_configure("subscript",offset=-4)
label_m.insert("insert","m","","V","subscript")
label_m.configure(state="disabled")
label_m.grid(row=1,column=0,sticky=Tk.W)
label_M=Tk.Text(frame_calculator,width=2,height=1.6,borderwidth=0,background=frame_calculator.cget("background"))
label_M.tag_configure("subscript",offset=-4)
label_M.insert("insert","M","","V","subscript")
label_M.configure(state="disabled")
label_M.grid(row=1,column=1,sticky=Tk.W)

entry_m=Tk.Entry(frame_calculator,width=5,textvariable=var_m_V)
entry_m.grid(row=2,column=0,sticky=Tk.W)
entry_M=Tk.Entry(frame_calculator,width=5,textvariable=var_M_V)
entry_M.grid(row=2,column=1,sticky=Tk.W)


label_distance = Tk.Label(frame_calculator,text="Distance (pc):")
label_distance.grid(row=1,column=2,sticky=Tk.E)
label_distance_value = Tk.Label(frame_calculator,textvariable=var_distance,width=10)
label_distance_value.grid(row=2,column=2,sticky=Tk.E)
button_distance = Tk.Button(frame_calculator,text="Calculate Distance",command=distance_calculator)
button_distance.grid(row=3,column=0,columnspan=2,sticky=Tk.W)


label_age = Tk.Label(frame_calculator,text="Lifetime (Gyr):")
label_age.grid(row=1,column=3,sticky=Tk.W)
label_age_value = Tk.Label(frame_calculator,textvariable=var_age,width=10)
label_age_value.grid(row=2,column=3,sticky=Tk.W)
button_age = Tk.Button(frame_calculator,text="Calculate Age",command=age_calculator)
button_age.grid(row=3,column=2,columnspan=2,sticky=Tk.W)

separator = Tk.Frame(mainframe,width=600,height=2,bg=SEPARATOR_COLOR,bd=1, relief=Tk.SUNKEN).grid(row=4,pady=2)

frame_message = Tk.Frame(mainframe)
frame_message.grid(row=5,column=0,stick=Tk.W)

label_message = Tk.Label(frame_message,textvariable=var_message)
label_message.grid(row=0)




## ----------
## Buttons/Menus
## ----------

def busy(msg="Working...",sleep=0):
    var_message.set(msg)
    root.config(cursor="watch")
    root.update()#_idletasks() #need to work through queued items
    if sleep!=0:
        time.sleep(sleep)
        notbusy()

def notbusy():
    var_message.set("")
    root.config(cursor="")


def popup_about():
    title="About"
    text=["Cornell University Department of Astronomy",
          "ASTR 1102 Lab: Star Clusters and Galactic Nebulae",
          "Original program by Martha P. Haynes,",
          "Jim B. Marshall, and Marc Parmet, 1987",
          "Python code by Michael Lam 2013",
          "",
          "M5, M67, Cepheid cluster data obtained through VizieR",
          "M5:  Sandquist+ 2004, Layden+ 2005, Rees 1993",
          "M45: SIMBAD4 rel 1.207  -  2013.07.31CEST05:25:05",
          "and references within",
          "M67: Yadav+ 2008",
          "Cepheid: Rabidoux+ 2010",
          "ZAMS: http://stev.oapd.inaf.it/YZVAR/, Bertelli+ 2009"]
    d = window_popup(root,title,text,WIDTH=50)
    root.wait_window(d.top)

def popup_commands():
    title="Commands"
    text=["Click on a mode to select a graph",
          "",
          "Overplot ZAMS: Plots a red ZAMS on a CMD",
          "Shift magnitude: Shifts the magnitude of the ZAMS",
          "Cepheid Period Estimate: Folds Cepheid data by this amount",
          "",
          "Calculate Distance, Lifetime: requires both m and M",
          "Remember that lifetime != current age for most stars"]
    d = window_popup(root,title,text,WIDTH=50)
    root.wait_window(d.top)
    
def popup_equations():
    title="Useful Equations"
    text=["Distance modulus:",
          "image_distance_modulus",
          "Cepheid II Period-Magnitude Relation (Matsunaga et al. 2006):",
          "image_period_magnitude"]

    d = window_popup(root,title,text,WIDTH=50)#,photo)
    root.wait_window(d.top)


class window_popup:
    def __init__(self,parent,title,txt,WIDTH=40):
        top = self.top = Tk.Toplevel(parent)
        top.title(title)
        top.geometry('+150+250')
        top.bind("<Return>",lambda event:self.ok())
        for i in range(len(txt)):
            if txt[i][:5]=="image":
                photo = eval(txt[i])
                label=Tk.Label(top,image=photo)
                label.image = photo # keep a reference!
                label.pack()
            else:
                Tk.Label(top,anchor=Tk.W,width=WIDTH,text=txt[i]).pack()
        b = Tk.Button(top,text="OK",command=self.ok)
        b.pack()
        b.focus_set()
    def ok(self):
        self.top.destroy()



def destroy(event):
    sys.exit()


def superdo(event):
    update_CMD()
    distance_calculator()
    age_calculator()
    

## Bindings
#root.bind("<Return>",superdo)
root.bind("<Escape>", destroy)
root.bind("<Control-q>", destroy)
root.bind("<F1>",lambda event: popup_about())
root.bind("<F2>",lambda event: popup_commands())
root.bind("<F3>",lambda event: popup_equations())
root.bind("<F10>",destroy)
root.bind("<Control-Key-1>",lambda event: update_CMD("M5"))
root.bind("<Control-Key-2>",lambda event: update_CMD("M45"))
root.bind("<Control-Key-3>",lambda event: update_CMD("M67"))
root.bind("<Control-Key-4>",lambda event: update_CMD("ZAMS"))
root.bind("<Control-Key-5>",lambda event: update_CMD("Cepheid"))




menubar = Tk.Menu(root)

filemenu = Tk.Menu(menubar, tearoff=0)
filemenu.add_command(label="Exit",accelerator="Esc", command=root.quit)
menubar.add_cascade(label="File", menu=filemenu)

modemenu = Tk.Menu(menubar, tearoff=0)
modemenu.add_command(label="M5",accelerator="Ctrl+1", command=lambda: update_CMD("M5"))
modemenu.add_command(label="M45",accelerator="Ctrl+2", command=lambda: update_CMD("M45"))
modemenu.add_command(label="M67",accelerator="Ctrl+3", command=lambda: update_CMD("M67"))
modemenu.add_separator()
modemenu.add_command(label="ZAMS",accelerator="Ctrl+4", command=lambda: update_CMD("ZAMS"))
modemenu.add_separator()
modemenu.add_command(label="Cepheid",accelerator="Ctrl+5", command=lambda: update_CMD("Cepheid"))
menubar.add_cascade(label="Mode", menu=modemenu)

helpmenu = Tk.Menu(menubar, tearoff=0)
helpmenu.add_command(label="About",accelerator="F1", command=popup_about)
helpmenu.add_command(label="Commands",accelerator="F2", command=popup_commands)
helpmenu.add_command(label="Useful Equations",accelerator="F3", command=popup_equations)
menubar.add_cascade(label="Help", menu=helpmenu)

# display the menu
root.config(menu=menubar)

update_CMD()

#root.configure(cursor=("@/usr/X11R6/include/X11/bitmaps/star","/usr/X11R6/include/X11/bitmaps/starMask", "white", "black"))

root.mainloop()
#Tk.mainloop()
