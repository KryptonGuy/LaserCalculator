from tkinter import *
import math
import numpy as np
import matplotlib.pyplot as plt
from PIL import ImageTk, Image

global h, c
h =  6.632 * 1e-34
c = 3 * 10e8
pi = math.pi
# system wavelength
wl = 1.03e-6
# default system M^2 value
M2 = 1.2
# default propagation distance between the telescope and the coupling lens!
pd = 4000

def wav_feq(wav):
    feq = c/float(wav)
    return feq

def e2_e1(feq):
    return h * feq

def mode_num(wav, l, n = 1):
    return (2 * n * l)/wav

def mode_sep_feq(l , n=1):
    return c / (2 * n * l)

def width_of_curve(delv, feq):
    return (c * delv)/ (feq * feq)

def relative_band_width(delv, feq):
    return delv / feq

def rayleighrang(d, wav):
    return (d * d) / wav

def semi_angle_diffraction(d, wav):
    return (1.22 * wav)/ d

def intensity(power , diameter):
    return (4 * power) / (math.pi * diameter * diameter)

def reflectivity(n1, n2 = 1):
    temp = (n1- n2)/(n1 + n2)
    return temp * temp

def beam_waist(focal, dist, wav):
    return (2 * focal * wav)/(math.pi * dist)

def interaction_time(radius, velocity):
    return (2 * radius) / velocity    

def beamdiameter():
    pass
 
def rz(s, P, v):
    return (12.528 * math.pow(s, 0.542))/(math.pow(P, 0.528) * math.pow(v, 0.322))

def ra(s, P, v):
    return (2.018 * math.pow(s, 0.670))/(math.pow(P, 0.451) * math.pow(v, 0.330))

def bpp(M2):
    '''
    unit: beam radius:mm angle:rad
    bpp(1) : return diffraction limit lambda/pi   
    '''        
    BPP = M2 * (wl*1.e3/pi)

    return BPP
    
def gaussianlens(z,w,f):
    '''
    all units are in mm.    
    
    '''
    w1 = np.sqrt(w**2/((1-z/f)**2+(pi*w**2/(wl*1e3)*f)**2))
    
    z1 = f - (f-z)/((1-z/f)**2+(pi*w**2/(wl*1e3)*f)**2)
    
    return z1,w1
     
# define a function to calculate the beam diameter for a random z
   
    
    
def gaussiantelescope(z,w,f1,f2,d):
          
    z1prime,w1 = gaussianlens(z,w,f1)
     
    z2 = f1+f2+d-z1prime
     
    z2prime,w2 = gaussianlens(z2,w1,f2)

    theta2 = bpp(M2)/w2
    
    w_after_chirpmirrors = gaussianbeam(w2,pd-z2prime)
    
    return z2prime,w2,theta2,w_after_chirpmirrors
    
    
     

def gaussian(w,z):
    
    wz = w*np.sqrt(1+((wl*1e3)*z/(pi*w**2))**2)
    
    return wz

def gaussianbeam(w,z):
    '''
    for given pulse parameters, calculate the beam radius at certain place
    '''
    w2 = w*np.sqrt(1.+((wl*1e-3)*z/(w**2*pi))**2)

#    r = z*(1.+((wl*1e3)*z/(w**2*pi))**-2)**2    
    
    return w2

def plotgaussian(w):
    '''
    Plot gaussian beam profile within +-3z_r by given wrist radius
    '''
    z_r = pi*w**2/1.03e-3
    
    z = np.linspace(-3*z_r,3*z_r,100)
    w2 = w*np.sqrt(1.+((wl*1e3)*z/(w**2*pi))**2)

    plt.plot(z,w2,color='blue')
    #plt.show()
    plt.grid(True)
    plt.plot(z,-w2,color='blue')
    plt.axvline(x=z_r,color='black',linestyle="--")
    plt.axvline(x=-z_r,color='black',linestyle="--")
    #plt.show(False)
   
def bppfiber(NA,d):  
    '''
    NA and MFD!
    d unit in um
    BPP unit: mm rad
    '''
    BPP = natoangel(NA) * d/2 * 1e-3
    return BPP

def natoangel(NA):
    
    theta = np.arcsin(NA)     
    
    return theta

def angel(w):
    '''
    unit:mrad
    '''
    theta = (wl*1e3)/(pi*w)*1000

    print (str("%.4f" %theta) + ' mrad')
    return theta    
    
    

def coupling(w,MFD):
    
    f = lensfocus(2*w,MFD)
        
    z,w = gaussianlens(f,w,f)
    
    return z,w,f


def lensfocus(Dbeam63,MFD):
    '''
    Dbeam is the diameter of D63(mm)
    f in mm
    MFD is the Mode Field Diameter of fiber(micron)
    '''
    f = Dbeam63*np.sqrt(2)*MFD*0.7625
    
    return f

def cost(speed, power, gas):
    return ((609.72) + (140 * power) + (gas)) / speed

def material(density, vap_temp, surr_temp, lat_fusion, lat_vap, heat_cap):
    return density * (lat_fusion + lat_vap + heat_cap(vap_temp - surr_temp))

def power_req(density, vap_temp, surr_temp, lat_fusion, lat_vap, heat_cap, speed):
    power = speed * material(density, vap_temp, surr_temp, lat_fusion, lat_vap, heat_cap)
    return power

def speed(power, density, vap_temp, surr_temp, lat_fusion, lat_vap, heat_cap):
    speed_req = power / material(density, vap_temp, surr_temp, lat_fusion, lat_vap, heat_cap)
    return speed_req

def modlding_material(avg_kerf, density, Cp, Lf, Tm, Ts, n, Lv, m):
    return ((avg_kerf / 1000) *  density * ((Cp * (Tm - Ts) + (Lf * 1000) + (m*Lv))))/n

def mol_power(avg_kerf, density, Cp, Lf, Tm, Ts, n, t, v, Lv, m):
    power = (t / 1000) * (v / 60) * modlding_material(avg_kerf, density, Cp, Lf, Tm, Ts, n, Lv, m)
    return power

def mol_velocity(power, avg_kerf, density, Cp, Lf, Tm, Ts, n, t, Lv, m):
    velocity = power / ((t/1000) * modlding_material(avg_kerf, density, Cp, Lf, Tm, Ts, n, Lv, m))
    return velocity * 60

def show(frame, x1=0, y1=0):
    global current
    current.place_forget()
    frame.place(x = x1, y=y1)
    current = frame

def putlabel(label, value, unit):
    label.configure(text = f"{value} {unit}", font = ("Verdana", 16))

def entry(entry_label):
    try:
        return float(entry_label.get())
    except:
        return 0    


class Example(Frame):
    def __init__(self, master, *pargs):
        Frame.__init__(self, master, *pargs)

        self.image = Image.open("2.jpg")
        self.img_copy= self.image.copy()

        self.background_image = ImageTk.PhotoImage(self.image)

        self.background = Label(self, image=self.background_image)
        self.background.pack(fill=BOTH, expand=YES)
        self.background.bind('<Configure>', self._resize_image)

    def _resize_image(self,event):

        new_width = event.width
        new_height = event.height

        self.image = self.img_copy.resize((new_width, new_height))

        self.background_image = ImageTk.PhotoImage(self.image)
        self.background.configure(image =  self.background_image)




root = Tk()

#Main
root.title("Laser Calculator")

try:
    root.iconbitmap("laser_black.ico")
except:
    pass

root.geometry('800x500')
font = ("Verdana", 8)

try:
    e = Example(root)
    e.pack(fill=BOTH, expand=YES)
except:
    pass
#E2-E1
frame_e2_e1 = LabelFrame(root)

btn_e2_e1 = Button(root, text = "Absorption of\nPhoton", command= lambda:show(frame_e2_e1, 270, 200))
btn_e2_e1.place(x=0, y=0)

feq_e2_e1 = Label(frame_e2_e1, text="Frequency (Hz)", fg="black", font=font)
feq_e2_e1.grid(row=0, column=0)

feq_entry_e2 = Entry(frame_e2_e1)
feq_entry_e2.grid(row=0, column=1)


btn_cal_e2_e1 = Button(frame_e2_e1, text = "Calculate", command= lambda:putlabel(cal_e2_e1, e2_e1(entry(feq_entry_e2)), "J"))
btn_cal_e2_e1.grid(row=1, column=1)

cal_e2_e1 = Label(frame_e2_e1)
cal_e2_e1.grid(row=2, column=1)

global current
current = frame_e2_e1
frame_e2_e1.place(x=270, y=200)

#Mode Number

frame_mod_num = LabelFrame(root)

btn_mod_num = Button(root, text = "Mode\nNumber", command= lambda:show(frame_mod_num, 270, 200))
btn_mod_num.place(x=84, y=0)

wav_mod_num = Label(frame_mod_num, text="Wavelength (m)", fg="black", font=font)
wav_mod_num.grid(row=0, column=0)

l_mod_num = Label(frame_mod_num, text="length of optical cavity (m)", fg="black", font=font)
l_mod_num.grid(row=1, column=0)

wav_entry_mod = Entry(frame_mod_num)
wav_entry_mod.grid(row=0, column=1)

l_entry_mod = Entry(frame_mod_num)
l_entry_mod.grid(row=1, column=1)

btn_mod_num = Button(frame_mod_num, text = "Calculate", command= lambda:putlabel(cal_mod_num, mode_num(entry(wav_entry_mod), entry(l_entry_mod)), " "))
btn_mod_num.grid(row=2, column=1)

cal_mod_num = Label(frame_mod_num)
cal_mod_num.grid(row=3, column=1)

#mode_sep_feq

frame_mode_sep = LabelFrame(root)

btn_mode_sep = Button(root, text = "Mode Seperation \n Frequency", command= lambda:show(frame_mode_sep, 270, 200))
btn_mode_sep.place(x=139, y=0)

l_mode_sep = Label(frame_mode_sep, text="Length of Optical\nCavity of Laser (m)", fg="black", font=font)
l_mode_sep.grid(row=0, column=0)

l_entry_sep = Entry(frame_mode_sep)
l_entry_sep.grid(row=0, column=1)

btn_mode_sep = Button(frame_mode_sep, text = "Calculate", command= lambda:putlabel(cal_mode_sep, mode_sep_feq(entry(l_entry_sep)), "Hz"))
btn_mode_sep.grid(row=1, column=1)

cal_mode_sep = Label(frame_mode_sep)
cal_mode_sep.grid(row=2, column=1)

#width_of_curve

frame_width_of_curve = LabelFrame(root)

btn_width_of_curve = Button(root, text = "Width of\nGain Curve", command= lambda:show(frame_width_of_curve, 270, 200))
btn_width_of_curve.place(x=243, y=0)

delv_width_of_curve = Label(frame_width_of_curve, text=" Δv (GHz)", fg="black", font=font)
delv_width_of_curve.grid(row=0, column=0)

feq_width_of_curve = Label(frame_width_of_curve, text="Frequency (Hz)", fg="black", font=font)
feq_width_of_curve.grid(row=1, column=0)

delv_entry_width = Entry(frame_width_of_curve)
delv_entry_width.grid(row=0, column=1)

feq_entry_width = Entry(frame_width_of_curve)
feq_entry_width.grid(row=1, column=1)

btn_width_of_curve = Button(frame_width_of_curve, text = "Calculate", command= lambda:putlabel(cal_width_of_curve, width_of_curve(float(delv_entry_width.get()), float(feq_entry_width.get())), "m"))
btn_width_of_curve.grid(row=2, column=1)

cal_width_of_curve = Label(frame_width_of_curve)
cal_width_of_curve.grid(row=3, column=1)

#relative_band_width

frame_relative_band_width= LabelFrame(root)
btn_relative_band_width = Button(root, text = "Relative\nBand Width", command= lambda:show(frame_relative_band_width, 270, 200))
btn_relative_band_width.place(x=312, y=0)

delv_relative_band_width = Label(frame_relative_band_width, text=" Δv (GHz)", fg="black", font=font)
delv_relative_band_width.grid(row=0, column=0)

feq_relative_band_width = Label(frame_relative_band_width, text="Frequency (Hz)", fg="black", font=font)
feq_relative_band_width.grid(row=1, column=0)

delv_band_entry_band = Entry(frame_relative_band_width)
delv_band_entry_band.grid(row=0, column=1)

feq_band_entry_band = Entry(frame_relative_band_width)
feq_band_entry_band.grid(row=1, column=1)

btn_relative_band_width = Button(frame_relative_band_width, text = "Calculate", command= lambda:putlabel(cal_relative_band_width, relative_band_width(float(delv_band_entry_band.get()), float(feq_band_entry_band.get())), " "))
btn_relative_band_width.grid(row=2, column=1)

cal_relative_band_width = Label(frame_relative_band_width)
cal_relative_band_width.grid(row=3, column=1)

#rayleighrange

frame_rayleighrange= LabelFrame(root)
btn_rayleighrange = Button(root, text = "Rayleigh\nRange", command= lambda:show(frame_rayleighrange, 270, 200))
btn_rayleighrange.place(x=385, y=0)

d_rayleighrange = Label(frame_rayleighrange, text=" Diameter of\nGaussian Beam (mm)", fg="black", font=font)
d_rayleighrange.grid(row=0, column=0)

wav_rayleighrange = Label(frame_rayleighrange, text="Wavelength (mm)", fg="black", font=font)
wav_rayleighrange.grid(row=1, column=0)

d_entry_dir = Entry(frame_rayleighrange)
d_entry_dir.grid(row=0, column=1)

wav_entry_dir = Entry(frame_rayleighrange)
wav_entry_dir.grid(row=1, column=1)

btn_rayleighrange = Button(frame_rayleighrange, text = "Calculate", command= lambda:putlabel(cal_rayleighrange, rayleighrang(float(d_entry_dir.get()), float(wav_entry_dir.get())), "mm"))
btn_rayleighrange.grid(row=2, column=1)

cal_rayleighrange = Label(frame_rayleighrange)
cal_rayleighrange.grid(row=3, column=1)

#semi_angle_diffraction

frame_semi_angle= LabelFrame(root)
btn_semi_angle = Button(root, text = "Semi Angle\nof Diffraction", command= lambda:show(frame_semi_angle, 270, 200))
btn_semi_angle.place(x=441, y=0)

d_semi_angle = Label(frame_semi_angle, text=" Diameter of\nGaussian Beam (m)", fg="black", font=font)
d_semi_angle.grid(row=0, column=0)

wav_semi_angle = Label(frame_semi_angle, text="Wavelength (m)", fg="black", font=font)
wav_semi_angle.grid(row=1, column=0)

d_entry_semi = Entry(frame_semi_angle)
d_entry_semi.grid(row=0, column=1)

wav_entry_semi = Entry(frame_semi_angle)
wav_entry_semi.grid(row=1, column=1)

btn_semi_angle = Button(frame_semi_angle, text = "Calculate", command= lambda:putlabel(cal_semi_angle, semi_angle_diffraction(float(d_entry_semi.get()), float(wav_entry_semi.get())), "rad"))
btn_semi_angle.grid(row=2, column=1)

cal_semi_angle = Label(frame_semi_angle)
cal_semi_angle.grid(row=3, column=1)

#intensity

frame_intensity= LabelFrame(root)
btn_intensity = Button(root, text = "Intensity\n ", command= lambda:show(frame_intensity, 270, 200))
btn_intensity.place(x=522, y=0)

power_intensity = Label(frame_intensity, text=" Power (W)", fg="black", font=font)
power_intensity.grid(row=0, column=0)

diameter_intensity = Label(frame_intensity, text="Diameter (m)", fg="black", font=font)
diameter_intensity.grid(row=1, column=0)

power_entry_inten = Entry(frame_intensity)
power_entry_inten.grid(row=0, column=1)

diameter_entry_inten = Entry(frame_intensity)
diameter_entry_inten.grid(row=1, column=1)

btn_intensity = Button(frame_intensity, text = "Calculate", command= lambda:putlabel(cal_intensity, intensity(float(power_entry_inten.get()), float(diameter_entry_inten.get())), "W/m²"))
btn_intensity.grid(row=2, column=1)

cal_intensity = Label(frame_intensity)
cal_intensity.grid(row=3, column=1)

#beam_waist

frame_beam_waist = LabelFrame(root)

btn_beam_waist = Button(root, text = "Beam\nWaist", command= lambda:show(frame_beam_waist, 270, 200))
btn_beam_waist.place(x=578, y=0)

focal_beam_waist = Label(frame_beam_waist, text=" Focal Length(mm)", fg="black", font=font)
focal_beam_waist.grid(row=0, column=0)

distance_beam_waist = Label(frame_beam_waist, text=" Distance (mm)", fg="black", font=font)
distance_beam_waist.grid(row=1, column=0)

wav_beam_waist = Label(frame_beam_waist, text="Wavelength (mm)", fg="black", font=font)
wav_beam_waist.grid(row=2, column=0)

focal_entry_beam_waist = Entry(frame_beam_waist)
focal_entry_beam_waist.grid(row=0, column=1)

distance_entry_beam_waist = Entry(frame_beam_waist)
distance_entry_beam_waist.grid(row=1, column=1)

wav_entry_beam_waist = Entry(frame_beam_waist)
wav_entry_beam_waist.grid(row=2, column=1)

btn_beam_waist = Button(frame_beam_waist, text = "Calculate", command= lambda:putlabel(cal_beam_waist, beam_waist(float(focal_entry_beam_waist.get()), float(distance_entry_beam_waist.get()), float(wav_entry_beam_waist.get())), "mm"))
btn_beam_waist.grid(row=3, column=1)

cal_beam_waist = Label(frame_beam_waist)
cal_beam_waist.grid(row=4, column=1)

#Interaction Time       

frame_interaction= LabelFrame(root)
btn_interaction = Button(root, text = "Interaction\nTime", command= lambda:show(frame_interaction, 270, 200))
btn_interaction.place(x=619, y=0)

radius_interaction = Label(frame_interaction, text=" Radius (m)", fg="black", font=font)
radius_interaction.grid(row=0, column=0)

velocity_interaction = Label(frame_interaction, text="Velocity (m/s)", fg="black", font=font)
velocity_interaction.grid(row=1, column=0)

radius_entry_interaction = Entry(frame_interaction)
radius_entry_interaction.grid(row=0, column=1)

velocity_entry_interaction = Entry(frame_interaction)
velocity_entry_interaction.grid(row=1, column=1)

btn_interaction = Button(frame_interaction, text = "Calculate", command= lambda:putlabel(cal_interaction, interaction_time(float(radius_entry_interaction.get()), float(velocity_entry_interaction.get())), " s"))
btn_interaction.grid(row=2, column=1)

cal_interaction = Label(frame_interaction)
cal_interaction.grid(row=3, column=1)

#Molding

frame_molding = LabelFrame(root)

btn_molding = Button(root, text = "Laser\n", command= lambda:show(frame_molding, 200, 150))
btn_molding.place(x=687, y=0)

kerf_molding = Label(frame_molding, text=" Average kerf width (mm)", fg="black", font=font)
kerf_molding.grid(row=0, column=0)

heat_fusion_molding = Label(frame_molding, text="Latent heat of fusion (kJ/kg)", fg="black", font=font)
heat_fusion_molding.grid(row=1, column=0)

heat_vapour_molding = Label(frame_molding, text="Latent heat of vaporization (kJ/kg)", fg="black", font=font)
heat_vapour_molding.grid(row=2, column=0)

density_molding = Label(frame_molding, text="Density (kg/m3)", fg="black", font=font)
density_molding.grid(row=3, column=0)

current_temp_molding = Label(frame_molding, text="Specific Heat", fg="black", font=font)
current_temp_molding.grid(row=4, column=0)

fraction_melt_molding = Label(frame_molding, text="Fraction of the melt vaporized ", fg="black", font=font)
fraction_melt_molding.grid(row=5, column=0)

melt_temp_molding = Label(frame_molding, text="Melting Point of material (Celsius)", fg="black", font=font)
melt_temp_molding.grid(row=6, column=0)

cup_molding = Label(frame_molding, text="Cupling Coefficient", fg="black", font=font)
cup_molding.grid(row=7, column=0)

thickness_molding = Label(frame_molding, text="Thickness (mm)", fg="black", font=font)
thickness_molding.grid(row=8, column=0)

kerf_entry_molding = Entry(frame_molding)
kerf_entry_molding.grid(row=0, column=1)

fusion_entry_molding = Entry(frame_molding)
fusion_entry_molding.grid(row=1, column=1)

vapour_entry_molding = Entry(frame_molding)
vapour_entry_molding.grid(row=2, column=1)

density_entry_molding = Entry(frame_molding)
density_entry_molding.grid(row=3, column=1)

specific_entry_molding = Entry(frame_molding)
specific_entry_molding.grid(row=4, column=1)

fraction_entry_molding = Entry(frame_molding)
fraction_entry_molding.grid(row=5, column=1)

melt_temp_entry_molding = Entry(frame_molding)
melt_temp_entry_molding.grid(row=6, column=1)

cup_entry_molding = Entry(frame_molding)
cup_entry_molding.grid(row=7, column=1)

thickness_entry_molding = Entry(frame_molding)
thickness_entry_molding.grid(row=8, column=1)

#Power
def radio(frame1, frame2):
    try:
        frame1.grid_forget()
    except:
        pass
    frame2.grid(row=10, column=1)

power = Radiobutton(frame_molding, command=lambda:radio(newframe2, newframe1), text="Power", value = 1)
power.grid(row=9, column=0)

newframe1 = LabelFrame(frame_molding)

velocity_molding = Label(newframe1, text="velocity (min/s)", fg="black", font=font)
velocity_molding.grid(row=0, column=0)

velocity_molding = Entry(newframe1)
velocity_molding.grid(row=0, column=1)

btn_molding = Button(newframe1, text = "Calculate", command= lambda:putlabel(cal_molding, mol_power(avg_kerf = float(kerf_entry_molding.get()),density = float(density_entry_molding.get()), Cp = float(specific_entry_molding.get()), Lf = float(fusion_entry_molding.get()), Tm = float(melt_temp_entry_molding.get()), Ts = 27.0, n = float(cup_entry_molding.get()), t = float(thickness_entry_molding.get()), v = float(velocity_molding.get()), Lv = float(vapour_entry_molding.get()), m = float(fraction_entry_molding.get())), " Watts"))
btn_molding.grid(row=2, column=1)

cal_molding = Label(newframe1)
cal_molding.grid(row=3, column=1)

#Time

velocity = Radiobutton(frame_molding, command=lambda:radio(newframe1, newframe2), text="Velocity", value = 2)
velocity.grid(row=9, column=1)

newframe2 = LabelFrame(frame_molding)
newframe2.grid(row=10, column=1)

power_molding = Label(newframe2, text="Power", fg="black", font=font)
power_molding.grid(row=8, column=0)

power_molding = Entry(newframe2)
power_molding.grid(row=8, column=1)

btn2_molding = Button(newframe2, text = "Calculate", command= lambda:putlabel(cal2_molding, mol_velocity(power = float(power_molding.get()), avg_kerf = float(kerf_entry_molding.get()),density = float(density_entry_molding.get()), Cp = float(specific_entry_molding.get()), Lf = float(fusion_entry_molding.get()), Tm = float(melt_temp_entry_molding.get()), Ts = 27.0, n = float(cup_entry_molding.get()), t = float(thickness_entry_molding.get()), Lv = float(vapour_entry_molding.get()), m = float(fraction_entry_molding.get())), " m/min"))
btn2_molding.grid(row=9, column=1)

cal2_molding = Label(newframe2)
cal2_molding.grid(row=10, column=1)

#Cost

frame_cost = LabelFrame(root)

btn_cost = Button(root, text = "Cost\n", command= lambda:show(frame_cost, 270, 200))
btn_cost.place(x=725, y=0)

speed_cost = Label(frame_cost, text=" Cutting Speed", fg="black", font=font)
speed_cost.grid(row=0, column=0)

power_cost = Label(frame_cost, text=" Laser Power", fg="black", font=font)
power_cost.grid(row=1, column=0)

gas_cost = Label(frame_cost, text="Cost of assist gas", fg="black", font=font)
gas_cost.grid(row=2, column=0)

speed_entry_cost = Entry(frame_cost)
speed_entry_cost.grid(row=0, column=1)

power_entry_cost = Entry(frame_cost)
power_entry_cost.grid(row=1, column=1)

gas_entry_cost = Entry(frame_cost)
gas_entry_cost.grid(row=2, column=1)

btn_cost = Button(frame_cost, text = "Calculate", command= lambda:putlabel(cal_cost, cost(float(speed_entry_cost.get()), float(power_entry_cost.get()), float(gas_entry_cost.get())), " "))
btn_cost.grid(row=3, column=1)

cal_cost = Label(frame_cost)
cal_cost.grid(row=4, column=1)

#btn_molding = Button(frame_molding, text = "calculate", command= lambda:putlabel(cal_beam_waist, beam_waist(float(focal_entry_beam_waist.get()), float(distance_entry_beam_waist.get()), float(wav_entry_beam_waist.get()))))
#btn_molding.grid(row=5, column=1)

#cal_molding = Label(frame_molding)
#cal_molding.grid(row=6, column=1)

root.mainloop()