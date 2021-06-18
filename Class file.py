
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate as interpol
import scipy.integrate as si

amu = 1.67377e-27
boltzmann_const = 1.38064852e-23
Electron_mass_kg = 9.10938e-31 #Mass of an electron [kg]
Electron_charge = 1.60217662e-19 #Elemetry charge [C]

dir = "C:/Users/tfmo1/OneDrive - University of Southampton/Documents/HET Code/"

class propellant:
    def __init__(self, Mass_kg, Cross_sectional_area_file):
        self.mass = Mass_kg
        self.mass_amu = Mass_kg/amu
        df = pd.DataFrame(pd.read_excel(dir + Cross_sectional_area_file+".xlsx"), \
                          columns=['Tev','Sigma'])
        self.cross_sectional_area_all = \
            [np.array(df['Tev']),np.array(df['Sigma'])]
    
    
    def cross_sectional_area(self,Electron_temperature):
        if Electron_temperature>2000:
            return 0
        Te_file = self.cross_sectional_area_all[0]
        Sigma_file = self.cross_sectional_area_all[1]
        xy = interpol.interp1d(Te_file,Sigma_file,fill_value="extrapolate")
        return(xy(Electron_temperature)*1e-20)

    def reaction_rate(self, Electron_temperature):
        def Electron_velocity_DF(w):
            Te = Electron_temperature
            a = np.sqrt(2/np.pi)
            b = (Electron_mass_kg/(Electron_charge * Te))**(3/2)
            c = w**2
            d = np.exp(-((Electron_mass_kg * (w**2))/(2 * Electron_charge * Te)))
            e = a * b * c * d
            return(e)
    
        
        def Boltzman_dis_average(w):
            EVDF = Electron_velocity_DF(w)
            Te_v = ((w**2) * Electron_mass_kg* np.pi) / (8 * Electron_charge)
            Sigma = self.cross_sectional_area(Te_v)
            return(w * EVDF * Sigma)
        RR,error = si.quad(Boltzman_dis_average,0,2e8)
        return(RR)
    
    def thermal_velocity(self,Temperature):
        return(np.sqrt((8 * boltzmann_const* Temperature)/(np.pi * self.mass)))


    def Scale_Target_thrust(self,Target_thrust,Target_voltage=300):
        C_T_1 = 829.725
        C_T_2 = 0.80994
        C_isp = 85.617
        C_Id = 907.3378
        C_P = 0.01872608
        C_hd =  0.256
        C_B_1 = 2.68882E-05
        C_B_2 = 1.378610642
    
    
        m_dot = Target_thrust / (C_T_1 * np.sqrt(Target_voltage))
        d = np.sqrt(Target_thrust /(C_T_2 * np.sqrt(Target_voltage)))
        I_d = C_Id * d**2
        P = I_d * Target_voltage
        h = C_hd * d
        
        A = d * h * np.pi
        v_n = self.thermal_velocity(800)
        n_n = m_dot / (A * v_n * self.mass)
        
        
        B_2 = C_B_2 * ((m_dot * np.sqrt(Target_voltage))/(d**2))
        L = C_B_1 * ((np.sqrt(Target_voltage))/(B_2))
        
        Isp = round(Target_thrust / (9.81*m_dot),2)
        d = round(d * 1e3,2)
        h = round(h * 1e3,2)
        L = round(L * 1e3,2)
        Target_thrust = round(Target_thrust * 1e3,2)
        
        m_dot = round(m_dot * 1e6,2)
        m_dot_sccm = m_dot/( 7.43583e-4 * self.mass /amu)
        B_2 = round(B_2,5)
        print(f"Power of {P} W, mean diam of {d} mm, channel width of {h} mm, channel length of {L} mm, Thrust of {Target_thrust} mN, Specific impulse of {Isp}, Field strength of {B_2} T, and mass flow rate of {m_dot} mg/s ({m_dot_sccm} sccm).")

        
    def Scale_Target_power(self,Target_power,Target_voltage):
        C_T_1 = 829.725
        C_T_2 = 0.80994
        C_isp = 85.617
        C_Id = 907.3378
        C_P = 0.01872608
        C_hd =  0.256
        C_B_1 = 2.68882E-05
        C_B_2 = 1.378610642


        I_d = Target_power / Target_voltage
        d = np.sqrt(I_d/C_Id)
        T = C_T_2 * np.sqrt(Target_voltage) * d**2
        m_dot = T / (C_T_1 * np.sqrt(Target_voltage))
        P = I_d * Target_voltage
        h = C_hd * d
        
        A = d * h * np.pi
        n_n = 1.6e19
        v_n = self.thermal_velocity(800)

        m_dot = n_n * (A * v_n * self.mass)
        #n_n = m_dot / (A * v_n * self.mass)
        
        B_2 = C_B_2 * ((m_dot * np.sqrt(Target_voltage))/(d**2))
        L = C_B_1 * ((np.sqrt(Target_voltage))/(B_2))
        
        Isp = round(T / (9.81*m_dot),2)
        d = round(d * 1e3,2)
        h = round(h * 1e3,2)
        L = round(L * 1e3,2)
        Target_thrust = round(T * 1e3,2)
        m_dot = round(m_dot * 1e6,2)
        m_dot_sccm = m_dot/( 7.43583e-4 * self.mass /amu)
        B_2 = round(B_2,5)
        
        print(f"Power of {P} W, mean diam of {d} mm, channel width of {h} mm, channel length of {L} mm, Thrust of {Target_thrust} mN, Specific impulse of {Isp}, Field strength of {B_2} T, and mass flow rate of {m_dot} mg/s ({m_dot_sccm} sccm).")




Xenon = propellant(131.293*amu,"Xenon1")
Krypton = propellant(83.798*amu,"Krypton1")
Argon = propellant(39.948*amu,"Argon1wetzel")
Neon = propellant(20.1797*amu,"Neon1wetzel")
Oxy = propellant(2*15.999*amu,"Oxygen1")
Nitro = propellant(2*14.007*amu,"Oxygen1")

def Scale_Target_thrust(Target_thrust,Target_voltage=300):
    C_T_1 = 829.725
    C_T_2 = 0.80994
    C_isp = 85.617
    C_Id = 907.3378
    C_P = 0.01872608
    C_hd =  0.256
    C_B_1 = 2.68882E-05
    C_B_2 = 1.378610642

    m_dot = Target_thrust / (C_T_1 * np.sqrt(Target_voltage))
    d = np.sqrt(Target_thrust /(C_T_2 * np.sqrt(Target_voltage)))
    I_d = C_Id * d**2
    P = I_d * Target_voltage
    h = C_hd * d
    
    A = d * h * np.pi
    v_n = 360
    n_n = m_dot / (A * v_n * 131.293*amu)
    
    B_2 = C_B_2 * ((m_dot * np.sqrt(Target_voltage))/(d**2))
    L = C_B_1 * ((np.sqrt(Target_voltage))/(B_2))
    
    Isp = round(Target_thrust / (9.81*m_dot),2)
    d = round(d * 1e3,2)
    h = round(h * 1e3,2)
    L = round(L * 1e3,2)
    Target_thrust = round(Target_thrust * 1e3,2)
    m_dot = round(m_dot * 1e6,2)
    B_2 = round(B_2,5)
    
    print(f"Power of {P} W, mean diam of {d} mm, channel width of {h} mm, channel length of {L} mm, Thrust of {Target_thrust} mN, Specific impulse of {Isp}, Field strength of {B_2} T, and mass flow rate of {m_dot} mg/s.")



def Scale_Target_power(Target_power,Target_voltage):
    C_T_1 = 829.725
    C_T_2 = 0.80994
    C_isp = 85.617
    C_Id = 907.3378
    C_P = 0.01872608
    C_hd =  0.256
    C_B_1 = 2.68882E-05
    C_B_2 = 1.378610642


    I_d = Target_power / Target_voltage
    d = np.sqrt(I_d/C_Id)
    T = C_T_2 * np.sqrt(Target_voltage) * d**2
    m_dot = T / (C_T_1 * np.sqrt(Target_voltage))
    P = I_d * Target_voltage
    h = C_hd * d
    
    A = d * h * np.pi
    v_n = 260
    n_n = m_dot / (A * v_n * 131.293*amu)
    
    B_2 = C_B_2 * ((m_dot * np.sqrt(Target_voltage))/(d**2))
    L = C_B_1 * ((np.sqrt(Target_voltage))/(B_2))
    
    Isp = round(T / (9.81*m_dot),2)
    d = round(d * 1e3,2)
    h = round(h * 1e3,2)
    L = round(L * 1e3,2)
    Target_thrust = round(T * 1e3,2)
    m_dot = round(m_dot * 1e6,2)
    B_2 = round(B_2,5)
    
    print(f"Power of {P} W, mean diam of {d} mm, channel width of {h} mm, channel length of {L} mm, Thrust of {Target_thrust} mN, Specific impulse of {Isp}, Field strength of {B_2} T, and mass flow rate of {m_dot} mg/s.")



print( Xenon.Scale_Target_power(100,300) )
print( Krypton.Scale_Target_power(100,300) )
print( Argon.Scale_Target_power(100,300) )
print( Neon.Scale_Target_power(100,300) )
print( Oxy.Scale_Target_power(100,300) )
print( Nitro.Scale_Target_power(100,300) )