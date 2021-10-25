import cantera as ct
from numpy.polynomial import polynomial as P
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
import sys
import csv
import os
from datetime import *
import time

# using time diffrence in saving csv file each time we run a loop
current_day = datetime.now()
Time = str(current_day.strftime('%H_%M_%S'))

# auto generating dir in respect to date for easily compiling of daily data
data = str(date.today())
if not os.path.exists(str(data)):
    os.makedirs(str(data))

# Numbers of case
case1 = '0.002538447, O2:0.208384625, N2:0.783923112'
case2 = '0.002538450, O2:0.111114625, N2:0.000222112'
case3 = '0.002538450, N2:0.008384625, O2:0.000923112'

T1 = 388.0
T2 = 398.0
T3 = 425.0

P1= 0.766301
P2 = 0.99999
P3 = 1.026586

Tcases = [T1, T2, T3]
cases = [case1, case2, case3]
Pcase = [P1, P2, P3]

for case in cases:
    for TC in Tcases:
        for PC in Pcase:
            if case == case1 and TC == T1 and PC == P1 :



                print('Case One Mole :', case, 'Temperature One :', TC, 'Pressure: ', PC)
                time.sleep(5)
                kase = str(case)
                # Read in polynomial coefficients for the expansion trace
                a_input = np.loadtxt('expfit_rev.dat')
                a = a_input.tolist()

                t = 0.0
                n = 0.0
                n_steps = 5000
                dt = 1e-05
                stroke = 0.1422
                cl = 0.0075               
                dia = 0.04
                vadd = 0.00002311
                taccel = 0.020501
                tdecel = 0.0069

                tcomp = 0.030
                tconst = tcomp - taccel - tdecel
                t2 = taccel + tconst
                volstart = (3.14/4.0 * dia**2)*(stroke+cl) + vadd
                vmax = stroke/(((taccel+tdecel)/2.0)+ tconst)

                # velocity for three stages of compression
                v_A = -3.14/4.0 *(dia**2)* (vmax * t/taccel)
                v_C = -3.14/4.0 *((dia**2) * vmax)
                v_D = -3.14/4.0 *(dia**2)*(vmax - vmax/tdecel * (t-t2))


                outfile = open(data+'/'+str(TC) +'.'+str(Time)+'comp.csv', 'w')
                csvfile = csv.writer(outfile)
                csvfile.writerow(['time(s)','temperature (K)','pressure(Bar)','volume (m3)','velocity (m/s)'])
                print('%10s %10s %10s  %10s %14s' % ('t [s]','T [K]','P [bar]','vol [m3]', 'vel [m/s]'))

                tout = []
                vout = []
                velout = []
                Tout = []
                Pout = []
                vtdc = []

                #Take in the chemical kinetics mechanism
                gas = ct.Solution('PRF_171.cti')
                air = ct.Solution('air.xml')
                gas.TPX = TC, PC*ct.one_atm,   'C2H5O2:{Mole}'.format(Mole=case)
                print('Temperature:', TC, 'MOLE: ', case, 'PRESSUE: ', PC)
                time.sleep(10)
                r = ct.IdealGasReactor(gas, energy='on', volume = volstart) 
                env = ct.Reservoir(contents=gas, name='environment')
                wall = ct.Wall(r, env,  velocity=v_A)
                sim = ct.ReactorNet([r])


                for n in range(n_steps): 
                    if t <= taccel:
                        wall.set_velocity(v_A)
                        t += dt
                        sim.advance(t)
                        disp = vmax * (t**2/2.0/taccel)
                        v_A = -3.14/4.0 *(dia**2)* (vmax * t/taccel)
                        volume = volstart - (3.14/4.0 *(dia**2)*(disp))
                        tout.append(t)
                        Tout.append(r.T)
                        Pout.append(1.e-5*r.thermo.P)
                        vout.append(r.volume)
                        velout.append(1000*v_A)
                        csvfile = csv.writer(outfile)
                        csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, 1000*v_A])
                        print('%10.3f %10.3f %10.3f %10.3f %10.3f' % (t, r.T, 1.e-5*r.thermo.P, r.volume, v_A))
                        
                        
                    
                    elif t > taccel and t <= t2:
                        wall.set_velocity(v_C)
                        t += dt
                        sim.advance(t)
                        disp = vmax * taccel/2.0 + (vmax *(t-taccel))
                        v_C = -3.14/4.0 *((dia**2) * vmax)
                        volume = volstart - (3.14/4.0 *(dia**2)*(disp))
                        tout.append(t)
                        Tout.append(r.T)
                        Pout.append(1.e-5*r.thermo.P)
                        vout.append(r.volume)
                        velout.append(1000*v_C)
                        csvfile = csv.writer(outfile)
                        csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, 1000*v_C])
                        print('%10.3e %10.3f %10.3f %10.3f %10.3f' % (t, r.T, 1.e-5*r.thermo.P, r.volume, v_C))
                    
                    elif t > t2 and t <= tcomp:
                        wall.set_velocity(v_D)
                        t += dt
                        sim.advance(t)
                        disp = (vmax * taccel/2.0) + (vmax * tconst)
                        disp = (disp + vmax *(t-t2)) - (vmax *(t-t2)**2/tdecel/2.0)
                        v_D = -3.14/4.0 *(dia**2)*(vmax - vmax/tdecel * (t-t2))
                        volume = volstart - (3.14/4.0 *(dia**2)*(disp))
                        tout.append(t)
                        Tout.append(r.T)
                        Pout.append(1.e-5*r.thermo.P)
                        vout.append(r.volume)
                        velout.append(1000*v_D)
                        csvfile = csv.writer(outfile)
                        csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, 1000*v_D])
                        print('%10.3e %10.3f %10.3f %10.3f %10.3f' % (t, r.T, 1.e-5*r.thermo.P, r.volume, v_D))
                        vtdc = volume
                        
                        

                    elif  t > tcomp:
                        sim.set_initial_time(t)
                        ta = (t - tcomp)
                        poly = a[0] + a[1]*ta + a[2]*ta**2 + a[3]*ta**3 + a[4]*ta**4 + a[5]*ta**5 + a[6]*ta**6 + a[7]*ta**7 + a[8]*ta**8 + a[9]*ta**9 + a[10]*ta**10 
                        volume =  P.polymul(vtdc, poly)    
                        dv = a[1] + 2*a[2]*ta + 3*a[3]*ta**2 + 4*a[4]*ta**3 +  5*a[5]*ta**4 +  6*a[6]*ta**5 + 7*a[7]*ta**6 + 8*a[8]*ta**7 +  9*a[9]*ta**8 +  10*a[10]*ta**9    
                        v_E = P.polymul(vtdc, dv)
                        wall.set_velocity(v_E)
                        t += dt
                        sim.advance(t) 
                        tout.append(t)
                        Tout.append(r.T)
                        Pout.append(1.e-5*r.thermo.P)
                        vout.append(r.volume)
                        velout.append(v_E)
                        csvfile = csv.writer(outfile)
                        csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, v_E])
                        print('%10.3f %10.3f %10.3f %10.3f %10.3f' %(t, r.T, 1.e-5*r.thermo.P, r.volume, v_E))
                        

                outfile.close()
                print('Output written to file result_vel.csv')
                print('Directory: '+os.getcwd())
                print (vmax)          
                print (volstart)
                print(vtdc)


                print('Directory: '+os.getcwd())
                fig = plt.figure()
                #plt.clf()
                plt.subplot(2,2,1)
                plt.plot(tout, vout)
                plt.xlabel('Time (s)')
                plt.ylabel('volume (m^3)')
                plt.grid(False)

                plt.subplot(2,2,2)
                plt.plot(tout, velout)
                plt.xlabel('Time (s)')
                plt.ylabel('velocity (m/s)')
                plt.grid(False)

                plt.subplot(2,2,3)
                plt.plot(tout, Tout)
                plt.xlabel('Time (s)')
                plt.ylabel('Temperature (K)')
                plt.grid(False)

                plt.subplot(2,2,4)
                plt.plot(tout, Pout)
                plt.xlabel('Time (s)')
                plt.ylabel('Pressure (bar)')
                plt.grid(False)
                plt.tight_layout()
                fig.savefig(str(TC) +'.'+str(Time)+'full_figure.png')
                plt.show()

            elif case == case2 and TC == T2 and PC == P2:

                print('Case Two Mole:', case2, 'Temperature Two: ', T2, 'Pressure Two: ',P2)
                time.sleep(5)

                # Read in polynomial coefficients for the expansion trace
                a_input = np.loadtxt('expfit_rev.dat')
                a = a_input.tolist()

                t = 0.0
                n = 0.0
                n_steps = 5000
                dt = 1e-05
                stroke = 0.1422
                cl = 0.0075               
                dia = 0.04
                vadd = 0.00002311
                taccel = 0.020501
                tdecel = 0.0069

                tcomp = 0.030
                tconst = tcomp - taccel - tdecel
                t2 = taccel + tconst
                volstart = (3.14/4.0 * dia**2)*(stroke+cl) + vadd
                vmax = stroke/(((taccel+tdecel)/2.0)+ tconst)

                # velocity for three stages of compression
                v_A = -3.14/4.0 *(dia**2)* (vmax * t/taccel)
                v_C = -3.14/4.0 *((dia**2) * vmax)
                v_D = -3.14/4.0 *(dia**2)*(vmax - vmax/tdecel * (t-t2))


                outfile = open(data+'/'+str(TC) +'.'+str(Time)+'comp.csv', 'w')
                csvfile = csv.writer(outfile)
                csvfile.writerow(['time(s)','temperature (K)','pressure(Bar)','volume (m3)','velocity (m/s)'])
                print('%10s %10s %10s  %10s %14s' % ('t [s]','T [K]','P [bar]','vol [m3]', 'vel [m/s]'))

                tout = []
                vout = []
                velout = []
                Tout = []
                Pout = []
                vtdc = []

                #Take in the chemical kinetics mechanism
                gas = ct.Solution('PRF_171.cti')
                air = ct.Solution('air.xml')
                gas.TPX = TC, PC*ct.one_atm,   'C2H5O2:{Mole}'.format(Mole=case)
                print('Case TWO :', TC, PC, case)
                time.sleep(7)
                r = ct.IdealGasReactor(gas, energy='on', volume = volstart) 
                env = ct.Reservoir(contents=gas, name='environment')
                wall = ct.Wall(r, env,  velocity=v_A)
                sim = ct.ReactorNet([r])


                for n in range(n_steps): 
                    if t <= taccel:
                        wall.set_velocity(v_A)
                        t += dt
                        sim.advance(t)
                        disp = vmax * (t**2/2.0/taccel)
                        v_A = -3.14/4.0 *(dia**2)* (vmax * t/taccel)
                        volume = volstart - (3.14/4.0 *(dia**2)*(disp))
                        tout.append(t)
                        Tout.append(r.T)
                        Pout.append(1.e-5*r.thermo.P)
                        vout.append(r.volume)
                        velout.append(1000*v_A)
                        csvfile = csv.writer(outfile)
                        csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, 1000*v_A])
                        print('%10.3f %10.3f %10.3f %10.3f %10.3f' % (t, r.T, 1.e-5*r.thermo.P, r.volume, v_A))
                        
                        
                    
                    elif t > taccel and t <= t2:
                        wall.set_velocity(v_C)
                        t += dt
                        sim.advance(t)
                        disp = vmax * taccel/2.0 + (vmax *(t-taccel))
                        v_C = -3.14/4.0 *((dia**2) * vmax)
                        volume = volstart - (3.14/4.0 *(dia**2)*(disp))
                        tout.append(t)
                        Tout.append(r.T)
                        Pout.append(1.e-5*r.thermo.P)
                        vout.append(r.volume)
                        velout.append(1000*v_C)
                        csvfile = csv.writer(outfile)
                        csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, 1000*v_C])
                        print('%10.3e %10.3f %10.3f %10.3f %10.3f' % (t, r.T, 1.e-5*r.thermo.P, r.volume, v_C))
                    
                    elif t > t2 and t <= tcomp:
                        wall.set_velocity(v_D)
                        t += dt
                        sim.advance(t)
                        disp = (vmax * taccel/2.0) + (vmax * tconst)
                        disp = (disp + vmax *(t-t2)) - (vmax *(t-t2)**2/tdecel/2.0)
                        v_D = -3.14/4.0 *(dia**2)*(vmax - vmax/tdecel * (t-t2))
                        volume = volstart - (3.14/4.0 *(dia**2)*(disp))
                        tout.append(t)
                        Tout.append(r.T)
                        Pout.append(1.e-5*r.thermo.P)
                        vout.append(r.volume)
                        velout.append(1000*v_D)
                        csvfile = csv.writer(outfile)
                        csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, 1000*v_D])
                        print('%10.3e %10.3f %10.3f %10.3f %10.3f' % (t, r.T, 1.e-5*r.thermo.P, r.volume, v_D))
                        vtdc = volume
                        
                        

                    elif  t > tcomp:
                        sim.set_initial_time(t)
                        ta = (t - tcomp)
                        poly = a[0] + a[1]*ta + a[2]*ta**2 + a[3]*ta**3 + a[4]*ta**4 + a[5]*ta**5 + a[6]*ta**6 + a[7]*ta**7 + a[8]*ta**8 + a[9]*ta**9 + a[10]*ta**10 
                        volume =  P.polymul(vtdc, poly)    
                        dv = a[1] + 2*a[2]*ta + 3*a[3]*ta**2 + 4*a[4]*ta**3 +  5*a[5]*ta**4 +  6*a[6]*ta**5 + 7*a[7]*ta**6 + 8*a[8]*ta**7 +  9*a[9]*ta**8 +  10*a[10]*ta**9    
                        v_E = P.polymul(vtdc, dv)
                        wall.set_velocity(v_E)
                        t += dt
                        sim.advance(t) 
                        tout.append(t)
                        Tout.append(r.T)
                        Pout.append(1.e-5*r.thermo.P)
                        vout.append(r.volume)
                        velout.append(v_E)
                        csvfile = csv.writer(outfile)
                        csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, v_E])
                        print('%10.3f %10.3f %10.3f %10.3f %10.3f' %(t, r.T, 1.e-5*r.thermo.P, r.volume, v_E))
                        

                outfile.close()
                print('Output written to file result_vel.csv')
                print('Directory: '+os.getcwd())
                print (vmax)          
                print (volstart)
                print(vtdc)


                print('Directory: '+os.getcwd())
                fig = plt.figure()
                #plt.clf()
                plt.subplot(2,2,1)
                plt.plot(tout, vout)
                plt.xlabel('Time (s)')
                plt.ylabel('volume (m^3)')
                plt.grid(False)

                plt.subplot(2,2,2)
                plt.plot(tout, velout)
                plt.xlabel('Time (s)')
                plt.ylabel('velocity (m/s)')
                plt.grid(False)

                plt.subplot(2,2,3)
                plt.plot(tout, Tout)
                plt.xlabel('Time (s)')
                plt.ylabel('Temperature (K)')
                plt.grid(False)

                plt.subplot(2,2,4)
                plt.plot(tout, Pout)
                plt.xlabel('Time (s)')
                plt.ylabel('Pressure (bar)')
                plt.grid(False)
                plt.tight_layout()
                fig.savefig(str(TC) +'.'+str(Time)+'full_figure.png')
                plt.show()

            elif case == case3 and TC == T3 and PC == P3:

                        print('Case Three Mole:', case3, 'Temperature Three: ', T2, 'Pressure Three', P3)
                        time.sleep(5)

                        # Read in polynomial coefficients for the expansion trace
                        a_input = np.loadtxt('expfit_rev.dat')
                        a = a_input.tolist()

                        t = 0.0
                        n = 0.0
                        n_steps = 5000
                        dt = 1e-05
                        stroke = 0.1422
                        cl = 0.0075               
                        dia = 0.04
                        vadd = 0.00002311
                        taccel = 0.020501
                        tdecel = 0.0069

                        tcomp = 0.030
                        tconst = tcomp - taccel - tdecel
                        t2 = taccel + tconst
                        volstart = (3.14/4.0 * dia**2)*(stroke+cl) + vadd
                        vmax = stroke/(((taccel+tdecel)/2.0)+ tconst)

                        # velocity for three stages of compression
                        v_A = -3.14/4.0 *(dia**2)* (vmax * t/taccel)
                        v_C = -3.14/4.0 *((dia**2) * vmax)
                        v_D = -3.14/4.0 *(dia**2)*(vmax - vmax/tdecel * (t-t2))


                        outfile = open(data+'/'+str(TC) +'.'+str(Time)+'comp.csv', 'w')
                        csvfile = csv.writer(outfile)
                        csvfile.writerow(['time(s)','temperature (K)','pressure(Bar)','volume (m3)','velocity (m/s)'])
                        print('%10s %10s %10s  %10s %14s' % ('t [s]','T [K]','P [bar]','vol [m3]', 'vel [m/s]'))

                        tout = []
                        vout = []
                        velout = []
                        Tout = []
                        Pout = []
                        vtdc = []

                        #Take in the chemical kinetics mechanism
                        gas = ct.Solution('PRF_171.cti')
                        air = ct.Solution('air.xml')
                        gas.TPX = TC, PC*ct.one_atm,   'C2H5O2:{Mole}'.format(Mole=case)
                        print('Case Three :', TC, PC, case)
                        time.sleep(7)
                        r = ct.IdealGasReactor(gas, energy='on', volume = volstart) 
                        env = ct.Reservoir(contents=gas, name='environment')
                        wall = ct.Wall(r, env,  velocity=v_A)
                        sim = ct.ReactorNet([r])


                        for n in range(n_steps): 
                            if t <= taccel:
                                wall.set_velocity(v_A)
                                t += dt
                                sim.advance(t)
                                disp = vmax * (t**2/2.0/taccel)
                                v_A = -3.14/4.0 *(dia**2)* (vmax * t/taccel)
                                volume = volstart - (3.14/4.0 *(dia**2)*(disp))
                                tout.append(t)
                                Tout.append(r.T)
                                Pout.append(1.e-5*r.thermo.P)
                                vout.append(r.volume)
                                velout.append(1000*v_A)
                                csvfile = csv.writer(outfile)
                                csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, 1000*v_A])
                                print('%10.3f %10.3f %10.3f %10.3f %10.3f' % (t, r.T, 1.e-5*r.thermo.P, r.volume, v_A))
                                
                                
                            
                            elif t > taccel and t <= t2:
                                wall.set_velocity(v_C)
                                t += dt
                                sim.advance(t)
                                disp = vmax * taccel/2.0 + (vmax *(t-taccel))
                                v_C = -3.14/4.0 *((dia**2) * vmax)
                                volume = volstart - (3.14/4.0 *(dia**2)*(disp))
                                tout.append(t)
                                Tout.append(r.T)
                                Pout.append(1.e-5*r.thermo.P)
                                vout.append(r.volume)
                                velout.append(1000*v_C)
                                csvfile = csv.writer(outfile)
                                csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, 1000*v_C])
                                print('%10.3e %10.3f %10.3f %10.3f %10.3f' % (t, r.T, 1.e-5*r.thermo.P, r.volume, v_C))
                            
                            elif t > t2 and t <= tcomp:
                                wall.set_velocity(v_D)
                                t += dt
                                sim.advance(t)
                                disp = (vmax * taccel/2.0) + (vmax * tconst)
                                disp = (disp + vmax *(t-t2)) - (vmax *(t-t2)**2/tdecel/2.0)
                                v_D = -3.14/4.0 *(dia**2)*(vmax - vmax/tdecel * (t-t2))
                                volume = volstart - (3.14/4.0 *(dia**2)*(disp))
                                tout.append(t)
                                Tout.append(r.T)
                                Pout.append(1.e-5*r.thermo.P)
                                vout.append(r.volume)
                                velout.append(1000*v_D)
                                csvfile = csv.writer(outfile)
                                csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, 1000*v_D])
                                print('%10.3e %10.3f %10.3f %10.3f %10.3f' % (t, r.T, 1.e-5*r.thermo.P, r.volume, v_D))
                                vtdc = volume
                                
                                

                            elif  t > tcomp:
                                sim.set_initial_time(t)
                                ta = (t - tcomp)
                                poly = a[0] + a[1]*ta + a[2]*ta**2 + a[3]*ta**3 + a[4]*ta**4 + a[5]*ta**5 + a[6]*ta**6 + a[7]*ta**7 + a[8]*ta**8 + a[9]*ta**9 + a[10]*ta**10 
                                volume =  P.polymul(vtdc, poly)    
                                dv = a[1] + 2*a[2]*ta + 3*a[3]*ta**2 + 4*a[4]*ta**3 +  5*a[5]*ta**4 +  6*a[6]*ta**5 + 7*a[7]*ta**6 + 8*a[8]*ta**7 +  9*a[9]*ta**8 +  10*a[10]*ta**9    
                                v_E = P.polymul(vtdc, dv)
                                wall.set_velocity(v_E)
                                t += dt
                                sim.advance(t) 
                                tout.append(t)
                                Tout.append(r.T)
                                Pout.append(1.e-5*r.thermo.P)
                                vout.append(r.volume)
                                velout.append(v_E)
                                csvfile = csv.writer(outfile)
                                csvfile.writerow([t, r.T, 1.e-5*r.thermo.P, r.volume, v_E])
                                print('%10.3f %10.3f %10.3f %10.3f %10.3f' %(t, r.T, 1.e-5*r.thermo.P, r.volume, v_E))
                                

                        outfile.close()
                        print('Output written to file result_vel.csv')
                        print('Directory: '+os.getcwd())
                        print (vmax)          
                        print (volstart)
                        print(vtdc)


                        print('Directory: '+os.getcwd())
                        fig = plt.figure()
                        #plt.clf()
                        plt.subplot(2,2,1)
                        plt.plot(tout, vout)
                        plt.xlabel('Time (s)')
                        plt.ylabel('volume (m^3)')
                        plt.grid(False)

                        plt.subplot(2,2,2)
                        plt.plot(tout, velout)
                        plt.xlabel('Time (s)')
                        plt.ylabel('velocity (m/s)')
                        plt.grid(False)

                        plt.subplot(2,2,3)
                        plt.plot(tout, Tout)
                        plt.xlabel('Time (s)')
                        plt.ylabel('Temperature (K)')
                        plt.grid(False)

                        plt.subplot(2,2,4)
                        plt.plot(tout, Pout)
                        plt.xlabel('Time (s)')
                        plt.ylabel('Pressure (bar)')
                        plt.grid(False)
                        plt.tight_layout()
                        fig.savefig(str(TC) +'.'+str(Time)+'full_figure.png')
                        plt.show()