import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import curve_fit
import time
from astropy.table import Table
import matplotlib.pyplot as plt

# Creating a function to save time, displacement and velocity in a file
def save(chtable,tval,dval,vval):
  time.sleep(1.5)
  # Save Menu
  print(":::::::::::::::::::::::::  SAVE MENU  :::::::::::::::::::::::::",end='\n\n')

  print("NOTE: SAVE SCOPE",end='\n\n')
  print('****************************************************************************************************************************************',end='\n\n')
  print("> This menu is for saving the result table (time, displacement, velocity) only.")
  print("> Printed parameters above are not included in the saved file.",end='\n\n')
  print('****************************************************************************************************************************************',end='\n\n')
  time.sleep(1.5)

  print(">> Features Available:-")
  print("1. Save to a new file / Overwrite an existing file")
  print("2. Append results to an existing / original file")
  print("3. Do not save (return to main menu)",end='\n\n')
  print('****************************************************************************************************************************************',end='\n\n')
  time.sleep(1.5)

  print("NOTE: SAVE FORMAT",end='\n')
  print('****************************************************************************************************************************************',end='\n\n')
  print("Only .csv file format is supported for saving results in this program.",end='\n\n')
  print('****************************************************************************************************************************************',end='\n\n')
  time.sleep(1.5)

  while True:   # Loop continues till a correct choice is provided by the user
    svch=input("Enter a choice from the save menu: ").strip()
    print('\n')
    print('****************************************************************************************************************************************',end='\n\n')

    # If user wants to save the table in a new file or overwrite the original file
    if svch=='1':
      fname=input("Enter the name of a new file / existing file to save the data (in.csv): ")

      if not fname.lower().endswith('.csv'): # Ensuring .csv extension
        print("ONLY .csv format is supported for saving results in this program !!!")
        continue

      # Creating a new file and saving data / overwriting an existing file
      chtable.write(fname,format='csv',overwrite=True)

        # Ensuring save / overwrite was successful
      if Table.read(fname,format='csv')==chtable:
        print(f"The table has been SUCCESSFULLY saved in the file {fname} !!!")
        break
      else:  # If unsuccessful the save loop repeats
        print("The save was UNSUCCESSFUL!!! Please TRY AGAIN !!!")

    # If the user wants to append the columns to a table in a pre existing file
    elif svch=='2':
      fname=input("Enter the name of the file in which the columns are to be appended: ")

      if not fname.lower().endswith('.csv'):
        print("ONLY .csv format is supported for saving results in this program !!!")
        continue

      try:                                         # Ensuring that the file name provided exists
        owtable=Table.read(fname,format='csv')
      except FileNotFoundError:
        print(f"The file name {fname} provided does NOT EXIST !!! Please ENSURE that the file names provided DO EXIST !")
        continue

      if len(Table.read(fname,format='csv'))==0:   # Ensuring the table is non empty
        print(f"The table in the file {fname} provided is EMPTY !!! For such cases please select the FIRST choice !")
        continue

      tcol=input("Enter a name for the column consisting of the specific time stamps: ")
      xcol=input("Enter a name for the column consisting of displacement values: ")
      vcol=input("Enter a name for the column consisting of velocity values: ")
      print('\n')

      # Appending columns into the table
      owtable[tcol]=tval
      owtable[xcol]=dval
      owtable[vcol]=vval

      print(f"The updated table with the appended columns to be saved in the file {fname} is:",end='\n\n')
      print(owtable)

      owtable.write(fname,format='csv',overwrite=True)    # Overwriting the file with the updated table

      if Table.read(fname,format='csv')==owtable:
        print(f"The table has been SUCCESSFULLY saved in the file {fname} !!!")
        break
      else:
        print("The save was UNSUCCESSFUL!!! Please TRY AGAIN !!!")

    # If the user does not want to save the data
    elif svch=='3':
      print("Are you sure? The program will return back to the main menu and the data will NOT be SAVED !!!",end='\n\n')

      # Ensuring whether the user wants to really not save the data
      check=input(("If you want to continue with save, type (y / Y / yes / YES): "))
      if check in ['y','Y','yes','YES']: # If yes then save loop continues
        continue
      else:                              # Otherwise the save loop breaks
        break

    # If the user goes with an invalid choice from save menu
    else:
      print("INVALID choice !!! Please enter a VALID choice from Save Menu !")

###################################################################################################################################################################################

# Creating a function which prompts the user until a valid float or int is entered, based on dtype
def only_num(inpl,dtype):
  # Repeatedly ask until a valid int / float is entered
  while True:
    val=input(inpl)

    if dtype=='int':
      try:
        return int(val)     # Accept only integer input
      except ValueError:
        # Show error and retry
        print(f"{val} is NOT a VALID input !!!")
        print("Input involving ONLY INTEGER data type is ALLOWED !!!",end='\n\n')
        time.sleep(1.5)
        continue

    else:
      try:
        return float(val)   # Accept only float input
      except ValueError:
        print(f"{val} is NOT a VALID input !!!")
        print("Input involving ONLY FLOATING POINT data type is ALLOWED !!!",end='\n\n')
        time.sleep(1.5)
        continue

###################################################################################################################################################################################

# Creating a function to display a table in either summarized or full view based on user input
def display(display_string,table_to_display):
  # Table Display Menu
  print(":::::::::::::::::::::::::  TABLE DISPLAY MENU  :::::::::::::::::::::::::",end='\n\n')
  print("1. Summarized View  -  Default view showing only head and tail.")
  print("2. Full View        -  Entire table with all rows and columns.")
  time.sleep(1.5)

  dispch=input("Enter a valid choice from the table display menu: ").strip()
  print('\n')
  print('*************************************************************************************************************************************************',end='\n\n')
  time.sleep(1.5)

  if dispch!='2':
    print(display_string + '(in default view)',end='\n\n')
    print(table_to_display,end='\n\n')                                          # Display the default summary view (Astropy's standard table print)
    print('*************************************************************************************************************************************************',end='\n\n')
  else:
    print(display_string + '(in full view)',end='\n\n')
    table_to_display.pprint_all()                                               # Display the full table using Astropy's pprint_all() to show all rows and columns
    print('\n')
    print('*************************************************************************************************************************************************',end='\n\n')

###################################################################################################################################################################################

# Start of FDO
loop=0            # Creating a variable named loop
while True:       # This variable figures whether the loop about to begin is fresh or is  restarting
  if loop!=0:
    print("The program will now return back to MAIN MENU !!!",end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')
  else:
    print("Welcome to FDO (Forced Damped Oscillator) Simulator!",end='\n\n')

    print("Some IMPORTANT PRECAUTIONS before using FDO Simulator.",end='\n\n')
    time.sleep(0.5)

    print("1--NOTE: UNIT CONSISTENCY",end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')
    print("> This program does not enforce units. Ensure your inputs are consistent.")
    print("> All outputs reflect the units implied by user inputs - no automatic conversion is performed.",end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')
    time.sleep(1.5)

    print("2--NOTE: ASSUMPTION OF TRANSIENT DYNAMICS",end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')
    print("> All computations incorporate both transient and steady-state dynamics.")
    print("> Results may be unreliable if input data excludes the transient regime.",end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')

  time.sleep(1.5)
  # Main Menu
  loop+=1
  print(":::::::::::::::::::::::::  MAIN MENU  :::::::::::::::::::::::::",end='\n\n')
  print("1. Solve equation by providing physical parameters.")
  print("2. Estimate parameters from displacement-time data.")
  print("3. Visualize oscillator behavior graphically.")
  print("4. Terminate the session.",end='\n\n')
  time.sleep(1.5)

  # Asking user their choice from main menu
  ch=input("Please enter your choice from MAIN MENU: ").strip()
  print('\n')
  print('****************************************************************************************************************************************',end='\n\n')

###################################################################################################################################################################################
###################################################################################################################################################################################

  # If user goes with the choice of solving equations via physical parameters
  if ch=='1':
    print("NOTE: INPUT CONSTRAINTS FOR FDO SIMULATION",end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')
    print("> Must be non-negative: natural angular frequency of the system (ω₀), damping coefficient (β=b/2m), angular frequency of driving force (ω),")
    print("                        initial time (t₀), final time (tf).")
    print("> Can be positive or negative: initial displacement (x₀), initial velocity (v₀), initial acceleration imparted by the driving force (F₀/m).",end='\n\n')
    print("> Please make required sign changes to ensure the equation behaves correctly.",end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')
    time.sleep(1.5)

    # Asking user to provide the parameters so as to solve the FDO Ordinary Differential Equation
    print("Enter system parameters:",end='\n\n')
    time.sleep(0.5)

    # A check to raise error when a strictly positive variable points to either zero or a negative number
    try:
      errdisp=''   # Variable to store the quantity for which incorrect value was entered
      errval=0     # Variable to store the incorrect value which was entered by the user

      w0=B=w=t0=tf=0 # Setting all the physical parameter variables to zero
      varcstl=[w0,B,w,t0,tf]
      # List consisting of physical quantities to be displayed to the user upon an incorrect value
      errdispl=['ω₀','β','ω','t₀','tf']
      # User display list (displayed during value input)
      displayl=["Natural angular frequency of the system (ω₀): ","Damping coefficient (β=b/2m): ","Angular frequency of driving force (ω): ",
                "Enter the initial time stamp from where the evaluation of the system is to be initiated from (t₀): ",
                "Enter the final time stamp where the evaluation of the system is to be terminated (tf): "]

      for i in range (len(varcstl)):
        varcstl[i]=only_num(displayl[i],'float')
        print("\n")

        if varcstl[i]<0:
          errval=varcstl[i]
          errdisp=errdispl[i]
          raise ValueError

    except ValueError:
      print(f"ERROR: The VALUE of {errdisp} provided is {errval}, which is LESSER THAN 0 !!!")
      print("Please make required sign changes to ensure the equation behaves correctly!",end='\n\n')
      continue

    # Asking users to enter the value of parameters which can be either positive or negative
    a0 = only_num("Initial acceleration imparted by the driving force (F₀/m): ",'float')
    print("\n")
    x0 = only_num("Initial displacement (x₀): ",'float')
    print("\n")
    v0 = only_num("Initial velocity (v₀): ",'float')
    print("\n")
    print('****************************************************************************************************************************************',end='\n\n')


    # We define f to be the function consisting of the system of differential
    # equations that needs to be solved

    f=lambda t,z: [z[1],-2*B*(z[1])-(w0**2)*z[0]+(a0)*np.cos(w*t)]
    print('\n')

    # Ask the user how they want the time points in the solution:
    # If 'random', solve_ivp will choose adaptive internal time points (no t_eval given)
    # This simulates a natural numerical solution without predefined sampling.
    # If 'specific', the user provides exact time stamps via t_eval (e.g., matching real data)
    # This gives the user control over how the solution is sampled or compared.

    evalloop=0     # evalloop has the same function as the variable loop
    while True:    # Loop continues till the time the user inputs a valid choice
      if evalloop!=0:
        print("The program will now return back to the Evaluation Menu !")
        print('****************************************************************************************************************************************',end='\n\n')

      time.sleep(1.5)
      # Evaluation Menu
      evalloop+=1
      print(":::::::::::::::::::::::::  EVALUATION MENU  :::::::::::::::::::::::::",end='\n\n')
      print("--- Choose mode of evaluation---",end='\n\n')
      print("1 - Compute displacement at random stamps in the range [t₀, tf]")
      print("2 - Compute displacement at specific time points")
      print("3 - Exit from evaluation process",end='\n\n')
      time.sleep(1.5)

      # Asking user their choice from the evaluation menu
      chev=input("Enter a choice for mode of evaluation:").strip()
      print('\n')
      print('****************************************************************************************************************************************',end='\n\n')

###################################################################################################################################################################################

      # If the user wants to evaluate the system at random time stamps
      if chev=='1':
        # Solving the differential equation via solve_ivp() function
        sol=solve_ivp(f,(t0,tf),[x0,v0])

        # Creating a table to represent the data in the tabular form
        randomt=Table()
        randomt['Time Stamps']=sol.t
        randomt['Displacement']=sol.y[0]
        randomt['Velocity']=sol.y[1]

        # Printing the table
        display("The table consisting of random time stamps, displacement, and velocity at those time stamps is: ",randomt)

###################################################################################################################################################################################

      # If the user wants to evaluate the system at specific time stamps
      elif chev=='2':
        print(":::::::::::::::::::::::::  TIME SOURCE MENU  :::::::::::::::::::::::::",end='\n\n')
        print("1. Evaluation at specific time stamps from a CSV file")
        print("2. Enter specific time stamps manually",end='\n\n')
        time.sleep(1.5)

        # Asking user for their choice from time source menu
        chsp=input("Enter a choice from the Time Source Menu: ").strip()
        print('\n')
        print('****************************************************************************************************************************************',end='\n\n')

        tl=[]  # List containing the specific time stamps

        # If the user wants to extract the time stamps from a .csv file
        if chsp=='1':
          namef=input("Enter the name of the file consisting of the specific time intervals at which evaluation is to be done (in .csv): ")

          if not namef.lower().endswith('.csv'):          # Ensuring .csv extension
            print("ONLY .csv format is supported for saving results in this program !!!")
            continue

          try:                                            # Ensuring the file name being provided exists
            evalt=Table(namef,format='csv')
          except FileNotFoundError:
            print(f"The file name {namef} provided does NOT EXIST !!! Please ENSURE that the file names provided DO EXIST !")
            continue

          if len(Table.read(namef,format='csv'))==0:      # Ensuring non empty table
            print(f"The table in the file {namef} provided is EMPTY !!! ALWAYS ensure the file is NON-EMPTY !")
            continue

          print(f"The table extracted from the file {namef} is:",end='\n\n')
          print(evalt,end='\n\n')
          print('****************************************************************************************************************************************',end='\n\n')

          # Asking user the name of the column consiting of the time stamps
          tcolname=input("Enter the name of the column consisting of the specific time stamps for evaluation: ")

          # Ensuring the column name provided by the user exist and have a numeric data type
          if tcolname in evalt.colnames and isinstance(evalt[tcolname][0],(float,np.floating)):
            tl=list(evalt[tcolname])

          else: # If the column name is not present or the values are not of numeric data type
            if tcolname not in evalt.colnames:
              print(f"The column name {tcolname} does NOT EXIST in the file {namef} provided !!! ")
              continue               # The loop repeats and the user is sent back to the evaluation menu
            else:
              print(f"The column name {tcolname} provided does NOT contain numeric data type !!!")
              continue

        # If the user wants to manually enter the time stamps
        elif chsp=='2':
          n=only_num("Enter the no.of time stamps for which the system has to be evaluated for:",'int')
          print('\n')
          for i in range(n):
            val=only_num(f"Enter time stamp {i+1} : ",'float')
            tl.append(val)
            print('\n')

        # If the user selects an incorrect choice from the time source menu
        else:
          print("INVALID choice from TIME SOURCE MENU !!!")
          continue        # The loop repeats and goes back to the evaluation menu

        sol=solve_ivp(f,(t0,tf),[x0,v0],t_eval=tl)

        specifict=Table()
        specifict['Specific Time Stamps']=sol.t
        specifict['Displacement']=sol.y[0]
        specifict['Velocity']=sol.y[1]

        display("The table consisting of specific time stamps, displacement, and velocity at those time stamps is: ",specifict)

###################################################################################################################################################################################

      # If the user wants to return back to the main menu
      elif chev=='3':
        print("The program will return back to Main Menu.",end='\n\n')
        # Ensuring user really wants to terminate the evluation process
        sure=input("If you want to still continue with evaluation process, enter (y / Y / yes / YES): ")
        if sure in ['y','Y','yes','YES']:
          continue                                                                                       # If yes, user returns to evaluation menu
        else:
          break
                                                                                                         # Else, user returns to main menu
      if chev=='3':
        continue

###################################################################################################################################################################################

      # If the user goes with an invalid choice from the evaluation menu
      else:
        print("INVALID choice !!! Please enter a valid choice from Evaluation Menu !")
        continue

###################################################################################################################################################################################

      # POST-COMPUTATION INTERFACE
      # ---------------------------------------------------------------

      # Calculating the amplitude of the displacement that would be obtained under steady state behaviour
      A=(a0)/((((w0**2)-w**2)**2+(2*B*w)**2)**(0.5))

      # Calculating the phase difference δ
      delta=np.arctan2((2*B*w),((w0**2)-(w**2)))

      # Calculating the resonance angular frequency
      wres= ((w0**2)-(2*(B**2)))**(0.5)

      # Calculating the quality factor
      Q=w0/(2*B)

      # Printing the list of derived parameters
      print("LIST OF DERIVED PARAMETERS:-",end='\n\n')
      print("1) Phase difference (δ): ",delta,end='\n\n')
      print("2) Amplitude of steady-state oscillation (A): ",A,end='\n\n')
      print("3) Natural angular frequency (ω₀): ",w0,end='\n\n')
      print("4) Resonance angular frequency (ω_res): ",wres,end='\n\n')
      print("5) Quality factor (Q): ",Q,end='\n\n')
      print("6) Maximum velocity during evaluation: ",A*w,end='\n\n')
      print('****************************************************************************************************************************************',end='\n\n')

      # Saving the result table
      if chev=='1':
        save(randomt,sol.t,sol.y[0],sol.y[1])

      elif chev=='2':
        save(specifict,sol.t,sol.y[0],sol.y[1])

      else:
        continue

###################################################################################################################################################################################
###################################################################################################################################################################################

  # If the user goes with the choice of estimating parameters from displacement-time data
  elif ch=='2':
    print(":::::::::::::::::::::::::  DATA SOURCE MENU  :::::::::::::::::::::::::",end='\n\n')
    print("1. Load time-displacement data from a CSV file")
    print("2. Enter time-displacement data manually",end='\n\n')
    time.sleep(1.5)

    # Asking user for their choice from data source menu
    chds=input("Enter a choice from the Data Source Menu: ").strip()
    print('\n')
    print('****************************************************************************************************************************************',end='\n\n')

    # If the user goes with the choice of loading displacement-time data from a .csv file
    if chds=='1':
      # Asking user for the name of the file consisting of the data sets regarding the time stamps and displacement
      fn=input("Please enter the name of the file containing the time and displacement data for the system (in.csv): ")
      print('\n')

      if not fn.lower().endswith('.csv'):
        print("ONLY .csv format is supported for saving results in this program !!!")
        continue

      try:
        parat=Table.read(fn,format='csv')
      except FileNotFoundError:
        print(f"The file name {fn} provided does NOT EXIST !!! Please ENSURE that the file names provided DO EXIST !")
        continue

      if len(Table.read(fn,format='csv'))==0:
        print(f"The table in the file {fn} provided is EMPTY !!! ALWAYS ensure the file is NON-EMPTY !")
        continue

      # To ensure that the no.of data points is never less than 3 - ensuring reliable parameter estimation
      if len(Table.read(fn,format='csv'))<3:
        print("ERROR: AT LEAST 3 time-displacement points ARE REQUIRED for parameter estimation !!!")
        print(f"The no.of time-displacement points were found to be {len(Table.read(fn,format='csv'))} which is less than 3 !")
        print("ALWAYS ensure that the no.of data points are GREATER THAN OR EQUAL to 3 for reliable paramter estimation ! ")
        continue

      print(f"The table extracted from the file {fn} is:",end='\n\n')
      print(parat,end='\n\n')
      print('****************************************************************************************************************************************',end='\n\n')

      x=input("Enter the name of the column consisting of the data regarding system's displacement: ")
      print('\n')
      t=input("Enter the name of the column consisting of the data regarding the time at which system's displacement was recorded: ")
      print('\n')

      if x in parat.colnames and t in parat.colnames and isinstance(parat[x][0],(float,np.floating)) and isinstance(parat[t][0],(float,np.floating)):
        tl,xl=parat[t],parat[x]

      else:
        if x not in parat.colnames:
          print(f"The column name {x} does NOT EXIST in the file {fn} provided !!!")
          continue
        elif t not in parat.colnames:
          print(f"The column name {t} does NOT EXIST in the file {fn} provided !!!")
          continue
        elif parat[x][0] is not isinstance(parat[x][0],(float,np.floating)):
          print(f"The column name {x} provided DOES NOT contain numeric data type !!!")
          continue
        else:
          print(f"The column name {t} provided DOES NOT contain numeric data type !!!")
          continue

    # If the user goes with the choice of manually inserting the displacement time data
    elif chds=='2':
      tli,xli=[],[]    # Creating empty lists tli and xli which would store the time and displacement data respectively

      dp=only_num("Enter the no.of data points available for paramter estimation: ",'int')

      if dp<3:
        print("ERROR: AT LEAST 3 time-displacement points ARE REQUIRED for parameter estimation !!!")
        print(f"The no.of time-displacement points were found to be {dp} which is less than 3 !")
        print("ALWAYS ensure that the no.of data points are GREATER THAN OR EQUAL to 3 for reliable paramter estimation ! ")
        continue

      for i in range(dp):
        tdp=only_num(f"Enter time stamp value {i+1}: ",'float')
        xdp=only_num(f"Enter displacemnt at {tdp} s in (m/s): ",'float')
        tli.append(tdp)
        xli.append(xdp)

      # Craeting numpy arrays of the lists tli and xli
      tl=np.array(tli,dtype=float)
      xl=np.array(xli,dtype=float)

    else:
      print("INVALID choice from Data Source Menu !!!")
      continue


    # The equation providing the relationship between the displacement and time is given as-
    fun=lambda t,A,w,d,C,B,wd: (A*(np.cos((w*t)-d)))+(C*(np.exp(-B*t))*np.cos((wd)*t))

    # Finding the suitable paramters for the equation provided above which best suits the data set provided
    para,pcov=curve_fit(fun,tl,xl)

    amp=para[0]   # Amplitude of the system during steady state behaviour
    afdf=para[1]  # Angular frequency of the driving force
    delta=para[2] # Phase difference δ
    C_=para[3]    # Amplitude of the system during the transient behaviour
    B_=para[4]    # Damping parameter
    afdp=para[5]  # Damped angular frequency

    # Calculating the natural angular frequency of the system
    afnt=(((afdp)**2)+((B_)**2))**(0.5)

    # Calculating the initial acceleration of the system provided by the driving force
    ai=(amp)*(((((afnt**2)-(afdf**2))**2)+((afdf*(B_)*2)**2))**0.5)

    # Calculating the resonance angular frequency
    wres= ((afnt**2)-(2*((B_)**2)))**(0.5)

    # Calculating the quality factor
    if B_!=0:
      Q=afnt/(2*(B_))
    else:
      Q="UNDEFINED (damping parameter is zero)"

    time.sleep(1.5)
    # Printing the time-dependent equation of motion and the estimated physical paramters
    print('****************************************************************************************************************************************',end='\n\n')

    print("Given below is the time-dependent equation of motion (with transient effects), based on estimated parameters from the provided data:",end='\n\n')
    print(f"x(t) = ({amp} · cos({afdf} · t - {delta})) + {C_} · (e^(-{B_} · t)) · cos({afdp} · t)",end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')
    time.sleep(1.5)

    print("The natural angular frequency of the system is: ",afnt,end='\n\n')
    print("The damping coefficient ratio (b/m) of the system is: ",2*(B_),end='\n\n')
    print("Resonance angular frequency ω_res: ",wres,end='\n\n')
    print("Quality factor Q: ",Q,end='\n\n')
    print("The initial acceleration of the system (F0/m) provided by the driving force is: ",ai,end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')
    time.sleep(1.5)

    print("Note: These values are determined WITHOUT knowing the MASS of the system !")
    print('****************************************************************************************************************************************',end='\n\n')
    print("To find absolute quantities, please multiply the ratio values by the actual mass (m) to obtain physical quantities.",end='\n\n')
    print("1) Spring constant: k=m*(ω0^2), where ω₀ is the natural angular frequency of the system.")
    print("2) Damping coefficient: b=m*(b/m)")
    print("3) Driving force amplitude: F0=m*(F0/m)",end='\n\n')
    print('****************************************************************************************************************************************',end='\n\n')
    time.sleep(1.5)

###################################################################################################################################################################################
###################################################################################################################################################################################

  # If the user goes with the choice of visualizing the oscillator behaviour graphically
  elif ch=='3':

    # Creating a function which asks user for figure metadeta
    def common(plot_type):
      global xlabel,ylabel,ptitle,label,lblcheck
      xlabel,ylabel,ptitle='','',''
      label=None
      lblcheck=0

      print(f"NOTE: The following entries apply to → {plot_type} plot.", end='\n\n')
      time.sleep(1)

      xlabel=input("Enter a label for X-axis: ")
      ylabel=input("Enter a label for Y-axis: ")
      ptitle=input("Enter the title for the plot (leave blank for none): ")

      asklgnd=input("Enter (y / Y / yes / YES) to include legend in the plot: ")
      if asklgnd in ['y','Y','yes','YES']:
        lblcheck=1
        label=input("Enter label for the curve (to appear in legend): ")

###################################################################################################################################################################################

    # Creating a function which displays and saves the plot
    def showsave(num):
      # Displaying the graph
      if num!=1:
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(1)
        print('\n')
        print('****************************************************************************************************************************************',end='\n\n')

      # Creating a loop with a two-tier verification, asking user whether to save the plot or not
      while True:
        # Asking user whether the plot is to be save
        pltsave=input("Enter (y / Y / yes / YES), if you want to save the above graph: ")
        # If the user wants to save the plot
        if pltsave in ['y','Y','yes','YES']:
          print("SETTINGS: Save Settings – ACADEMIC DEFAULTS APPLIED")
          print("DPI = 300 | BBOX_INCHES = 'tight' | TRANSPARENT = False")
          print('****************************************************************************************************************************************',end='\n\n')
          time.sleep(1.5)

          print("NOTE: Only .png and .pdf formats are supported. Please name your file accordingly.", end='\n\n')
          print('****************************************************************************************************************************************',end='\n\n')
          time.sleep(1.5)

          while True:
            # Asking user for the name of the file where the plot will be saved
            pltf=input("Enter a name for the file where the plot is to be saved (include .png or .pdf extension): ")

            # Ensuring that the file has only .png or .pdf extension
            if not (pltf.lower().endswith('.png') or pltf.lower().endswith('.pdf')):
              print(f"The file name {pltf} does NOT end with a .png or .pdf extension !!!")
              time.sleep(0.5)
              print("Please try again !!!",end='\n\n')
              time.sleep(0.5)
              continue
            else:
              break
          # Saving the plot
          plt.savefig(pltf,dpi=300,bbox_inches='tight',transparent=False)
          print("Plot has been SAVED SUCCESSFULLY !!!",end='\n\n')

        # If the use does not want to save the plot
        else:
          print("Are you sure you don't want to save the plot? The plot would be discarded and shall not be retrievable.")
          pltsure=input("If you don't want to loose progress, enter (y / Y / yes / YES): ")   # Reaffirming user's decision
          # If the user wants to continue with saving the plot
          if pltsure in ['y','Y','yes','YES']:
            time.sleep(1)
            continue    # The loop starts again
          # If the user wants to discard the plot
          else:
            break       # The loop breaks

###################################################################################################################################################################################

    vizloop=0
    while True:
      if vizloop!=0:
        print("The program is now returning back to the VISUALIZATION MENU !!!",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
      else:
        print("NOTE: ACADEMIC DEFAULTS APPLIED",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
        print("> Certain VISUAL SETTINGS are FIXED for CLARITY and NOT USER-CONFIGURABLE.")
        print("> It will be SHOWN when a PLOT is SELECTED from the MENU.")
        print("> To MODIFY them, EDIT the SOURCE CODE DIRECTLY.",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')

      time.sleep(1.5)
      # Visualization Menu
      vizloop+=1
      print(":::::::::::::::::::::::::  VISUALIZATION MENU  :::::::::::::::::::::::::",end='\n\n')
      print("1. Dynamic Response (x(t), v(t), a(t))")
      print("2. Amplitude Response A(ω)")
      print("3. Phase Lag Response δ(ω)")
      print("4. Exit Visualization Menu",end='\n\n')
      time.sleep(1.5)

      # Asking user for their choice from the visualization menu
      vizch=input("Enter a choice from Visualization Menu: ").strip()
      print('\n')
      print('****************************************************************************************************************************************',end='\n\n')

      # Displaying the default settings to the user
      if vizch=='1' or vizch=='2' or vizch=='3':
        print("SETTINGS: Line Plot – ACADEMIC DEFAULTS APPLIED", end='\n\n')
        print("COLOR = 'black' | LINESTYLE = '-' | MARKER = OFF | GRID = ON", end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
        time.sleep(1.5)

###################################################################################################################################################################################

      # If the user wants to visualize dynamic response
      if vizch=='1':
        # Display input constraints and domain-specific information
        print("NOTE: INPUT CONSTRAINTS FOR FDO VISUALIZATION",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
        print("> Must be non-negative: steady state amplitude (A), angular frequency of driving force (ω), damping coefficient of the system (β=b/2m),")
        print("                        damped angular frequency of the system (ωd), phase lag (δ), initial time (t₀), final time (tf).")
        print("> Can be positive or negative: transient state amplitude (C).",end='\n\n')
        print("> Please make required sign changes to ensure the equation behaves correctly.",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
        time.sleep(1.5)

        print("Enter system parameters:",end='\n\n')
        time.sleep(0.5)

        try:
          errdisp=''
          errval=0

          A=w=B=wd=t0=tf=0
          varcstl=[A,w,wd,B,t0,tf]

          errdispl=['A','ω','ωd','β','t₀','tf']

          displayl=["Steady state amplitude of the system (A): ","Angular frequency of driving force (ω): ","Damped angular frequency of the system (ωd): ",
                    "Damping coefficient of the system (β=b/2m): ","Enter the initial time stamp from where the visualization of the system is to be initiated from (t₀): ",
                    "Enter the final time stamp where the visualization of the system is to be terminated (tf): "]

          for i in range (len(varcstl)):
            varcstl[i]=only_num(displayl[i],'float')
            print("\n")

            if varcstl[i]<0:
              errval=varcstl[i]
              errdisp=errdispl[i]
              raise ValueError

        except ValueError:
          print(f"ERROR: The VALUE of {errdisp} provided is {errval}, which is LESSER THAN 0 !!!")
          print("Please make required sign changes to ensure the equation behaves correctly !",end='\n\n')
          continue

        # Ensuring that the no.of evaluation point is always ≥ 2, to plot the graph
        while True:
          try:
            steps=only_num("Enter number of evaluation points between t_start and t_end (≥ 2): ",'int')
            if steps<2:
              raise ValueError
            else:
              break
          except ValueError:
            print(f"ERROR: INVALID evaluation point count: {steps} !!! ATLEAST 2 steps ARE required !!!")
            continue

        # Enusring that the phase lag is always in the range of [0,π]
        while True:
          try:
            plag=only_num("Phase lag (δ) (in radians) [0 ≤ δ ≤ π]: ",'float')
            if not 0<=plag<=(np.pi):
              raise ValueError
            else:
              break
          except ValueError:
            print(f"Please ENTER CORRECT values !!! The δ provided is OUT of [0,π] range !!!")
            continue

        c=only_num("Transient state amplitude of the system (C): ",'float')

        # Create a uniformly spaced time array to evaluate the system's motion from t₀ to t_f
        t=np.linspace(t0,tf,steps)
        # Time-dependent equation of motion generated from the physical paramterers provided by the user
        eqn=(A*(np.cos((w*t)-plag)))+(c*(np.exp(-B*t))*np.cos((wd)*t))

        dfdx=np.gradient(eqn,t)                                                 # First-order time derivative of displacement (velocity)
        d2fdx2=np.gradient(dfdx,t)                                              # Second-order time derivative of displacement (acceleration)


        print("The time-dependent equation of motion based upon the paramters received is:",end='\n\n')
        print(f"x(t) = ({A} · cos({w} · t - {plag})) + {c} · (e^(-{B} · t)) · cos({wd} · t)",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
        time.sleep(1.5)

        # Create vertically stacked subplots for displacement, velocity, and acceleration
        # Sharing the same time axis (x-axis) across all three plots
        fig,ax=plt.subplots(3,1,figsize=(4,9),sharex=True)

        # Defining titles for each subplot
        comdisp=['Displacement-Time','Velocity-Time','Acceleration-Time']
        # Store the corresponding y-data for each plot
        yplot=[eqn,dfdx,d2fdx2]
        #Creating lists to store metadeta for each plot (labels and titles)
        xll,yll,ptl,lbl=[],[],[],[]

        # Generating plots one by one using common() to get user input for metadata
        for i in range(3):
          common(comdisp[i])

          # Checking if legend is requested by user (set via global 'lblcheck' in common())
          # With legend
          if lblcheck==1:
            ax[i].plot(t,yplot[i],label=label,color='black',linestyle='-')
            lbl.append(label)
          # Without legend
          else:
            ax[i].plot(t,yplot[i],color='black',linestyle='-')
            lbl.append('')

          # Store user input for later reuse in the save block
          xll.append(xlabel)
          yll.append(ylabel)
          ptl.append(ptitle)

        # Apply labels, titles, and legend/grid settings for all 3 axes
        for index,axis in enumerate(ax):
          axis.set_title(ptl[index])
          axis.set_xlabel(xll[index])
          axis.set_ylabel(yll[index])
          if lbl[index]!='':
            axis.legend(True)
          axis.grid(True)

        # Format layout, display the combined subplot figure
        plt.tight_layout()
        plt.show(block=False)
        plt.pause(1)
        print('\n')
        print('****************************************************************************************************************************************',end='\n\n')

        # Asking user whether to save each subplot separately or save the full figure
        svind=input("Enter (y / Y / yes / YES) to save each subplot (x, v, a) as a separate image: ")

        # If the user wants to save individual plots: one figure per curve
        if svind in ['y','Y','yes','YES']:
          for k in range(3):
            if lbl[k]!='':
              plt.plot(t,yplot[k],label=lbl[k],color='black',linestyle='-')
              plt.legend()
            else:
              plt.plot(t,yplot[k],color='black',linestyle='-')

            plt.title(ptl[k])
            plt.xlabel(xll[k])
            plt.ylabel(yll[k])
            plt.grid(True)
            showsave(0)        # Save individual figure and close

        # If the user wants to save the full 3-subplot figure as a single file
        else:
          fig,ax=plt.subplots(3,1,figsize=(4,9),sharex=True)
          for k in range(3):
            if lbl[k]!='':
              ax[k].plot(t,yplot[k],label=lbl[k],color='black',linestyle='-')
            else:
              ax[k].plot(t,yplot[k],color='black',linestyle='-')
          for index,axis in enumerate(ax):
            axis.set_title(ptl[index])
            axis.set_xlabel(xll[index])
            axis.set_ylabel(yll[index])
            if lbl[index]!='':
              axis.legend(True)
            axis.grid(True)
          showsave(1)          # Save the full combined figure and close

###################################################################################################################################################################################

      # If the user wants to visualize amplitude (vizch=='2') or phase lag (vizch=='3') response
      elif vizch=='2' or vizch=='3':
        print("NOTE: INPUT CONSTRAINTS FOR FDO VISUALIZATION",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
        print("> Must be non-negative: natural angular frequency of the system (ω₀), angular frequency of driving force (ω), damping coefficient of the system (β=b/2m)")
        print("> Driving angular frequency (ω) values MUST be ≥ 0, with ω_start < ω_end, and number of evaluation points ≥ 2.")

        if vizch=='2':
          print("> Can be positive or negative: initial acceleration imparted by the driving force (F₀/m).",end='\n\n')

        print("> Please make required sign changes to ensure the equation behaves correctly.",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
        time.sleep(1.5)

        # Informing the user regarding x-axis display flexibility (ω vs f)
        print("NOTE: X-AXIS DISPLAY CONTROL",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
        print("> This section uses angular frequency ω (rad/s) for all internal calculations.")
        print("> You may still choose to *display* results in linear frequency f (Hz) if preferred.")
        print("> Simply provide ω_start and ω_end as inputs when asked, and visualization units will adjust accordingly.",end='\n\n')
        print('****************************************************************************************************************************************',end='\n\n')
        time.sleep(1.5)

        while True:
          try:
            B=only_num("Damping coefficient of the system (β=b/2m): ",'float')
            if B<0:
              raise ValueError
            else:
              break
          except ValueError:
            print(f"ERROR: The VALUE of β provided is {B}, which is LESSER THAN 0 !!!")
            print("Please make required sign changes to ensure the equation behaves correctly !",end='\n\n')
            continue

        while True:
          try:
            w0=only_num("Natural angular frequency of the system (ω₀) (in rad/s): ",'float')
            if w0<0:
              raise ValueError
            else:
              break
          except ValueError:
            print(f"ERROR: The VALUE of ω₀ provided is {w0}, which is LESSER THAN 0 !!!")
            print("Please make required sign changes to ensure the equation behaves correctly !",end='\n\n')
            continue

        # Asking user for initial acceleration only for amplitude plot (vizch == 2)
        if vizch=='2':
          a0=only_num("Initial acceleration imparted by the driving force (F₀/m): ",'float')

        # Prompt for angular frequency range
        print("Enter ω (in rad/s) range:",end='\n\n')
        while True:
          time.sleep(0.5)

          try:
            ws=only_num("Enter starting value for the ω range: ",'float')
            if ws<0:
              raise ValueError
          except ValueError:
            print(f"Starting value of ω range was provided as {ws} which is <0 !!! ω MUST ALWAYS be ≥ 0 !!!")
            print("Please try again !!!")
            continue

          try:
            we=only_num("Enter ending value for the ω range: ",'float')
            if we<0:
              raise ValueError
          except ValueError:
            print(f"Ending value of ω range was provided as {we} which is <0 !!! ω MUST ALWAYS be ≥ 0 !!!")
            print("Please try again !!!")
            continue

          try:
            if not ws<we:
              raise ValueError
          except ValueError:
            print("The starting value of ω range CANNOT be LESS THAN ending value of ω range !!!")
            print("Please try again !!!")
            continue

          break

        while True:
          try:
            steps=only_num("Enter number of evaluation points (≥ 2): ",'int')
            if steps<2:
              raise ValueError
            else:
              break
          except ValueError:
            print(f"ERROR: INVALID evaluation point count: {steps} !!! ATLEAST 2 steps ARE required !!!")
            continue

        # Generate angular frequency array and compute equivalent linear frequency array
        wrange=np.linspace(ws,we,steps)
        frange=(wrange)/(2*np.pi)

        # Asking user whether to display linear frequency (f) (Hz) instead of angular frequency (ω) (rad/s) on x-axis
        hzx=input("Enter (y / Y / yes / YES), to display linear frequency (f) (Hz) instead of angular frequency (ω) (rad/s) on x-axis: ")

        ##################################################################################################################################

        # If the user wants to visualize amplitude response
        if vizch=='2':
          # Calculating amplitude response A(ω)
          A=(a0)/((((w0**2)-(wrange)**2)**2+(2*B*(wrange))**2)**(0.5))
          common('Amplitude Response A(ω)')

          # Plot with x-axis in Hz
          if hzx in ['y','Y','yes','YES']:
            if lblcheck==1:
              plt.plot(frange,A,label=label,color='black',linestyle='-')
              plt.legend()
            else:
              plt.plot(frange,A,color='black',linestyle='-')

          # Plot with x-axis in rad/s
          else:
            if lblcheck==1:
              plt.plot(wrange,A,label=label,color='black',linestyle='-')
              plt.legend()
            else:
              plt.plot(wrange,A,color='black',linestyle='-')

        ##################################################################################################################################

        # If the user wants to visualize phase lag response
        else:
          # Calculating phase lag response δ(ω)
          plag=np.arctan2((2*B*(wrange)),((w0**2)-((wrange)**2)))
          common('Phase Lag Response δ(ω)')

          if hzx in ['y','Y','yes','YES']:
            if lblcheck==1:
              plt.plot(frange,plag,label=label,color='black',linestyle='-')
              plt.legend()
            else:
              plt.plot(frange,plag,color='black',linestyle='-')

          else:
            if lblcheck==1:
              plt.plot(wrange,plag,label=label,color='black',linestyle='-')
              plt.legend()
            else:
              plt.plot(wrange,plag,color='black',linestyle='-')

        ##################################################################################################################################

        # Applying the plot metadata
        plt.title(ptitle)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.grid(True)

        # Saving the final plot
        showsave(0)

###################################################################################################################################################################################

      # If the user wants to return back to the main menu
      elif vizch=='4':
        print("The program will return back to Main Menu.",end='\n\n')
        # Ensuring user really wants to terminate the visualization process
        vizsure=input("If you want to still continue with visualization process, enter (y / Y / yes / YES): ")

        if vizsure in ['y','Y','yes','YES']:                                    # If yes, user returns back to visualization menu
          continue
        else:                                                                   # Else, user returns back to main menu
          break

###################################################################################################################################################################################

      # If the user goes with an invalid choice from the visualization menu
      else:
        print("INVALID choice from Visualization Menu !!! Please enter a VALID choice !")

###################################################################################################################################################################################
###################################################################################################################################################################################

  # If the user goes with the choice of terminating the session
  elif ch=='4':
    print("Are you sure? Your session will be terminated and you would exit the simulator.")

    # Ensuring that the user really wants to terminate the current session
    mainsure=input("If you still want to continue with this session, please enter (y / Y / yes / YES): ")
    if mainsure in ['y','Y','yes','YES']:    # If the user still wants to continue with the session, the user is returned back to the main menu
      continue
    else:                                    # If the user wants to terminate the session
      print("Thank you for using the FDO (Forced Damped Oscillator) Simulator !!")
      print("May your oscillations decay gracefully, and your resonance never be destructive !! :)",end='\n\n')
      time.sleep(1)
      print("Terminating your session...")
      time.sleep(1.5)
      print(".....",end='\n\n')
      time.sleep(1.5)
      print("Session termianted SUCCESSFULLY ! This system has now returned to equilibrium !")
      break

###################################################################################################################################################################################
###################################################################################################################################################################################

  # If the user goes with an invalid choice from the main menu
  else:
    print("INVALID choice from Main Menu !!! Please enter a VALID choice !")

###################################################################################################################################################################################
###################################################################################################################################################################################

