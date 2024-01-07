import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots


def extract_mpms_data(df, specimen_name):
    """
    Extracts and separates MPMS (Magnetic Property Measurement System) data 
    for a specific specimen from a dataframe.

    This function filters data for a given specimen and separates it based on 
    different MagIC measurement method codes. It specifically looks for data 
    corresponding to 'LP-FC' (Field Cooled), 'LP-ZFC' (Zero Field Cooled),
    'LP-CW-SIRM:LP-MC' (Room Temperature SIRM measured upon cooling), and 
    'LP-CW-SIRM:LP-MW' (Room Temperature SIRM measured upon Warming).

    Parameters:
        df (pandas.DataFrame): The dataframe containing MPMS measurement data.
        specimen_name (str): The name of the specimen to filter data for.

    Returns:
        tuple: A tuple containing four pandas.DataFrames:
            - fc_data: Data filtered for 'LP-FC' method.
            - zfc_data: Data filtered for 'LP-ZFC' method.
            - rtsirm_cool_data: Data filtered for 'LP-CW-SIRM:LP-MC' method.
            - rtsirm_warm_data: Data filtered for 'LP-CW-SIRM:LP-MW' method.

    Example:
        >>> measurements_df = pd.DataFrame({...})
        >>> fc, zfc, rtsirm_cool, rtsirm_warm = extract_mpms_data(measurements_df, 'Specimen_1')
    """

    specimen_df = df[df['specimen'] == specimen_name]

    fc_data = specimen_df[specimen_df['method_codes'].str.contains('LP-FC')]
    zfc_data = specimen_df[specimen_df['method_codes'].str.contains('LP-ZFC')]
    rtsirm_cool_data = specimen_df[specimen_df['method_codes'].str.contains('LP-CW-SIRM:LP-MC')]
    rtsirm_warm_data = specimen_df[specimen_df['method_codes'].str.contains('LP-CW-SIRM:LP-MW')]

    return fc_data, zfc_data, rtsirm_cool_data, rtsirm_warm_data


def get_plotly_marker(matplotlib_marker):
    """
    Maps a Matplotlib marker style to its Plotly equivalent.

    Parameters:
        matplotlib_marker (str): A marker style in Matplotlib format.

    Returns:
        str: Corresponding marker style in Plotly format.
    """
    marker_dict = {
        '.': 'circle',
        'o': 'circle',
        '^': 'triangle-up',
        's': 'square',
        '*': 'star',
        '+': 'cross',
        'x': 'x',
        'D': 'diamond',
        '|': 'line-ns',
        '_': 'line-ew',
    }
    return marker_dict.get(matplotlib_marker, 'circle')  # Default to 'circle' if not found


def plot_mpms_data(fc_data, zfc_data, rtsirm_cool_data, rtsirm_warm_data, 
                   fc_color='#1f77b4', zfc_color='#ff7f0e', rtsirm_cool_color='#2ca02c', rtsirm_warm_color='#d62728',
                   fc_marker='.', zfc_marker='.', rtsirm_cool_marker='.', rtsirm_warm_marker='.',
                   symbol_size=10, use_plotly=False, plot_derivative=False):
    """
    Plots MPMS data and optionally its derivatives for Field Cooled, Zero Field Cooled, RTSIRM Cooling, and RTSIRM Warming using either Matplotlib or Plotly.

    Parameters:
        fc_data, zfc_data, rtsirm_cool_data, rtsirm_warm_data (DataFrame): DataFrames containing the MPMS data.
        fc_color, zfc_color, rtsirm_cool_color, rtsirm_warm_color (str): HEX color codes for each plot.
        fc_marker, zfc_marker, rtsirm_cool_marker, rtsirm_warm_marker (str): Marker symbols for each plot.
        symbol_size (int): Size of the markers in matplotlib, symbol size in plotly is fixed.
        use_plotly (bool): If True, uses Plotly for plotting. Otherwise, uses Matplotlib.
        plot_derivative (bool): If True, plots the derivative of the magnetization data.
    """
    if plot_derivative:
        fc_derivative = thermomag_derivative(fc_data['meas_temp'], fc_data['magn_mass'])
        zfc_derivative = thermomag_derivative(zfc_data['meas_temp'], zfc_data['magn_mass'])
        rtsirm_cool_derivative = thermomag_derivative(rtsirm_cool_data['meas_temp'], rtsirm_cool_data['magn_mass'])
        rtsirm_warm_derivative = thermomag_derivative(rtsirm_warm_data['meas_temp'], rtsirm_warm_data['magn_mass'])


    if use_plotly:
        # Create subplot layout
        fig = make_subplots(rows=1, cols=2, subplot_titles=("FC and ZFC Data", "RTSIRM Cooling and Warming Data"))

        # Adding FC and ZFC data to the first subplot
        fig.add_trace(go.Scatter(x=fc_data['meas_temp'], y=fc_data['magn_mass'], mode='markers+lines', 
                                 name='FC', marker=dict(color=fc_color, symbol=get_plotly_marker(fc_marker))),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=zfc_data['meas_temp'], y=zfc_data['magn_mass'], mode='markers+lines', 
                                 name='ZFC', marker=dict(color=zfc_color, symbol=get_plotly_marker(zfc_marker))),
                      row=1, col=1)

        # Adding RTSIRM Cooling and Warming data to the second subplot
        fig.add_trace(go.Scatter(x=rtsirm_cool_data['meas_temp'], y=rtsirm_cool_data['magn_mass'], mode='markers+lines', 
                                 name='RTSIRM cooling', marker=dict(color=rtsirm_cool_color, symbol=get_plotly_marker(rtsirm_cool_marker))),
                      row=1, col=2)
        fig.add_trace(go.Scatter(x=rtsirm_warm_data['meas_temp'], y=rtsirm_warm_data['magn_mass'], mode='markers+lines', 
                                 name='RTSIRM warming', marker=dict(color=rtsirm_warm_color, symbol=get_plotly_marker(rtsirm_warm_marker))),
                      row=1, col=2)

        # Update xaxis and yaxis labels
        fig.update_xaxes(title_text="T (K)", row=1, col=1)
        fig.update_xaxes(title_text="T (K)", row=1, col=2)
        fig.update_yaxes(title_text="M (Am2/kg)", row=1, col=1)
        fig.update_yaxes(title_text="M (Am2/kg)", row=1, col=2)

        fig.update_layout(title="MPMS Data")
        fig.show()

    else:
        # Matplotlib plotting
        fig = plt.figure(figsize=(9,4.5))

        # Plotting FC and ZFC data or derivatives
        ax0 = fig.add_subplot(1,2,1)
        if plot_derivative:
            ax0.plot(fc_derivative['T'], fc_derivative['dM_dT'], marker=fc_marker, linestyle='-', color=fc_color, markersize=symbol_size, label='FC Derivative')
            ax0.plot(zfc_derivative['T'], zfc_derivative['dM_dT'], marker=zfc_marker, linestyle='-', color=zfc_color, markersize=symbol_size, label='ZFC Derivative')
        else:
            ax0.plot(fc_data['meas_temp'], fc_data['magn_mass'], marker=fc_marker, linestyle='-', color=fc_color, markersize=symbol_size, label='FC')
            ax0.plot(zfc_data['meas_temp'], zfc_data['magn_mass'], marker=zfc_marker, linestyle='-', color=zfc_color, markersize=symbol_size, label='ZFC')
        ax0.set_xlim(0, 300)
        ax0.set_ylabel('M (Am2/kg)' if not plot_derivative else 'dM/dT')
        ax0.legend()
        ax0.grid(True)

        # Plotting RTSIRM Cooling and Warming data or derivatives
        ax1 = fig.add_subplot(1,2,2)
        if plot_derivative:
            ax1.plot(rtsirm_cool_derivative['T'], rtsirm_cool_derivative['dM_dT'], marker=rtsirm_cool_marker, linestyle='-', color=rtsirm_cool_color, markersize=symbol_size, label='RTSIRM Cooling Derivative')
            ax1.plot(rtsirm_warm_derivative['T'], rtsirm_warm_derivative['dM_dT'], marker=rtsirm_warm_marker, linestyle='-', color=rtsirm_warm_color, markersize=symbol_size, label='RTSIRM Warming Derivative')
        else:
            ax1.plot(rtsirm_cool_data['meas_temp'], rtsirm_cool_data['magn_mass'], marker=rtsirm_cool_marker, linestyle='-', color=rtsirm_cool_color, markersize=symbol_size, label='RTSIRM cooling')
            ax1.plot(rtsirm_warm_data['meas_temp'], rtsirm_warm_data['magn_mass'], marker=rtsirm_warm_marker, linestyle='-', color=rtsirm_warm_color, markersize=symbol_size, label='RTSIRM warming')
        ax1.set_xlim(0, 300)
        ax1.set_ylabel('M (Am2/kg)' if not plot_derivative else 'dM/dT')
        ax1.set_xlabel('T (K)')
        ax1.legend()
        ax1.grid(True)

        fig.tight_layout()
        plt.show()
        
        
def thermomag_derivative(temps, mags):
    """
    Calculates the derivative of magnetization with respect to temperature.

    This function computes the first derivative of magnetization (M) with respect to 
    temperature (T). It takes into account the changes in magnetization and temperature
    to produce a derivative curve, which is essential in thermomagnetic analysis.

    Parameters:
        temps (pd.Series): A pandas Series representing the temperatures at which 
                           magnetization measurements were taken.
        mags (pd.Series): A pandas Series representing the magnetization measurements.

    Returns:
        pd.DataFrame: A pandas DataFrame with two columns:
                      'T' - Midpoint temperatures for each temperature interval.
                      'dM_dT' - The derivative of magnetization with respect to temperature.

    Example:
        >>> temps = pd.Series([100, 150, 200, 250, 300])
        >>> mags = pd.Series([1.2, 1.5, 1.8, 2.0, 2.2])
        >>> result = thermomag_derivative(temps, mags)
    """
    temps = temps.reset_index(drop=True)
    mags = mags.reset_index(drop=True)
    
    dT = temps.diff()
    dM = mags.diff()

    dM_dT = dM/dT
    dM_dT_real = dM_dT[1:]
    dM_dT_real.reset_index(drop = True, inplace=True)

    temps_dM_dT = []
    for n in range(len(temps)-1):
        temps_dM_dT.append(temps[n] + dT[n+1]/2)
    temps_dM_dT = pd.Series(temps_dM_dT)

    dM_dT_df = pd.DataFrame({'T':temps_dM_dT,'dM_dT':dM_dT_real})
    
    return dM_dT_df


def verwey_estimate(temps, mags, t_range_background_min = 50,
           t_range_background_max = 250,
           excluded_t_min=75,
           excluded_t_max = 150,
           poly_deg = 3):
    
    dM_dT_df = thermomag_derivative(temps, mags)
        
    fig = plt.figure(figsize=(7,7))
    ax0 = fig.add_subplot(1,1,1)
    ax0.plot(dM_dT_df['T'],dM_dT_df['dM_dT'],'.-',color='red',label='total')

    temps_filtered_indices = [i for i in np.arange(len(temps_dM_dT)) if ((float(temps_dM_dT[i]) > float(t_range_background_min)) and (float(temps_dM_dT[i])  < float(excluded_t_min)) ) or ((float(temps_dM_dT[i]) > float(excluded_t_max)) and (float(temps_dM_dT[i])  < float(t_range_background_max)))]
    temps_filtered = dM_dT_df['T'][temps_filtered_indices]
    dM_dT_filtered = dM_dT_df['dM_dT'][temps_filtered_indices]

    # ax0.plot(temps_filtered,dM_dT_filtered,'.-',color='grey',label='total (filtered)')

    poly_background_fit = np.polyfit(temps_filtered, dM_dT_filtered, poly_deg)
    dM_dT_filtered_polyfit = np.poly1d(poly_background_fit)(temps_filtered)

    residuals = dM_dT_filtered - dM_dT_filtered_polyfit
    ss_tot = np.sum((dM_dT_filtered - np.mean(dM_dT_filtered)) ** 2)
    ss_res = np.sum(residuals ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    temps_background_indices = [i for i in np.arange(len(temps_dM_dT)) if ((float(temps_dM_dT[i]) > float(t_range_background_min)) and (float(temps_dM_dT[i])  < float(t_range_background_max)))]
    temps_background_indices = [i for i in np.arange(len(temps_dM_dT)) if ((float(temps_dM_dT[i]) > float(t_range_background_min)) and (float(temps_dM_dT[i])  < float(t_range_background_max)))]
    temps_background = dM_dT_df['T'][temps_background_indices]
    temps_background.reset_index(drop = True, inplace=True)
    dM_dT_background = dM_dT_df['dM_dT'][temps_background_indices]

    dM_dT_polyfit = np.poly1d(poly_background_fit)(temps_background)

    mgt_dM_dT = dM_dT_polyfit - dM_dT_background 
    mgt_dM_dT.reset_index(drop = True, inplace=True)

    max_mgt_dM_dT = max(mgt_dM_dT)
    max_mgt_dM_dT_index = mgt_dM_dT.argmax()
    verwey_estimate = temps_background[max_mgt_dM_dT_index]

    ax0.plot(temps_background,dM_dT_polyfit,'.-',color='green',label='background fit')
    ax0.plot(temps_background,mgt_dM_dT,'.-',color='blue',label='magnetic (background fit minus total)')
    ax0.scatter(verwey_estimate,max_mgt_dM_dT,color='black',label='Verwey temperature estimate')
    ax0.set_ylabel('dM/dT (Am2/kg/K)')
    ax0.set_xlabel('T (K)')
    ax0.legend(loc='lower right')
    ax0.grid(True)
    ax0.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    plt.show()

    print('The T range for background fit is: ' + str(t_range_background_min) + ' K to ' + str(t_range_background_max) + ' K')
    print('The excluded T range is: ' + str(excluded_t_min) + ' K to ' + str(excluded_t_max) + ' K')
    print('The polynomial degree for background fit is: ' + str(poly_deg))
    print('The r-squared value for the background fit is: ' + str(round(r_squared,3)))
    print('The Verwey temperature estimate is: ' + str(round(verwey_estimate,1)) + ' K')
    
    return verwey_estimate, r_squared

# Function that needs to be further developed that could be used to read in the data from the compact MagIC format
# def read_and_process_compact_format(file_path):
#     # Initialize variables to store asterisk lines and data lines
#     asterisk_lines = []
#     data_lines = []

#     # Read the file and separate lines starting with '*' and others
#     with open(file_path, 'r') as file:
#         for line in file:
#             if line.strip().startswith('*'):
#                 asterisk_lines.append(line.strip().strip('*').split('\t'))
#             else:
#                 data_lines.append(line)

#     # Create a DataFrame from the data lines
#     data_string = ''.join(data_lines)
#     data_df = pd.read_csv(StringIO(data_string), sep='\t', skiprows=1)

#     # Process asterisk lines to add them as columns in the DataFrame
#     for item in asterisk_lines:
#         if len(item) == 2:  # Assuming each asterisk line has exactly two elements
#             data_df[item[0]] = item[1]
#         else:
#             print(f"Warning: Line '{item}' does not have two elements and will be skipped")

#     return data_df