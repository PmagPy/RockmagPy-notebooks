import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import ipywidgets as widgets
from IPython.display import display
import math
import scipy.special as sc
from scipy.optimize import brent
from numpy.linalg import svd, lstsq, solve

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
            - fc_data: Data filtered for 'LP-FC' method if available, otherwise an empty DataFrame.
            - zfc_data: Data filtered for 'LP-ZFC' method if available, otherwise an empty DataFrame.
            - rtsirm_cool_data: Data filtered for 'LP-CW-SIRM:LP-MC' method if available, otherwise an empty DataFrame.
            - rtsirm_warm_data: Data filtered for 'LP-CW-SIRM:LP-MW' method if available, otherwise an empty DataFrame.

    Example:
        >>> fc, zfc, rtsirm_cool, rtsirm_warm = extract_mpms_data(measurements_df, 'Specimen_1')
    """

    specimen_df = df[df['specimen'] == specimen_name]

    fc_data = specimen_df[specimen_df['method_codes'].str.contains('LP-FC', na=False)]
    zfc_data = specimen_df[specimen_df['method_codes'].str.contains('LP-ZFC', na=False)]
    rtsirm_cool_data = specimen_df[specimen_df['method_codes'].str.contains('LP-CW-SIRM:LP-MC', na=False)]
    rtsirm_warm_data = specimen_df[specimen_df['method_codes'].str.contains('LP-CW-SIRM:LP-MW', na=False)]

    return fc_data, zfc_data, rtsirm_cool_data, rtsirm_warm_data

def extract_hysteresis_data(df, specimen_name):
    """
    Extracts and separates Hysteresis data 
    for a specific specimen from a dataframe.

    This function filters data for a given specimen and separates it based on 
    different MagIC measurement method codes. It specifically looks for data 
    corresponding to 'LP-HYS' Hysteresis Data

    Parameters:
        df (pandas.DataFrame): The dataframe containing MPMS measurement data.
        specimen_name (str): The name of the specimen to filter data for.

    Returns:
        tuple: A tuple containing four pandas.DataFrames:
            - fc_data: Data filtered for 'LP-FC' method if available, otherwise an empty DataFrame.
            - zfc_data: Data filtered for 'LP-ZFC' method if available, otherwise an empty DataFrame.
            - rtsirm_cool_data: Data filtered for 'LP-CW-SIRM:LP-MC' method if available, otherwise an empty DataFrame.
            - rtsirm_warm_data: Data filtered for 'LP-CW-SIRM:LP-MW' method if available, otherwise an empty DataFrame.

    Example:
        >>> fc, zfc, rtsirm_cool, rtsirm_warm = extract_mpms_data(measurements_df, 'Specimen_1')
    """

    specimen_df = df[df['specimen'] == specimen_name]

    hyst_data = specimen_df[specimen_df['method_codes'].str.contains('LP-HYS', na=False)]
	
    data_type = determine_magn_data_type(hyst_data)

    print('Data Extracted!')

    return hyst_data, data_type

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
                   fc_color='#1f77b4', zfc_color='#ff7f0e', rtsirm_cool_color='#17becf', rtsirm_warm_color='#d62728',
                   fc_marker='.', zfc_marker='.', rtsirm_cool_marker='.', rtsirm_warm_marker='.',
                   symbol_size=10, use_plotly=False, plot_derivative=False, return_figure=False,
                   drop_first=False, drop_last=False):
    """
    Plots MPMS data and optionally its derivatives for Field Cooled, Zero Field Cooled, RTSIRM Cooling, and RTSIRM Warming using either Matplotlib or Plotly.

    Parameters:
        fc_data, zfc_data, rtsirm_cool_data, rtsirm_warm_data (DataFrame): DataFrames containing the MPMS data.
        fc_color, zfc_color, rtsirm_cool_color, rtsirm_warm_color (str): HEX color codes for each plot.
        fc_marker, zfc_marker, rtsirm_cool_marker, rtsirm_warm_marker (str): Marker symbols for each plot.
        symbol_size (int): Size of the markers in matplotlib, symbol size in plotly is fixed.
        use_plotly (bool): If True, uses Plotly for plotting. Otherwise, uses Matplotlib.
        plot_derivative (bool): If True, plots the derivative of the magnetization data.
        
    Returns:
        fig: The matplotlib.figure.Figure object containing the plot (only when using Matplotlib).
    """
    if drop_first:
        fc_data = fc_data[1:]
        zfc_data = zfc_data[1:]
        rtsirm_cool_data = rtsirm_cool_data[1:]
        rtsirm_warm_data = rtsirm_warm_data[1:]
        
    if drop_last:
        fc_data = fc_data[:-1]
        zfc_data = zfc_data[:-1]
        rtsirm_cool_data = rtsirm_cool_data[:-1]
        rtsirm_warm_data = rtsirm_warm_data[:-1]
        
    if plot_derivative:
        fc_derivative = thermomag_derivative(fc_data['meas_temp'], fc_data['magn_mass'])
        zfc_derivative = thermomag_derivative(zfc_data['meas_temp'], zfc_data['magn_mass'])
        rtsirm_cool_derivative = thermomag_derivative(rtsirm_cool_data['meas_temp'], rtsirm_cool_data['magn_mass'])
        rtsirm_warm_derivative = thermomag_derivative(rtsirm_warm_data['meas_temp'], rtsirm_warm_data['magn_mass'])

    if use_plotly:
        rows, cols = (2, 2) if plot_derivative else (1, 2)
        fig = make_subplots(rows=rows, cols=cols)
        
        # Add original data traces
        fig.add_trace(go.Scatter(x=fc_data['meas_temp'], y=fc_data['magn_mass'], mode='markers+lines', name='FC', marker=dict(color=fc_color)), row=1, col=1)
        fig.add_trace(go.Scatter(x=zfc_data['meas_temp'], y=zfc_data['magn_mass'], mode='markers+lines', name='ZFC', marker=dict(color=zfc_color)), row=1, col=1)
        fig.add_trace(go.Scatter(x=rtsirm_cool_data['meas_temp'], y=rtsirm_cool_data['magn_mass'], mode='markers+lines', name='RTSIRM Cooling', marker=dict(color=rtsirm_cool_color)), row=1, col=2)
        fig.add_trace(go.Scatter(x=rtsirm_warm_data['meas_temp'], y=rtsirm_warm_data['magn_mass'], mode='markers+lines', name='RTSIRM Warming', marker=dict(color=rtsirm_warm_color)), row=1, col=2)

        # Add derivative data traces if required
        if plot_derivative:
            fig.add_trace(go.Scatter(x=fc_derivative['T'], y=fc_derivative['dM_dT'], mode='markers+lines', 
                                    name='FC Derivative', marker=dict(color=fc_color, symbol=get_plotly_marker(fc_marker))),
                                    row=2, col=1)
            fig.add_trace(go.Scatter(x=zfc_derivative['T'], y=zfc_derivative['dM_dT'], mode='markers+lines', 
                                    name='ZFC Derivative', marker=dict(color=zfc_color, symbol=get_plotly_marker(zfc_marker))),
                        row=2, col=1)
            fig.add_trace(go.Scatter(x=rtsirm_cool_derivative['T'], y=rtsirm_cool_derivative['dM_dT'], mode='markers+lines', 
                                    name='RTSIRM Cooling Derivative', marker=dict(color=rtsirm_cool_color, symbol=get_plotly_marker(rtsirm_cool_marker))),
                                    row=2, col=2)
            fig.add_trace(go.Scatter(x=rtsirm_warm_derivative['T'], y=rtsirm_warm_derivative['dM_dT'], mode='markers+lines', 
                                    name='RTSIRM Warming Derivative', marker=dict(color=rtsirm_warm_color, symbol=get_plotly_marker(rtsirm_warm_marker))),
                                    row=2, col=2)
        
        # Update layout and axis titles
        # Set y-axis label for the first row to 'M (Am2/kg)'
        fig.update_yaxes(title_text="M (Am2/kg)", row=1, col=1)
        fig.update_yaxes(title_text="M (Am2/kg)", row=1, col=2)

        # Set y-axis label for the second row to 'dM/dT' if plot_derivative is True
        if plot_derivative:
            fig.update_yaxes(title_text="dM/dT", row=2, col=1)
            fig.update_yaxes(title_text="dM/dT", row=2, col=2)

        fig.update_xaxes(title_text="T (K)", row=1, col=1)
        if plot_derivative:
            # Ensure the x-axis label is applied to the bottom row, regardless of whether it's the first or second row
            fig.update_xaxes(title_text="T (K)", row=2, col=1)
            fig.update_xaxes(title_text="T (K)", row=2, col=2)
            fig.update_layout(height=900)
        else:
            fig.update_xaxes(title_text="T (K)", row=1, col=2)
            fig.update_layout(height=450)

        fig.update_layout(title="MPMS Data and Derivatives" if plot_derivative else "MPMS Data")
        
        fig.show()
            
    else:
        # Matplotlib plotting
        if plot_derivative:
            fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
        else:
            fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
            
        # Plot original data
        if not plot_derivative:
            axs[0].plot(fc_data['meas_temp'], fc_data['magn_mass'], color=fc_color, marker=fc_marker, linestyle='-', markersize=symbol_size, label='FC')
            axs[0].plot(zfc_data['meas_temp'], zfc_data['magn_mass'], color=zfc_color, marker=zfc_marker, linestyle='-', markersize=symbol_size, label='ZFC')
            axs[1].plot(rtsirm_cool_data['meas_temp'], rtsirm_cool_data['magn_mass'], color=rtsirm_cool_color, marker=rtsirm_cool_marker, linestyle='-', markersize=symbol_size, label='RTSIRM Cooling')
            axs[1].plot(rtsirm_warm_data['meas_temp'], rtsirm_warm_data['magn_mass'], color=rtsirm_warm_color, marker=rtsirm_warm_marker, linestyle='-', markersize=symbol_size, label='RTSIRM Warming')
            for ax in axs:
                ax.set_xlabel("Temperature (K)")
                ax.set_ylabel("Magnetization (Am$^2$/kg)")
                ax.legend()
                ax.grid(True)
                ax.set_xlim(0, 300)
        
        elif plot_derivative:
            axs[0, 0].plot(fc_data['meas_temp'], fc_data['magn_mass'], color=fc_color, marker=fc_marker, linestyle='-', markersize=symbol_size, label='FC')
            axs[0, 0].plot(zfc_data['meas_temp'], zfc_data['magn_mass'], color=zfc_color, marker=zfc_marker, linestyle='-', markersize=symbol_size, label='ZFC')
            axs[0, 1].plot(rtsirm_cool_data['meas_temp'], rtsirm_cool_data['magn_mass'], color=rtsirm_cool_color, marker=rtsirm_cool_marker, linestyle='-', markersize=symbol_size, label='RTSIRM Cooling')
            axs[0, 1].plot(rtsirm_warm_data['meas_temp'], rtsirm_warm_data['magn_mass'], color=rtsirm_warm_color, marker=rtsirm_warm_marker, linestyle='-', markersize=symbol_size, label='RTSIRM Warming')
            for ax in axs[0, :]:
                ax.set_xlabel("Temperature (K)")
                ax.set_ylabel("Magnetization (Am$^2$/kg)")
                ax.legend()
                ax.grid(True)
                ax.set_xlim(0, 300)
                
        if plot_derivative:
            axs[1, 0].plot(fc_derivative['T'], fc_derivative['dM_dT'], marker=fc_marker, linestyle='-', color=fc_color, markersize=symbol_size, label='FC Derivative')
            axs[1, 0].plot(zfc_derivative['T'], zfc_derivative['dM_dT'], marker=zfc_marker, linestyle='-', color=zfc_color, markersize=symbol_size, label='ZFC Derivative')
            axs[1, 1].plot(rtsirm_cool_derivative['T'], rtsirm_cool_derivative['dM_dT'], marker=rtsirm_cool_marker, linestyle='-', color=rtsirm_cool_color, markersize=symbol_size, label='RTSIRM Cooling Derivative')
            axs[1, 1].plot(rtsirm_warm_derivative['T'], rtsirm_warm_derivative['dM_dT'], marker=rtsirm_warm_marker, linestyle='-', color=rtsirm_warm_color, markersize=symbol_size, label='RTSIRM Warming Derivative')
            
            for ax in axs[1, :]:
                ax.set_xlabel("Temperature (K)")
                ax.set_ylabel("dM/dT")
                ax.legend()
                ax.grid(True)
                ax.set_xlim(0, 300)

        fig.tight_layout()
        plt.show()
        
        if return_figure:
            return fig


def plot_hyst_data(hyst_data=None,
                   fields=None,
                   moments=None,
                   data_type=None,
                   hyst_color='#1f77b4',
                   hyst_marker='.',
                   symbol_size=5, use_plotly=False, return_figure=False, double_plot=False, **kwargs):
    """
    Plots hysteresis data using either Matplotlib or Plotly.

    Parameters:
        hyst_data (DataFrame): DataFrames containing the hysteresis data.
        hyst_color: HEX color codes for each plot.
        hyst_marker (str): Marker symbols for each plot.
        symbol_size (int): Size of the markers in matplotlib, symbol size in plotly is fixed.
        use_plotly (bool): If True, uses Plotly for plotting. Otherwise, uses Matplotlib.
        
    Returns:
        fig: The matplotlib.figure.Figure object containing the plot (only when using Matplotlib).
    """

    data_units_dict_mplt = {'magn_mass': 'Am$^2$/kg', 'magn_moment': 'Am$^2$', 'magn_volume': 'A/m'}
    data_units_dict_plotly = {'magn_mass': 'Am2/kg', 'magn_moment': 'Am2', 'magn_volume': 'A/m'}
    
    # for multiple plots
    n_plots = double_plot + 1
    fields1 = kwargs.get('fields1', None)
    moments1 = kwargs.get('moments1', None)

    if use_plotly:
        rows, cols = (1, n_plots)
        fig = make_subplots(rows=rows, cols=cols)
        
        # Add original data traces
        if hyst_data is not None:
            fig.add_trace(go.Scatter(x=hyst_data['meas_field_dc'], y=hyst_data[data_type], mode='markers+lines', name='Hyst', marker=dict(color=hyst_color)), row=1,    col=1)
        elif fields is not None:
            fig.add_trace(go.Scatter(x=fields, y=moments, mode='markers+lines', name='Hyst', marker=dict(color=hyst_color)), row=1,    col=1)

        if double_plot == True:
            fig.add_trace(go.Scatter(x=fields1, y=moments1, mode='markers+lines', name='Hyst', marker=dict(color=hyst_color)), row=1,    col=2)

        
        # Update layout and axis titles
        # Set y-axis label for the first row to 'M (Am2/kg)'
        fig.update_yaxes(title_text='M '+ data_units_dict_plotly[data_type], row=1, col=1)
		
        if double_plot == True:
            fig.update_yaxes(title_text='M '+ data_units_dict_plotly[data_type], row=1, col=2)

        fig.update_xaxes(title_text="Field (Tesla)", row=1, col=1)


        fig.update_layout(title="Hysteresis Data")
        
        fig.show()
            
    else:
        # Matplotlib plotting
        fig, axs = plt.subplots(nrows=1, ncols=n_plots, figsize=(6 * n_plots, 4))
            
        # Plot original data
        if n_plots == 2:
            ax0 = axs[0]
        else:
            ax0 = axs


        if hyst_data is not None:
            ax0.plot(hyst_data['meas_field_dc'], hyst_data[data_type], color=hyst_color, marker=hyst_marker, linestyle='-', markersize=symbol_size, label='Loop')
        elif fields is not None:
            ax0.plot(fields, moments, color=hyst_color, marker=hyst_marker, linestyle='-', markersize=symbol_size, label='Loop')

		
        #for ax in axs:
        ax0.set_xlabel("Field (Tesla)")
        ax0.set_ylabel("Magnetization " + data_units_dict_mplt[data_type])
        ax0.legend()
        ax0.grid(True)
            #ax.set_xlim(0, 300)
    
        if double_plot == True:			
            axs[1].plot(fields1, moments1, color=hyst_color, marker=hyst_marker, linestyle='-', markersize=symbol_size, label='Loop')
            axs[1].set_xlabel("Field (Tesla)")
            axs[1].set_ylabel("Magnetization " + data_units_dict_mplt[data_type])
            axs[1].legend()
            axs[1].grid(True)


        fig.tight_layout()
        
        if return_figure:
            return fig, axs
        else:
            plt.show()


def make_mpms_plots(measurements):
    """
    Create a UI for specimen selection and dynamically update MPMS plots based on the selected
    specimen and plot library choice. This version adds event handlers to ensure updates occur
    upon initial selection.

    Parameters:
    experiments : pandas.DataFrame
        The dataframe containing experiment data with columns including 'specimen' and 'method_codes'.
    measurements : pandas.DataFrame
        The dataframe containing measurement data used for plotting MPMS data.
    """
    # Filter to get specimens with desired method codes
    experiments = measurements.groupby(['specimen', 'method_codes']).size().reset_index().iloc[:, :2]
    filtered_experiments = experiments[experiments['method_codes'].isin(['LP-FC', 'LP-ZFC'])]
    specimen_options = filtered_experiments['specimen'].unique().tolist()

    # Dropdown for specimen selection
    specimen_dropdown = widgets.Dropdown(
        options=specimen_options,
        description='Specimen:',
        value=specimen_options[0]
    )

    # Radio buttons for plot library choice
    plot_choice = widgets.RadioButtons(
        options=[('matplotlib', False), ('plotly', True)],
        description='Plot with:',
        disabled=False
    )

    # Interactive output container
    out = widgets.Output()

    def update_mpms_plots(specimen_name, use_plotly):
        """
        Update MPMS plots based on the selected specimen and plotting library choice.
        """
        with out:
            out.clear_output(wait=True)
            fc_data, zfc_data, rtsirm_cool_data, rtsirm_warm_data = extract_mpms_data(measurements, specimen_name)
            plot_mpms_data(fc_data, zfc_data, rtsirm_cool_data, rtsirm_warm_data, use_plotly=use_plotly, plot_derivative=True)

    def on_specimen_change(change):
        update_mpms_plots(change['new'], plot_choice.value)

    def on_plot_choice_change(change):
        update_mpms_plots(specimen_dropdown.value, change['new'])

    specimen_dropdown.observe(on_specimen_change, names='value')
    plot_choice.observe(on_plot_choice_change, names='value')

    # Initial plot to ensure something is displayed right away
    update_mpms_plots(specimen_dropdown.value, plot_choice.value)

    # Display UI components
    display(specimen_dropdown, plot_choice, out)
           

def make_hyst_plots(measurements):
    """
    Create a UI for specimen selection and dynamically update Hysteresis loop plots based on the selected
    specimen and plot library choice. This version adds event handlers to ensure updates occur
    upon initial selection.

    Parameters:
    experiments : pandas.DataFrame
        The dataframe containing experiment data with columns including 'specimen' and 'method_codes'.
    measurements : pandas.DataFrame
        The dataframe containing measurement data used for plotting Hysteresis loop data.
    """
    # Filter to get specimens with desired method codes
    experiments = measurements.groupby(['specimen', 'method_codes']).size().reset_index().iloc[:, :2]
    filtered_experiments = experiments[experiments['method_codes'].isin(['LP-HYS',])]
    specimen_options = filtered_experiments['specimen'].unique().tolist()

    # Dropdown for specimen selection
    specimen_dropdown = widgets.Dropdown(
        options=specimen_options,
        description='Specimen:',
        value=specimen_options[0]
    )

    # Radio buttons for plot library choice
    plot_choice = widgets.RadioButtons(
        options=[('matplotlib', False), ('plotly', True)],
        description='Plot with:',
        disabled=False
    )

	
    # Interactive output container
    out = widgets.Output()

    def det_data_type(data):
        if data['magn_mass'].isnull().all() == False:
            data_type = 'magn_mass'
        elif data['magn_moment'].isnull().all() == False:
            data_type = 'mag_moment'
        else:
            data_type = 'magn_volume'
		
        return data_type

    def update_hyst_plots(specimen_name, use_plotly):#, data_type):
        """
        Update hysteresis loop based on the selected specimen and plotting library choice.
        """
        with out:
            out.clear_output(wait=True)
            hyst_data = extract_hysteresis_data(measurements, specimen_name)
            data_type = det_data_type(hyst_data)
			
            plot_hyst_data(hyst_data, use_plotly=use_plotly, data_type=data_type)

        return data_type

    def on_specimen_change(change):
        data_type = update_hyst_plots(change['new'], plot_choice.value)#, data_type_choice.value)
        return data_type

    def on_plot_choice_change(change):
        data_type = update_hyst_plots(specimen_dropdown.value, change['new'])#, data_type_choice.value)
        return data_type


    specimen_dropdown.observe(on_specimen_change, names='value')
    plot_choice.observe(on_plot_choice_change, names='value')

    # Initial plot to ensure something is displayed right away
    update_hyst_plots(specimen_dropdown.value, plot_choice.value)

    # Display UI components
    display(specimen_dropdown, plot_choice, out)


def determine_magn_data_type(data):
    '''
    Search through data and determine which magnetization column is populatat. This function will choose the first column from a tiered list that is not null. The options in preferential order are: mass normalized, moment, volume ormalized.
	
    Inputs: 
    data : pandas.DataFrame
	    dataframe with MagIC database columns names
    Outputs: 
        data_type : str
            column name of non-null 
	'''

    if data['magn_mass'].isnull().all() == False:
        data_type = 'magn_mass'
    elif data['magn_moment'].isnull().all() == False:
        data_type = 'mag_moment'
    else:
        data_type = 'magn_volume'

    return data_type

def make_hyst_plots2(measurements):
    """
    Create a UI for specimen selection and dynamically update Hysteresis loop plots based on the selected
    specimen and plot library choice. This version adds event handlers to ensure updates occur
    upon initial selection.

    Parameters:
    experiments : pandas.DataFrame
        The dataframe containing experiment data with columns including 'specimen' and 'method_codes'.
    measurements : pandas.DataFrame
        The dataframe containing measurement data used for plotting Hysteresis loop data.
		
	Output:
	specimen : ipywidgets.widgets.widget_selection.Dropdown
	    dropdown specimen widget variable
    """

    # choose your specimen of choice
    method, specimen, experiement = interactive_method_specimen_selection_hyst(measurements)
	
	# Radio buttons for plot library choice
    plot_choice = widgets.RadioButtons(
    options=[('matplotlib', False), ('plotly', True)],
    description='Plot with:',
    disabled=False)
	
    # Interactive output container
    out = widgets.Output()
	
	#plot the chosen function
    def update_hyst_plots(specimen_name, use_plotly):#, data_type):
        """
        Update hysteresis loop based on the selected specimen and plotting library choice.
        """
        with out:
            out.clear_output(wait=True)
            hyst_data, data_type = extract_hysteresis_data(measurements, specimen_name)
			
            plot_hyst_data(hyst_data, use_plotly=use_plotly, data_type=data_type)
	
	# functions to update plot
    def on_specimen_change(change):
        update_hyst_plots(change['new'], plot_choice.value)

    def on_plot_choice_change(change):
        update_hyst_plots(specimen.value, change['new'])
	
    # Monitor for changes in selection	
    specimen.observe(on_specimen_change, names='value')
    plot_choice.observe(on_plot_choice_change, names='value')
	
	# initial plotting of th efirst specimen
    update_hyst_plots(specimen.value, plot_choice.value)
	
	# elements to display
    display(plot_choice, out)
	
    return specimen
	

    
def thermomag_derivative(temps, mags, drop_first=False, drop_last=False):
    """
    Calculates the derivative of magnetization with respect to temperature and optionally
    drops the data corresponding to the highest and/or lowest temperature.

    Parameters:
        temps (pd.Series): A pandas Series representing the temperatures at which
                           magnetization measurements were taken.
        mags (pd.Series): A pandas Series representing the magnetization measurements.
        drop_last (bool): Optional; whether to drop the last row from the resulting
                          DataFrame. Defaults to False. Useful when there is an
                          artifact associated with the end of the experiment.
        drop_first (bool): Optional; whether to drop the first row from the resulting
                           DataFrame. Defaults to False. Useful when there is an
                          artifact associated with the start of the experiment.

    Returns:
        pd.DataFrame: A pandas DataFrame with two columns:
                      'T' - Midpoint temperatures for each temperature interval.
                      'dM_dT' - The derivative of magnetization with respect to temperature.
                      If drop_last is True, the last temperature point is excluded.
                      If drop_first is True, the first temperature point is excluded.
    """
    temps = temps.reset_index(drop=True)
    mags = mags.reset_index(drop=True)
    
    dT = temps.diff()
    dM = mags.diff()
    dM_dT = dM / dT
    dM_dT_real = dM_dT[1:]
    dM_dT_real.reset_index(drop=True, inplace=True)

    temps_dM_dT = [temps[n] + dT[n + 1] / 2 for n in range(len(temps) - 1)]
    temps_dM_dT = pd.Series(temps_dM_dT)

    dM_dT_df = pd.DataFrame({'T': temps_dM_dT, 'dM_dT': dM_dT_real})

    # Drop the last row if specified
    if drop_last:
        dM_dT_df = dM_dT_df[:-1].reset_index(drop=True)
    
    # Drop the first row if specified
    if drop_first:
        dM_dT_df = dM_dT_df[1:].reset_index(drop=True)
    
    return dM_dT_df


def zero_crossing(dM_dT_temps, dM_dT, make_plot=False, xlim=None):
    """
    Calculate the temperature at which the second derivative of magnetization with respect to 
    temperature crosses zero. This value provides an estimate of the peak of the derivative 
    curve that is more precise than the maximum value.

    The function computes the second derivative of magnetization (dM/dT) with respect to 
    temperature, identifies the nearest points around the maximum value of the derivative, 
    and then calculates the temperature at which this second derivative crosses zero using 
    linear interpolation.

    Parameters:
        dM_dT_temps (pd.Series): A pandas Series representing temperatures corresponding to
                                 the first derivation of magnetization with respect to temperature.
        dM_dT (pd.Series): A pandas Series representing the first derivative of 
                           magnetization with respect to temperature.
        make_plot (bool, optional): If True, a plot will be generated. Defaults to False.
        xlim (tuple, optional): A tuple specifying the x-axis limits for the plot. Defaults to None.

    Returns:
        float: The estimated temperature at which the second derivative of magnetization 
               with respect to temperature crosses zero.

    Note:
        The function assumes that the input series `dM_dT_temps` and `dM_dT` are related to 
        each other and are of equal length.
    """    
    
    max_dM_dT_temp = dM_dT_temps[dM_dT.idxmax()]
    
    d2M_dT2 = thermomag_derivative(dM_dT_temps, dM_dT)
    d2M_dT2_T_array = d2M_dT2['T'].to_numpy()
    max_index = np.searchsorted(d2M_dT2_T_array, max_dM_dT_temp)

    d2M_dT2_T_before = d2M_dT2['T'][max_index-1]
    d2M_dT2_before = d2M_dT2['dM_dT'][max_index-1]
    d2M_dT2_T_after = d2M_dT2['T'][max_index]
    d2M_dT2_after = d2M_dT2['dM_dT'][max_index]

    zero_cross_temp = d2M_dT2_T_before + ((d2M_dT2_T_after - d2M_dT2_T_before) / (d2M_dT2_after - d2M_dT2_before)) * (0 - d2M_dT2_before)

    if make_plot:
        fig = plt.figure(figsize=(12,4))
        ax0 = fig.add_subplot(1,1,1)
        ax0.plot(d2M_dT2['T'], d2M_dT2['dM_dT'], '.-', color='purple', label='magnetite (background fit minus measurement)')
        ax0.plot(d2M_dT2_T_before, d2M_dT2_before, '*', color='red')
        ax0.plot(d2M_dT2_T_after, d2M_dT2_after, '*', color='red')
        ax0.plot(zero_cross_temp, 0, 's', color='blue')
        label = f'{zero_cross_temp:.1f} K'
        ax0.text(zero_cross_temp+2, 0, label, color='blue', 
                verticalalignment='center', horizontalalignment='left',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        ax0.set_ylabel('d$^2$M/dT$^2$')
        ax0.set_xlabel('T (K)')
        ax0.grid(True)
        ax0.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
        if xlim is not None:
            ax0.set_xlim(xlim)
        plt.show()
        
    return zero_cross_temp


def verwey_estimate(temps, mags, 
                    t_range_background_min = 50,
                    t_range_background_max = 250,
                    excluded_t_min = 75,
                    excluded_t_max = 150,
                    poly_deg = 3,
                    plot_zero_crossing = False,
                    plot_title = None):
    """
    This function estimates the Verwey transition temperature and remanence loss of magnetite from MPMS data.
    
    Parameters:
    """
    
    temps.reset_index(drop=True, inplace=True)
    mags.reset_index(drop=True, inplace=True)

    dM_dT_df = thermomag_derivative(temps, mags)
    temps_dM_dT = dM_dT_df['T']

    temps_dM_dT_filtered_indices = [i for i in np.arange(len(temps_dM_dT)) if ((float(temps_dM_dT[i]) > float(t_range_background_min)) and (float(temps_dM_dT[i])  < float(excluded_t_min)) ) or ((float(temps_dM_dT[i]) > float(excluded_t_max)) and (float(temps_dM_dT[i])  < float(t_range_background_max)))]
    temps_dM_dT_filtered = dM_dT_df['T'][temps_dM_dT_filtered_indices]
    dM_dT_filtered = dM_dT_df['dM_dT'][temps_dM_dT_filtered_indices]

    poly_background_fit = np.polyfit(temps_dM_dT_filtered, dM_dT_filtered, poly_deg)
    dM_dT_filtered_polyfit = np.poly1d(poly_background_fit)(temps_dM_dT_filtered)

    residuals = dM_dT_filtered - dM_dT_filtered_polyfit
    ss_tot = np.sum((dM_dT_filtered - np.mean(dM_dT_filtered)) ** 2)
    ss_res = np.sum(residuals ** 2)
    r_squared = 1 - (ss_res / ss_tot)

    temps_dM_dT_background_indices = [i for i in np.arange(len(temps_dM_dT)) if ((float(temps_dM_dT[i]) > float(t_range_background_min)) and (float(temps_dM_dT[i])  < float(t_range_background_max)))]
    temps_dM_dT_background = dM_dT_df['T'][temps_dM_dT_background_indices]
    temps_dM_dT_background.reset_index(drop=True, inplace=True)
    dM_dT_background = dM_dT_df['dM_dT'][temps_dM_dT_background_indices]
    dM_dT_polyfit = np.poly1d(poly_background_fit)(temps_dM_dT_background)

    mgt_dM_dT = dM_dT_polyfit - dM_dT_background 
    mgt_dM_dT.reset_index(drop = True, inplace=True)

    temps_background_indices = [i for i in np.arange(len(temps)) if ((float(temps[i]) > float(t_range_background_min)) and (float(temps[i])  < float(t_range_background_max)))]
    temps_background = temps[temps_background_indices]

    poly_func = np.poly1d(poly_background_fit)
    background_curve = np.cumsum(poly_func(temps_background) * np.gradient(temps_background))

    last_background_temp = temps_background.iloc[-1]    
    last_background_mag = background_curve[-1]
    target_temp_index = np.argmin(np.abs(temps - last_background_temp))
    mags_value = mags[target_temp_index]
    background_curve_adjusted = background_curve + (mags_value - last_background_mag)

    mags_background = mags[temps_background_indices]
    mgt_curve = mags_background - background_curve_adjusted
    
    verwey_estimate = zero_crossing(temps_dM_dT_background, mgt_dM_dT, 
                                    make_plot=plot_zero_crossing, xlim=(excluded_t_min, excluded_t_max))
    
    remanence_loss = np.trapz(mgt_dM_dT, temps_dM_dT_background)
    
    fig = plt.figure(figsize=(12,5))
    ax0 = fig.add_subplot(1,2,1)
    ax0.plot(temps, mags, '.-', color='red', label='measurement')
    ax0.plot(temps_background, background_curve_adjusted, '.-', color='green', label='background fit')
    ax0.plot(temps_background, mgt_curve, '.-', color='blue', label='magnetite (meas. minus background)')
    verwey_y_value = np.interp(verwey_estimate, temps_background, mgt_curve)
    ax0.plot(verwey_estimate, verwey_y_value, '*', color='pink', markersize=10,
         markeredgecolor='black', markeredgewidth=1,
         label='Verwey estimate' + ' (' + str(round(verwey_estimate,1)) + ' K)')
    ax0.set_ylabel('M (Am$^2$/kg)')
    ax0.set_xlabel('T (K)')
    ax0.legend(loc='upper right')
    ax0.grid(True)
    ax0.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    if plot_title is not None:
        ax0.set_title(plot_title)

    ax1 = fig.add_subplot(1,2,2)
    ax1.plot(dM_dT_df['T'], dM_dT_df['dM_dT'], '.-', color='red', label='measurement')
    ax1.plot(temps_dM_dT_background, dM_dT_polyfit, '.-', color='green', label='background fit'+ ' (r$^2$ = ' + str(round(r_squared,3)) + ')' )
    ax1.plot(temps_dM_dT_background, mgt_dM_dT, '.-', color='blue', label='magnetite (background fit minus measurement)')
    verwey_y_value = np.interp(verwey_estimate, temps_dM_dT_background, mgt_dM_dT)
    ax1.plot(verwey_estimate, verwey_y_value, '*', color='pink', markersize=10,
         markeredgecolor='black', markeredgewidth=1,
         label='Verwey estimate' + ' (' + str(round(verwey_estimate,1)) + ' K)')
    rectangle = patches.Rectangle((excluded_t_min, ax1.get_ylim()[0]), excluded_t_max - excluded_t_min, 
                                  ax1.get_ylim()[1] - ax1.get_ylim()[0], 
                                  linewidth=0, edgecolor=None, facecolor='gray', 
                                  alpha=0.3)
    ax1.add_patch(rectangle)
    rect_legend_patch = patches.Patch(color='gray', alpha=0.3, label='excluded from background fit')
    handles, labels = ax1.get_legend_handles_labels()
    handles.append(rect_legend_patch)  # Add the rectangle legend patch
    ax1.legend(handles=handles, loc='lower right')
    ax1.set_ylabel('dM/dT (Am$^2$/kg/K)')
    ax1.set_xlabel('T (K)')
    ax1.grid(True)
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(0,0))
    if plot_title is not None:
        ax1.set_title(plot_title)
    #plt.show()

    return verwey_estimate, remanence_loss


def interactive_verwey_estimate(measurements, specimen_dropdown, method_dropdown):
    
    selected_specimen_name = specimen_dropdown.value
    selected_method = method_dropdown.value

    fc_data, zfc_data, rtsirm_cool_data, rtsirm_warm_data = extract_mpms_data(measurements, selected_specimen_name)
    if selected_method == 'LP-FC':
        temps = fc_data['meas_temp']
        mags = fc_data['magn_mass']
    elif selected_method == 'LP-ZFC':
        temps = zfc_data['meas_temp']
        mags = zfc_data['magn_mass']

    # Determine a fixed width for the descriptions to align the sliders
    description_width = '250px'  # Adjust this based on the longest description
    slider_total_width = '600px'  # Total width of the slider widget including the description

    description_style = {'description_width': description_width}
    slider_layout = widgets.Layout(width=slider_total_width)  # Set the total width of the slider widget

    # Update sliders to use IntRangeSlider
    background_temp_range_slider = widgets.IntRangeSlider(
        value=[60, 250], min=0, max=300, step=1,
        description='Background Temperature Range (K):',
        layout=slider_layout, style=description_style
    )

    excluded_temp_range_slider = widgets.IntRangeSlider(
        value=[75, 150], min=0, max=300, step=1,
        description='Excluded Temperature Range (K):',
        layout=slider_layout, style=description_style
    )

    poly_deg_slider = widgets.IntSlider(
        value=3, min=1, max=5, step=1,
        description='Background Fit Polynomial Degree:',
        layout=slider_layout, style=description_style
    )

    # Function to reset sliders to initial values
    def reset_sliders(b):
        background_temp_range_slider.value = (60, 250)
        excluded_temp_range_slider.value = (75, 150)
        poly_deg_slider.value = 3

    # Create reset button
    reset_button = widgets.Button(description="Reset to Default Values", layout=widgets.Layout(width='200px'))
    reset_button.on_click(reset_sliders)
    
    title_label = widgets.Label(value='Adjust Parameters for ' + selected_specimen_name + ' ' + selected_method + ' fit')

    # Add the reset button to the UI
    ui = widgets.VBox([ 
        title_label,
        background_temp_range_slider, 
        excluded_temp_range_slider, 
        poly_deg_slider,
        reset_button
    ])

    out = widgets.interactive_output(
        lambda background_temp_range, excluded_temp_range, poly_deg: verwey_estimate(
            temps, mags, 
            background_temp_range[0], background_temp_range[1], 
            excluded_temp_range[0], excluded_temp_range[1], 
            poly_deg
        ), {
            'background_temp_range': background_temp_range_slider,
            'excluded_temp_range': excluded_temp_range_slider,
            'poly_deg': poly_deg_slider,
        }
    )

    out.layout.height = '500px'

    display(ui, out)


def interactive_verwey_specimen_method_selection(measurements):
    """
    Creates and displays dropdown widgets for selecting a specimen and the corresponding
    available method codes (specifically 'LP-FC' and 'LP-ZFC') from a given DataFrame of measurements.
    This function filters the measurements to include only those with desired method codes,
    dynamically updates the method dropdown based on the selected specimen, and organizes
    the dropdowns vertically in the UI.

    Parameters:
        measurements (pd.DataFrame): The DataFrame containing measurement data with columns
                                     'specimen' and 'method_codes'. It is expected to have
                                     at least these two columns where 'specimen' identifies
                                     the specimen name and 'method_codes' contains the method
                                     codes associated with each measurement.

    Returns:
        tuple: A tuple containing the specimen dropdown widget (`ipywidgets.Dropdown`)
               and the method dropdown widget (`ipywidgets.Dropdown`). The specimen dropdown
               allows for the selection of a specimen, and the method dropdown updates to
               display only the methods available for the selected specimen. The initial
               selection in the specimen dropdown is set to the first specimen option.

    Note:
        The method dropdown is initially populated based on the methods available for the
        first selected specimen. The available methods are specifically filtered for 'LP-FC'
        and 'LP-ZFC' codes.
    """
    # Filter to get specimens with desired method codes
    experiments = measurements.groupby(['specimen', 'method_codes']).size().reset_index().iloc[:, :2]
    filtered_experiments = experiments[experiments['method_codes'].isin(['LP-FC', 'LP-ZFC'])]
    specimen_options = filtered_experiments['specimen'].unique().tolist()

    selected_specimen_name = specimen_options[0]  # Example initial selection

    # Dropdown for specimen selection
    specimen_dropdown = widgets.Dropdown(
        options=specimen_options,
        description='Specimen:',
        value=selected_specimen_name
    )

    # Method dropdown initialized with placeholder options
    method_dropdown = widgets.Dropdown(
        description='Method:',
    )

    # Function to update method options based on selected specimen
    def update_method_options(change):
        selected_specimen = change['new']
        # Filter experiments to get methods available for the selected specimen
        available_methods = filtered_experiments[filtered_experiments['specimen'] == selected_specimen]['method_codes'].tolist()
        # Update method dropdown options and reset its value
        method_dropdown.options = available_methods
        if available_methods:
            method_dropdown.value = available_methods[0]
        else:
            method_dropdown.value = None

    # Register the update function with specimen dropdown
    specimen_dropdown.observe(update_method_options, names='value')

    # Initially populate method dropdown based on the first selected specimen
    update_method_options({'new': selected_specimen_name})

    # Creating a UI layout using VBox to organize the dropdowns vertically
    ui_layout = widgets.VBox([specimen_dropdown, method_dropdown])

    # Display the UI layout
    display(ui_layout)
    
    return specimen_dropdown, method_dropdown

def y_polate(n, xarr, yarr, xval):
    # linear interpolation/extrapolation, returns best-fit y for a specified x (e.g. for Mr)
    intercept, slope = linefit1(n, xarr, yarr)
    ycal = intercept + slope * xval

    return ycal
def x_polate(n, xarr, yarr, yval):
    intercept, slope = linefit1(n, xarr, yarr)
    xcal = (yval-intercept) / slope

    return xcal

def linefit1(n, xarr, yarr):
    sum_x = sum(xarr)
    sum_y = sum(yarr)
    sum_xy = sum(x * y for x, y in zip(xarr, yarr))
    sum_x_squared = sum(x ** 2 for x in xarr)

    slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x ** 2)
    intercept = (sum_y - slope * sum_x) / n

    return intercept, slope

def linefit2(n, xarr, yarr):
    sumx = sumy = st2 = b = ssr = sst = 0
    for i in range(n):
        sumx += xarr[i]
        sumy += yarr[i]

    xavg = sumx / n
    yavg = sumy / n

    for i in range(n):
        xdelta = xarr[i] - xavg
        st2 += xdelta ** 2
        ydelta = yarr[i] - yavg
        sst += ydelta ** 2
        b += xdelta * yarr[i]

    if st2 > 0:
        b /= st2
    else:
        b = 9E9  # slope

    a = (sumy - sumx * b) / n  # intercept

    for i in range(n):
        # ssr += (yarr[i] - (a + b * xarr[i])) ** 2  # sum squared residuals
        ssr += ((a + b * xarr[i]) - yavg) ** 2  # sum squared due to regression

    if sst > 0:
        rsqr = ssr / sst
    else:
        rsqr = 1

    slope = b
    intercept = a

    return intercept, slope, rsqr

def branch_sub(n_loop, grid_moments):

    j = n_loop // 2
    moment_sub = np.ndarray([j])
    for i in range (j):
        moment_sub[i] = grid_moments[i]-grid_moments[n_loop-1-i]
    return moment_sub

def loop_grid(n_loop, polydegree, nsmooth, loop_fields, loop_moments, B_offset, M_offset):

    n2 = nsmooth // 2  # make nsmooth odd; n2=# on each side, + 1 in center
    xvals = np.ndarray([nsmooth])
    yvals = np.ndarray([nsmooth])
    xvals[:]=0
    yvals[:]=0

    grid_fields = np.ndarray([n_loop])
    grid_moments = np.ndarray([n_loop])
    grid_fields[:] = 0
    grid_moments[:] = 0


    field_step = loop_fields[0] / n_loop * 4
    grid_fields[0] = round(1000 * (loop_fields[0] - abs(B_offset))) / 1000  # nearest mT
    grid_fields[n_loop // 2] = -round(1000 * (loop_fields[0] - abs(B_offset))) / 1000  # nearest mT
    field_step = grid_fields[0] / (n_loop - 2) * 4
    for i in range(1, n_loop // 2):  #was +1, but this is 0 bassed array
        grid_fields[i] = grid_fields[i - 1] - field_step
        grid_fields[i + n_loop // 2] = grid_fields[i - 1 + n_loop // 2] + field_step


    j = 1
    grid_moments[0] = loop_moments[0] - M_offset

    for i in range(1, n_loop // 2):
        while not (loop_fields[j] - B_offset >= grid_fields[i]): # or (j <= 2):
            if j<=1:
                break
            if j > 1:
                j -= 1
        while not (loop_fields[j] - B_offset <= grid_fields[i]): #or (j >= n_loop // 2):
            if j >= n_loop // 2:
                break
            if j < n_loop // 2:
                j += 1

        k = 1
        xvals[1] = loop_fields[j] - B_offset
        xvals[2] = loop_fields[j - k] - B_offset
        yvals[1] = loop_moments[j] - M_offset
        yvals[2] = loop_moments[j - k] - M_offset

        while abs(xvals[1] - xvals[2]) <= field_step / 2:
            k += 1
            xvals[2] = loop_fields[j - k] - B_offset
            yvals[2] = loop_moments[j - k] - M_offset

        ycal = y_polate(2, xvals, yvals, grid_fields[i])  #, grid_moments[i])
        grid_moments[i]=ycal
        #print(grid_moments[i])
    j -= 1
    grid_moments[n_loop // 2] = loop_moments[n_loop // 2] - M_offset

    for i in range(n_loop // 2, n_loop):
        while not (loop_fields[j] - B_offset <= grid_fields[i]): # or (j <= n_loop // 2 + 1):
            if j <= n_loop // 2: #+ 1:
                break
            if j > n_loop // 2: #+ 1:
                j -= 1
        while not (loop_fields[j] - B_offset >= grid_fields[i]): # or (j >= n_loop):
            if j >= n_loop -1:
                break
            if j < n_loop -1:
                j += 1

        k = 1
        xvals[1] = loop_fields[j] - B_offset
        xvals[2] = loop_fields[j - k] - B_offset
        yvals[1] = loop_moments[j] - M_offset
        yvals[2] = loop_moments[j - k] - M_offset

        while abs(xvals[1] - xvals[2]) <= field_step / 2:
            k += 1
            xvals[2] = loop_fields[j - k] - B_offset
            yvals[2] = loop_moments[j - k] - M_offset

        ycal = y_polate(2, xvals, yvals, grid_fields[i])   #, grid_moments[i])
        grid_moments[i] = ycal
        #print(grid_moments[i])
    j = 2
    """IRMDB software uses only linear fitting for gridding at this point
    if polydegree > 1:
        for i in range(n2 + 1, n_loop // 2 - (n2 + 1) + 1):
            while not (loop_fields[j] - B_offset >= grid_fields[i]): # or (j <= 2):
                if j<=2:
                    break
                if j > 2:
                    j -= 1
            while not (loop_fields[j] - B_offset <= grid_fields[i]): # or (j >= n_loop // 2):
                if j >= n_loop // 2:
                    break
                if j < n_loop // 2:
                    j += 1

            for k in range(-n2, n2 + 1):
                x[k + n2 + 1] = loop_fields[j + k] - B_offset
                y[k + n2 + 1] = loop_moments[j + k] - M_offset

            # Call to polynomialfit1 function goes here
            print(x)
            r2 = float()
            polynomialfit1(nsmooth, 0, polydegree, x, y) #polynomialfit1(nsmooth, 0, polydegree, x, y, y, r2)


            y1 = 0
            x1 = 1

            for k in range(polydegree + 1):
                if k > 0:
                    x1 *= grid_fields[i]
                y1 += y[k + 1] * x1

            grid_moments[i] = y1

        j -= 1

        for i in range(n_loop // 2 + n2 + 1, n_loop - (n2 + 1) + 1):
            while not (loop_fields[j] - B_offset <= grid_fields[i]) or (j <= n_loop // 2 + 1):
                if j > n_loop // 2 + 1:
                    j -= 1
            while not (loop_fields[j] - B_offset >= grid_fields[i]) or (j >= n_loop):
                if j < n_loop:
                    j += 1

            try:
                for k in range(-n2, n2 + 1):
                    x[k + n2 + 1] = loop_fields[j + k] - B_offset
                    y[k + n2 + 1] = loop_moments[j + k] - M_offset
            except:
                x1 = 0

            if j + k <= n_loop:
                # Call to polynomialfit1 function goes here
                # polynomialfit1(nsmooth, 0, polydegree, x, y, y, r2)
                # The implementation of polynomialfit1 function is needed

                y1 = 0
                x1 = 1

                for k in range(polydegree + 1):
                    if k > 0:
                        x1 *= grid_fields[i]
                    y1 += y[k + 1] * x1

                grid_moments[i] = y1
    """            
    return grid_fields, grid_moments

def loop_test_linearity(n_looppoints, loop_fields, loop_moments):

    intercept, slope = linefit1(n_looppoints, loop_fields, loop_moments)
    n = n_looppoints // 2
    sst = 0
    ssr = 0
    SSD = 0
    SSLF = 0
    SSPE = 0
    ymean = 0
    ymean = sum(loop_moments) / (2 * n)

    for i in range(1, n+1):   # range(1, n + 1):
        b1 = loop_fields[i - 1]
        b2 = -loop_fields[i + n - 1]
        m1 = loop_moments[i - 1]
        m2 = -loop_moments[i + n - 1]
        mfit = b1 * slope + intercept
        sst += (m1 - ymean) ** 2
        ssr += (mfit - ymean) ** 2
        SSD += (mfit - m1) ** 2
        mfit = b2 * slope + intercept
        sst += (m2 - ymean) ** 2
        ssr += (mfit - ymean) ** 2
        SSD += (mfit - m2) ** 2
        # SSPE += ((m1 - m2) / 2) ** 2  Chat GPT botched the conversion
        SSPE += ((m1 - m2) ** 2) / 2

    SSLF = SSD - SSPE
    MSR = ssr
    MSD = SSD / (2 * n - 2)
    MSPE = SSPE / n
    MSLF = SSLF / (n - 2)
    FL = MSR / MSD
    FNL = MSLF / MSPE

    return FNL, slope, intercept

def loop_errorcal(n_looppoints, loopfields, moments):

    r2, boff1, moff1 = loop_Hshift_brent(n_looppoints, loopfields, moments, -loopfields[0] / 2, 0, loopfields[0] / 2, 1E-6)
    #result = brent(objective_function, brack=(-loopfields[0] / 2, 0, loopfields[0] / 2), tol=1E-6)

    M_offset = moff1
    B_offset = boff1

    if 1 + r2 > 0:
        M_sn = math.sqrt(1 + r2)  # =sqrt(1-R^2) = sqrt(SSD/SST) ~ noise/signal
    else:
        M_sn = 0

    return M_offset, B_offset, M_sn

def loop_Hshift_brent(n_looppoints, loopfields, moments, ax, bx, cx, tol):
   # def loop_R2(n_looppoints, loop_fields, loop_moments, H_shift, V_shift):
        # Define loop_R2 function here or import it from another module
   #     pass

    itmax = 100
    cgold = 0.3819660
    zeps = 1.0E-10

    a = ax if ax < cx else cx
    b = ax if ax > cx else cx
    v = bx
    w = v
    x26 = v   #changing x to x26 because seems like python does not like me using x
    e = 0.0
    #fx = -loop_R2(n_looppoints, loop_fields, loop_moments, x)
    fx=loop_R2v2(n_looppoints, loopfields, moments, x26)[0]
    fx = fx * -1
    #print(fx)
    fv = fx
    fw = fx
    d = 0
    goto_1 = True
    goto_2 = True
    for iter in range(1, itmax + 1):
        xm = 0.5 * (a + b)
        tol1 = tol * abs(x26) + zeps
        tol2 = 2.0 * tol1

        if abs(x26 - xm) <= (tol2 - 0.5 * (b - a)):
            break

        if abs(e) > tol1:
            r = (x26 - w) * (fx - fv)
            q = (x26 - v) * (fx - fw)
            p = (x26 - v) * q - (x26 - w) * r
            q = 2.0 * (q - r)

            if q > 0.0:
                p = -p
            q = abs(q)
            etemp = e
            e = d
            goto_1 = True

            if abs(p) >= abs(0.5 * q * etemp) or p <= q * (a - x26) or p >= q * (b - x26):
                goto_1 = True
            else:
                d = p / q
                u = x26 + d

                if (u - a) < tol2 or (b - u) < tol2:
                    d = sign(tol1, xm - x26)
                goto_2 = True
                goto_1 = False

        if goto_1: #or (not goto_1 and abs(d) >= tol1):
            if x26 >= xm:
                e = a - x26
            else:
                e = b - x26
            d = cgold * e

        if goto_2: #or (not goto_2 and abs(d) >= tol1):
            if abs(d) >= tol1:
                u = x26 + d
            else:
                u = x26 + sign(tol1, d)
            fu=loop_R2v2(n_looppoints, loopfields, moments, u)[0]
            fu = fu * -1
            if fu <= fx:
                if u >= x26:
                    a = x26
                else:
                    b = x26
                v = w
                fv = fw
                w = x26
                fw = fx
                x26 = u
                fx = fu
            else:
                if u < x26:
                    a = u
                else:
                    b = u
                if fu <= fw or w == x26:
                    v = w
                    fv = fw
                    w = u
                    fw = fu
                elif fu <= fv or v == x26 or v == 2:
                    v = u
                    fv = fu

    boff1 = x26
    r2 = fx
    moff1 = loop_R2v2(n_looppoints, loopfields, moments, u)[1]   #moff1 = loop_R2v2(n_looppoints, loop_fields, loop_moments, u)[1]
    return r2, boff1, moff1

def loop_R2v2(n_looppoints, loopfields, moments, x26):
    n = n_looppoints
    min2 = 9E9
    max2 = -min2
    loop_fields_Hshift = []
    for i in range(n):
        loop_fields_Hshift.append(loopfields[i] - x26)  #loop_fields.append(float(s1))

    for i in range(n):
        if loop_fields_Hshift[i] > max2:
            max2 = loop_fields_Hshift[i]

    for i in range(n):
        if loop_fields_Hshift[i] < min2:
            min2 = loop_fields_Hshift[i]

    min1 = -max2
    max1 = -min2
    n1 = n // 2
    i2 = 0
    x1 = np.zeros(n1)
    y1 = np.zeros(n1)

    n2 = 0
    for i in range(n1, n):
        x22 = loop_fields_Hshift[i]
        if min1 < x22 < max1:
            while -loop_fields_Hshift[i2] < x22:
                i2 += 1
            if i2 > 0:
                n2 += 1
                x1[n2] = moments[i]
                dx = (-loop_fields_Hshift[i2] - x22) / (-loop_fields_Hshift[i2] + loop_fields_Hshift[i2 - 1])
                dy = dx * (-moments[i2] + moments[i2 - 1])
                y = -moments[i2] - dy
                y1[n2] = -y

    intercept, slope, rsqr = linefit2(n2, x1, y1)
    V_shift = intercept / 2
    r2 = rsqr
    return r2, V_shift

def sign(a, b):
    return abs(a) if b > 0.0 else -abs(a)

def loop_errorcal2(n_looppoints, loop_fields, loop_moments, M_offset, B_offset):
    n = n_looppoints
    min2 = 9E9
    max2 = -min2
    ErrX = np.zeros(n_looppoints)
    ErrY = np.zeros(n_looppoints)
    for i in range(n):
        loop_fields[i] -= B_offset

    for i in range(n):
        if loop_fields[i] > max2:
            max2 = loop_fields[i]

    for i in range(n):
        if loop_fields[i] < min2:
            min2 = loop_fields[i]

    min1 = -max2
    max1 = -min2
    n1 = n // 2
    i2 = 1
    n2 = 0

    for i in range(n1, n):
        x = loop_fields[i]
        y0 = loop_moments[i]
        if (x - min1 > 1E-10) and (x < max1):
            while -loop_fields[i2] < x:
                i2 += 1
            if i2 > 0:
                n2 += 1
                ErrX[n2] = -x
                dx = (-loop_fields[i2] - x) / (-loop_fields[i2] + loop_fields[i2 - 1])
                dy = dx * (-loop_moments[i2] + loop_moments[i2 - 1])
                y = -loop_moments[i2] - dy
                ErrY[n2] = y0 - y
    return ErrX, ErrY


def loop_delta_M(n_looppoints, loopfields, moments):
    # Initialize variables
    n2 = n_looppoints // 2
    xarr = np.zeros(n2 + 1)  # array needs to be filled with zeros
    yarr = np.zeros(n2 + 1)
    mrh = np.zeros(n2 + 1)

    E_hys = 0
    Mrhmax = 0
    noise1 = 0

    # Main loop
    for i in range(1, n2 + 1):
        j = n_looppoints + 51 - i
        if j > n_looppoints:
            j = n_looppoints
        while not (loopfields[j - 1] <= loopfields[i - 1] or j == n_looppoints // 2):
            j -= 1

        xarr[0] = loopfields[j - 1]
        xarr[1] = loopfields[j]
        yarr[0] = moments[j - 1]
        yarr[1] = moments[j]

        if xarr[0] == xarr[1]:
            y1 = yarr[0]
        else:
            y1 = y_polate(2, xarr, yarr, loopfields[i - 1])

        mrh[i - 1] = (moments[i - 1] - y1) / 2

        if mrh[i - 1] > Mrhmax:
            Mrhmax = mrh[i - 1]

        if i > 1:
            E_hys -= (mrh[i - 1] + mrh[i - 2]) / 2 * (loopfields[i - 1] - loopfields[i - 2])
            noise1 += (mrh[i - 1] - mrh[i - 2]) ** 2

    noise1 = np.sqrt(noise1 / n2)
    drift1 = np.sum(mrh[:10])

    # Calculate Brh, Bih
    MS = np.mean(moments[2:6])
    j = 0
    while not (moments[j] - mrh[j] <= MS / 2 or j == n_looppoints // 2):
        j += 1
    xarr[0] = loopfields[j - 1]
    xarr[1] = loopfields[j - 2]
    yarr[0] = moments[j - 1] - mrh[j - 1]
    yarr[1] = moments[j - 2] - mrh[j - 2]
    x1 = x_polate(2, xarr, yarr, MS / 2)
    Bih = x1

    j = 0
    while not (mrh[j] >= Mrhmax / 2 or j == n_looppoints // 2):
        j += 1
    xarr[0] = loopfields[j - 1]
    xarr[1] = loopfields[j - 2]
    yarr[0] = mrh[j - 1]
    yarr[1] = mrh[j - 2]
    x1 = x_polate(2, xarr, yarr, Mrhmax / 2)
    Brh = x1

    return mrh, E_hys, Brh, Bih

def Qfactor(M_sn):
    """
    A quality factor Q is calculated for both the whole curve, Q, and the ferromagnetic part, Qf.
    Q is the decimal log of the signal/noise ratio, calculated from the mean square mismatch
    between symmetrically equivalent moments, and typically ranges between 0 (s/n~1) and ~3 (s/n~1000).
    A higher value indicates better quality data.
    
    Input: (s/n)**-1/2

    """
    
    
    if M_sn > 0:
        Q = math.log(M_sn, 10) * -1 # potentially off by a factor of 2
    else:
        Q = 10
    
    return Q














def interactive_specimen_selection(measurements):
    """
    Creates and displays a dropdown widget for selecting a specimen from a given DataFrame
    of measurements. 

    Parameters:
        measurements (pd.DataFrame): The DataFrame containing measurement data with a column
                                     'specimen'. It is expected to have at least this column
                                     where 'specimen' identifies the specimen name.

    Returns:
        ipywidgets.Dropdown: A dropdown widget allowing for the selection of a specimen.
                             The initial selection in the dropdown is set to the first
                             specimen option.
    """
    # Extract unique specimen names from the measurements DataFrame
    specimen_options = measurements['specimen'].unique().tolist()

    # Set the initial selection to the first specimen option, if available
    selected_specimen_name = specimen_options[0] if specimen_options else None

    # Create a dropdown for specimen selection
    specimen_dropdown = widgets.Dropdown(
        options=specimen_options,
        description='Specimen:',
        value=selected_specimen_name
    )

    # Display the dropdown widget
    display(specimen_dropdown)

    return specimen_dropdown


def interactive_specimen_experiment_selection(measurements):
    # make two drop down ipywidgets for the user to select the sample and the associated experiment
    specimen_dropdown = widgets.Dropdown(
        options = measurements['specimen'].unique(),
        description = 'specimen:',
        disabled = False,
    )

    experiment_dropdown = widgets.Dropdown(
        options = measurements['experiment'].unique(),
        description = 'Experiment:',
        disabled = False,
    )
    # make sure to set the default value of the experiment dropdown to the first experiment in the specimen dropdown
    experiment_dropdown.options = measurements[measurements['specimen']==specimen_dropdown.value]['experiment'].unique()

    # make sure to update the experiment dropdown based on the specimen selected
    def update_experiment(*args):
        experiment_dropdown.options = measurements[measurements['specimen']==specimen_dropdown.value]['experiment'].unique()

    specimen_dropdown.observe(update_experiment, 'value')

    # display the dropdowns
    display(specimen_dropdown, experiment_dropdown)
    
    return specimen_dropdown, experiment_dropdown


def interactive_method_specimen_selection_hyst(measurements):
    # make two drop down ipywidgets for the user to select the experiment and the associated specimen

    method_dropdown = widgets.Dropdown(
        options = measurements[measurements['method_codes'].str.contains("HYS")]['method_codes'].unique(),
        description = 'Method Code:',
        disabled = False,
    )
	
    specimen_dropdown = widgets.Dropdown(
        options = measurements[measurements['method_codes'].str.contains("HYS")]['specimen'].unique(),
        description = 'Specimen:',
        disabled = False,
    )
	
    experiment_dropdown = widgets.Dropdown(
        options = measurements[measurements['method_codes'].str.contains("HYS")]['experiment'].unique(),
        description = 'Experiment:',
        disabled = False,
    )	

    # make sure to set the default value of the specimen dropdown to the first specimen in the experiment dropdown
    specimen_dropdown.options = measurements[measurements['method_codes']==method_dropdown.value]['specimen'].unique()
    experiment_dropdown.options = measurements[measurements['specimen']==specimen_dropdown.value]['experiment'].unique()

    # make sure to update the experiment dropdown based on the specimen selected
    def update_specimen(*args):
        specimen_dropdown.options = measurements[measurements['method_codes']==method_dropdown.value]['specimen'].unique()

    def update_experiment(*args):
        experiment_dropdown.options = measurements[(measurements['specimen']==specimen_dropdown.value) & (measurements['method_codes']==method_dropdown.value)]['experiment'].unique()

    method_dropdown.observe(update_specimen, 'value')
    specimen_dropdown.observe(update_experiment, 'value')

    # display the dropdowns
    display(method_dropdown, specimen_dropdown, experiment_dropdown)
    
    return method_dropdown, specimen_dropdown, experiment_dropdown



def goethite_removal(rtsirm_warm_data, 
                     rtsirm_cool_data,
                     t_min=150, t_max=290, poly_deg=2,
                     rtsirm_cool_color='#17becf', rtsirm_warm_color='#d62728',
                     symbol_size=4, return_data=False):
    """
    Analyzes and visualizes the removal of goethite signal from Room Temperature Saturation
    Isothermal Remanent Magnetization (RTSIRM) warming and cooling data. The function fits
    a polynomial to the RTSRIM warming curve between specified temperature bounds to model
    the goethite contribution, then subtracts this fit from the original data. The corrected
    and uncorrected magnetizations are plotted, along with their derivatives, to assess the
    effect of goethite removal.

    Parameters:
        rtsirm_warm_data (pd.DataFrame): DataFrame containing 'meas_temp' and 'magn_mass' columns
                                         for RTSIRM warming data.
        rtsirm_cool_data (pd.DataFrame): DataFrame containing 'meas_temp' and 'magn_mass' columns
                                         for RTSIRM cooling data.
        t_min (int, optional): Minimum temperature for polynomial fitting. Default is 150.
        t_max (int, optional): Maximum temperature for polynomial fitting. Default is 290.
        poly_deg (int, optional): Degree of the polynomial to fit. Default is 2.
        rtsirm_cool_color (str, optional): Color code for plotting cooling data. Default is '#17becf'.
        rtsirm_warm_color (str, optional): Color code for plotting warming data. Default is '#d62728'.
        symbol_size (int, optional): Size of the markers in the plots. Default is 4.
        return_data (bool, optional): If True, returns the corrected magnetization data for both
                                      warming and cooling. Default is False.

    Returns:
        Tuple[pd.Series, pd.Series]: Only if return_data is True. Returns two pandas Series
                                     containing the corrected magnetization data for the warming
                                     and cooling sequences, respectively.
    """
    
    rtsirm_warm_temps = rtsirm_warm_data['meas_temp']
    rtsirm_warm_mags = rtsirm_warm_data['magn_mass']
    rtsirm_cool_temps = rtsirm_cool_data['meas_temp']
    rtsirm_cool_mags = rtsirm_cool_data['magn_mass']
    
    rtsirm_warm_temps.reset_index(drop=True, inplace=True)
    rtsirm_warm_mags.reset_index(drop=True, inplace=True)
    rtsirm_cool_temps.reset_index(drop=True, inplace=True)
    rtsirm_cool_mags.reset_index(drop=True, inplace=True)
    
    rtsirm_warm_temps_filtered_indices = [i for i in np.arange(len(rtsirm_warm_temps)) if ((float(rtsirm_warm_temps[i]) > float(t_min)) and (float(rtsirm_warm_temps[i])  < float(t_max)) )]
    rtsirm_warm_temps_filtered = rtsirm_warm_temps[rtsirm_warm_temps_filtered_indices]
    rtsirm_warm_mags_filtered = rtsirm_warm_mags[rtsirm_warm_temps_filtered_indices]
    
    geothite_fit = np.polyfit(rtsirm_warm_temps_filtered, rtsirm_warm_mags_filtered, poly_deg)
    rtsirm_warm_mags_polyfit = np.poly1d(geothite_fit)(rtsirm_warm_temps)
    rtsirm_cool_mags_polyfit = np.poly1d(geothite_fit)(rtsirm_cool_temps)
    
    rtsirm_warm_mags_corrected = rtsirm_warm_mags - rtsirm_warm_mags_polyfit
    rtsirm_cool_mags_corrected = rtsirm_cool_mags - rtsirm_cool_mags_polyfit
    
    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(12, 8))
    
    axs[0, 0].plot(rtsirm_warm_temps, rtsirm_warm_mags, color=rtsirm_warm_color, 
                   marker='o', linestyle='-', markersize=symbol_size, label='RTSIRM Warming')
    axs[0, 0].plot(rtsirm_cool_temps, rtsirm_cool_mags, color=rtsirm_cool_color, 
                   marker='o', linestyle='-', markersize=symbol_size, label='RTSIRM Cooling')
    axs[0, 0].plot(rtsirm_warm_temps, rtsirm_warm_mags_polyfit, color=rtsirm_warm_color, 
                   linestyle='--', label='goethite fit')
    axs[0, 1].plot(rtsirm_warm_temps, rtsirm_warm_mags_corrected, color=rtsirm_warm_color, 
                   marker='s', linestyle='-', markersize=symbol_size, label='RTSIRM Warming (goethite removed)')
    axs[0, 1].plot(rtsirm_cool_temps, rtsirm_cool_mags_corrected, color=rtsirm_cool_color, 
                   marker='s', linestyle='-', markersize=symbol_size, label='RTSIRM Cooling (goethite removed)')
    
    ax0 = axs[0, 0] 
    rectangle = patches.Rectangle((t_min, ax0.get_ylim()[0]), t_max - t_min, 
                            ax0.get_ylim()[1] - ax0.get_ylim()[0], 
                            linewidth=0, edgecolor=None, facecolor='gray', 
                            alpha=0.3)
    ax0.add_patch(rectangle)
    rect_legend_patch = patches.Patch(color='gray', alpha=0.3, label='excluded from background fit')
    handles, labels = ax0.get_legend_handles_labels()
    handles.append(rect_legend_patch)  # Add the rectangle legend patch
    
    for ax in axs[0, :]:
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Magnetization (Am$^2$/kg)")
        ax.legend()
        ax.grid(True)
        ax.set_xlim(0, 300)
             
    rtsirm_cool_derivative = thermomag_derivative(rtsirm_cool_data['meas_temp'], 
                                                       rtsirm_cool_data['magn_mass'], drop_first=True)
    rtsirm_warm_derivative = thermomag_derivative(rtsirm_warm_data['meas_temp'], 
                                                       rtsirm_warm_data['magn_mass'], drop_last=True)
    
    rtsirm_cool_derivative_corrected = thermomag_derivative(rtsirm_cool_data['meas_temp'], 
                                                       rtsirm_cool_mags_corrected, drop_first=True)
    rtsirm_warm_derivative_corrected = thermomag_derivative(rtsirm_warm_data['meas_temp'], 
                                                       rtsirm_warm_mags_corrected, drop_last=True)

    axs[1, 0].plot(rtsirm_cool_derivative['T'], rtsirm_cool_derivative['dM_dT'], 
                   marker='o', linestyle='-', color=rtsirm_cool_color, markersize=symbol_size, label='RTSIRM Cooling Derivative')
    axs[1, 0].plot(rtsirm_warm_derivative['T'], rtsirm_warm_derivative['dM_dT'], 
                   marker='o', linestyle='-', color=rtsirm_warm_color, markersize=symbol_size, label='RTSIRM Warming Derivative')        
    axs[1, 1].plot(rtsirm_cool_derivative_corrected['T'], rtsirm_cool_derivative_corrected['dM_dT'], 
                   marker='s', linestyle='-', color=rtsirm_cool_color, markersize=symbol_size, label='RTSIRM Cooling Derivative\n(goethite removed)')
    axs[1, 1].plot(rtsirm_warm_derivative_corrected['T'], rtsirm_warm_derivative_corrected['dM_dT'], 
                   marker='s', linestyle='-', color=rtsirm_warm_color, markersize=symbol_size, label='RTSIRM Warming Derivative\n(goethite removed)')  
    for ax in axs[1, :]:
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("dM/dT")
        ax.legend()
        ax.grid(True)
        ax.set_xlim(0, 300)

    fig.tight_layout()
    plt.show()

    if return_data:
        rtsirm_warm_adjusted = pd.DataFrame({'meas_temp': rtsirm_warm_temps, 'corrected_magn_mass': rtsirm_warm_mags_corrected})
        rtsirm_cool_adjusted = pd.DataFrame({'meas_temp': rtsirm_cool_temps, 'corrected_magn_mass': rtsirm_cool_mags_corrected})
        return rtsirm_warm_adjusted, rtsirm_cool_adjusted
    
    
def interactive_goethite_removal(measurements, specimen_dropdown):
    
    selected_specimen_name = specimen_dropdown.value

    fc_data, zfc_data, rtsirm_cool_data, rtsirm_warm_data = extract_mpms_data(measurements, selected_specimen_name)

    # Determine a fixed width for the descriptions to align the sliders
    description_width = '250px'  # Adjust this based on the longest description
    slider_total_width = '600px'  # Total width of the slider widget including the description

    description_style = {'description_width': description_width}
    slider_layout = widgets.Layout(width=slider_total_width)  # Set the total width of the slider widget

    # Update sliders to use IntRangeSlider
    temp_range_slider = widgets.IntRangeSlider(
        value=[150, 290], min=0, max=300, step=1,
        description='Geothite Fit Temperature Range (K):',
        layout=slider_layout, style=description_style
    )

    poly_deg_slider = widgets.IntSlider(
        value=2, min=1, max=3, step=1,
        description='Goethite Fit Polynomial Degree:',
        layout=slider_layout, style=description_style
    )

    # Function to reset sliders to initial values
    def reset_sliders(b):
        temp_range_slider.value = (150, 290)
        poly_deg_slider.value = 2

    # Create reset button
    reset_button = widgets.Button(description="Reset to Default Values", layout=widgets.Layout(width='200px'))
    reset_button.on_click(reset_sliders)
    
    title_label = widgets.Label(value='Adjust Parameters for ' + selected_specimen_name + ' ' + 'goethite' + ' fit')

    # Add the reset button to the UI
    ui = widgets.VBox([ 
        title_label,
        temp_range_slider, 
        poly_deg_slider,
        reset_button
    ])

    out = widgets.interactive_output(
        lambda temp_range, poly_deg: goethite_removal(
            rtsirm_warm_data, rtsirm_cool_data, 
            temp_range[0], temp_range[1],
            poly_deg
        ), {
            'temp_range': temp_range_slider,
            'poly_deg': poly_deg_slider,
        }
    )

    out.layout.height = '500px'

    display(ui, out)

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