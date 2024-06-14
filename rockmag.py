import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.patches as patches
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import ipywidgets as widgets
from IPython.display import display
from scipy.interpolate import interp1d


try:
    from lmfit import Parameters, Model  # for fitting
    from lmfit.models import SkewedGaussianModel
except ImportError:
    Parameters = None
    Model = None
    SkewedGaussianModel = None

try:
    import statsmodels.api as sm
    lowess = sm.nonparametric.lowess
except ImportError:
    sm = None
    lowess = None


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

    return hyst_data

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


def plot_hyst_data(hyst_data,
                   hyst_color='#1f77b4',
                   hyst_marker='.',
                   symbol_size=5, use_plotly=False, return_figure=False):
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

    if use_plotly:
        rows, cols = (1, 1)
        fig = make_subplots(rows=rows, cols=cols)
        
        # Add original data traces
        fig.add_trace(go.Scatter(x=hyst_data['meas_field_dc'], y=hyst_data['magn_mass'], mode='markers+lines', name='Hyst', marker=dict(color=hyst_color)), row=1,    col=1)
        
        # Update layout and axis titles
        # Set y-axis label for the first row to 'M (Am2/kg)'
        fig.update_yaxes(title_text="M (Am2/kg)", row=1, col=1)
        fig.update_yaxes(title_text="M (Am2/kg)", row=1, col=2)

        fig.update_xaxes(title_text="Field (Tesla)", row=1, col=1)


        fig.update_layout(title="Hysteresis Data")
        
        fig.show()
            
    else:
        # Matplotlib plotting
        fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(12, 4))
            
        # Plot original data

        axs[0].plot(hyst_data['meas_field_dc'], hyst_data['magn_mass'], color=hyst_color, marker=hyst_marker, linestyle='-', markersize=symbol_size, label='Loop')
        axs[1].plot(hyst_data['meas_field_dc'], hyst_data['magn_mass'], color=hyst_color, marker=hyst_marker, linestyle='-', markersize=symbol_size, label='Loop')
        for ax in axs:
            ax.set_xlabel("Field (Tesla)")
            ax.set_ylabel("Magnetization (Am$^2$/kg)")
            ax.legend()
            ax.grid(True)
            #ax.set_xlim(0, 300)
    

        fig.tight_layout()
        plt.show()
        
        if return_figure:
            return fig


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

    def update_hyst_plots(specimen_name, use_plotly):
        """
        Update hysteresis loop based on the selected specimen and plotting library choice.
        """
        with out:
            out.clear_output(wait=True)
            hyst_data = extract_hysteresis_data(measurements, specimen_name)
            plot_hyst_data(hyst_data, use_plotly=use_plotly)

    def on_specimen_change(change):
        update_hyst_plots(change['new'], plot_choice.value)

    def on_plot_choice_change(change):
        update_hyst_plots(specimen_dropdown.value, change['new'])

    specimen_dropdown.observe(on_specimen_change, names='value')
    plot_choice.observe(on_plot_choice_change, names='value')

    # Initial plot to ensure something is displayed right away
    update_hyst_plots(specimen_dropdown.value, plot_choice.value)

    # Display UI components
    display(specimen_dropdown, plot_choice, out)

    
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

def y_polate(n, xarr, yarr, xval, yval):
    # linear interpolation/extrapolation, returns best-fit y for a specified x (e.g. for Mr)
    intercept, slope = linefit1(n, xarr, yarr)
    ycal = intercept + slope * xval

    return ycal


def linefit1(n, xarr, yarr):
    sum_x = sum(xarr)
    sum_y = sum(yarr)
    sum_xy = sum(x * y for x, y in zip(xarr, yarr))
    sum_x_squared = sum(x ** 2 for x in xarr)

    slope = (n * sum_xy - sum_x * sum_y) / (n * sum_x_squared - sum_x ** 2)
    intercept = (sum_y - slope * sum_x) / n

    return intercept, slope

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
        while not (loop_fields[j] - B_offset <= grid_fields[i]) or (j >= n_loop // 2):
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

        ycal = y_polate(2, xvals, yvals, grid_fields[i], grid_moments[i])
        grid_moments[i]=ycal
        #print(grid_moments[i])
    j -= 1
    grid_moments[n_loop // 2] = loop_moments[n_loop // 2] - M_offset

    for i in range(n_loop // 2, n_loop):
        while not (loop_fields[j] - B_offset <= grid_fields[i]): # or (j <= n_loop // 2 + 1):
            if j <= n_loop // 2 + 1:
                break
            if j > n_loop // 2 + 1:
                j -= 1
        while not (loop_fields[j] - B_offset >= grid_fields[i]) or (j >= n_loop):
            if j < n_loop:
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

        ycal = y_polate(2, xvals, yvals, grid_fields[i], grid_moments[i])
        grid_moments[i] = ycal
        #print(grid_moments[i])
    j = 2

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

    return grid_fields, grid_moments






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

# define function for splitting the curves into warm and cool cycles
def split_warm_cool(experiment):
    '''
    splits the X-T curve into warm and cool cycles

    Parameters
    ----------
    experiment : pandas DataFrame
        the IRM experiment data exported into MagIC format

    Returns
    -------
    warm_T : list
        list of temperatures for the warm cycle
    warm_X : list
        list of susceptibilities for the warm cycle
    cool_T : list
        list of temperatures for the cool cycle
    '''
    Tlist = experiment['meas_temp'] # temperature list
    Xlist = experiment['susc_chi_mass'] # Chi list
    
    warmorcool = np.array(np.insert((np.diff(Tlist) > 0 )* 1, 0, 1))
#     print(warmorcool)
    warm_T = [Tlist[i] for i in range(len(warmorcool)) if warmorcool[i]==1]
    cool_T = [Tlist[i] for i in range(len(warmorcool)) if warmorcool[i]==0]
    warm_X = [Xlist[i] for i in range(len(warmorcool)) if warmorcool[i]==1]
    cool_X = [Xlist[i] for i in range(len(warmorcool)) if warmorcool[i]==0]

    return warm_T, warm_X, cool_T, cool_X

# define function for plotting the X-T curve
def plot_X_T(experiment, 
             temp_unit='C', 
             smooth_window=0,
             remove_holder=True, 
             min_temp=None, 
             max_temp=None):
    '''
    plot the high temperature X-T curve

    Parameters
    ----------
    experiment : pandas DataFrame
        the IRM experiment data exported into MagIC format
    temp_unit : str
        the unit of temperature, either 'K' or 'C'
    smooth_window : int
        the window size for smoothing the data
    remove_holder : bool
        whether to remove the holder signal
    min_temp : float
        the minimum temperature for the holder signal
    max_temp : float
        the maximum temperature for the holder signal

    Returns
    -------
    fig : plotly.graph_objs.Figure
        the plotly figure object
    '''
    warm_T, warm_X, cool_T, cool_X = split_warm_cool(experiment)

    # Create a plot for the 'heat' dataset with red points and a line
    fig_heat = go.Figure()
    if temp_unit == 'C':
        warm_T = [T-273.15 for T in warm_T]
        cool_T = [T-273.15 for T in cool_T]

    else:
        raise ValueError('temp_unit must be either "K" or "C"')
    
    if remove_holder:
        assert min_temp is not None, 'min_temp must be provided'
        assert min_temp < max_temp, 'min_temp must be less than max_temp'
        assert min_temp < max(warm_T), 'min_temp must be less than the maximum temperature in the warm cycle'
        assert max_temp is not None, 'max_temp must be provided'
        # now use the min max temp range to select the holder X data
        holder_warm_X = [X for X, T in zip(warm_X, warm_T) if T>=min_temp and T<=max_temp]
        holder_cool_X = [X for X, T in zip(cool_X, cool_T) if T>=min_temp and T<=max_temp]
        holder_warm_X_average = np.mean(holder_warm_X)
        holder_cool_X_average = np.mean(holder_cool_X)
        warm_X = [X - holder_warm_X_average for X in warm_X]
        cool_X = [X - holder_cool_X_average for X in cool_X]

    smoothed_warm_T, smoothed_warm_X, smoothed_warm_Tvars, smoothed_warm_Xvars = X_T_running_average(warm_T, warm_X, smooth_window)
    smoothed_cool_T, smoothed_cool_X, smoothed_cool_Tvars, smoothed_cool_Xvars = X_T_running_average(cool_T, cool_X, smooth_window)

    # Add heat data as a scatter plot with red points and lines
    fig_heat.add_trace(go.Scatter(
        x=warm_T,
        y=warm_X,
        mode='markers',
        marker=dict(color='red', opacity=0.5),
        name='Heating - zero corrected'
    ))
    fig_heat.add_trace(go.Scatter(
        x=smoothed_warm_T,
        y=smoothed_warm_X,
        mode='lines',
        line=dict(color='red'),
        name='Heating - zero corrected - smoothed'
    ))
    # Create a scatter plot for the 'cooling' dataset with blue points and a line
    fig_cool = go.Figure()

    # Add cooling data as a scatter plot with blue points and lines
    fig_cool.add_trace(go.Scatter(
        x=cool_T,
        y=cool_X,
        mode='markers',
        marker=dict(color='blue', opacity=0.5),
        name='Cooling - zero corrected'
    ))
    fig_cool.add_trace(go.Scatter(
        x=smoothed_cool_T,
        y=smoothed_cool_X,
        mode='lines',
        line=dict(color='blue'),
        name='Cooling - zero corrected - smoothed'
    ))

    # Combine the two figures
    fig = go.Figure(data=fig_heat.data + fig_cool.data)

    width = 900  # Adjust the width to your preference
    height = int(width / 1.618)  # Calculate height based on the golden ratio

    # Update the layout with a unified title and labels
    fig.update_layout(
            title={
            'text': f"{experiment['specimen'].unique()[0]}",
            'x': 0.5,  # Center the title
            'xanchor': 'center'
        },
        width=width,
        height=height,
        xaxis_title=f'Temperature (&deg;{temp_unit})',
        yaxis_title='<i>k</i> (m<sup>3</sup> kg<sup>-1</sup>)',
        showlegend=True,
        paper_bgcolor='white',  # Background color of the entire plot
        plot_bgcolor='white',   # Background color of the plotting area
        font=dict(
            family='Roboto, sans-serif',
            size=18,
            color="Black"
        ),
        xaxis=dict(
            tick0=0,  # Start at 0
            dtick=100,  # Tick every 100 units
            gridcolor='lightgray',  # Color of the grid lines
            gridwidth=1,  # Width of the grid lines
            showline=True,  # Show x-axis line
            linewidth=1,  # Width of the x-axis line
            linecolor='black'  # Color of the x-axis line
        ),
        yaxis=dict(
            # Automatic tick marks on the y-axis
            autorange= True,  # Reversed for easier reading
            tickformat='.1e',  # Scientific notation format
            gridcolor='lightgray',  # Color of the grid lines
            gridwidth=1,  # Width of the grid lines
            showline=True,  # Show y-axis line
            linewidth=1,  # Width of the y-axis line
            linecolor='black'  # Color of the y-axis line
        ),
            shapes=[
            # Frame around the plot
            dict(
                type='rect',
                xref='paper', yref='paper',
                x0=0, y0=0, x1=1, y1=1,
                line=dict(color='black', width=2)
            )
        ],
        margin=dict(
            l=50,  # Adjust left margin
            r=50,  # Adjust right margin
            b=50,  # Adjust bottom margin
            t=80   # Adjust top margin for title
        )
    )

    # add a new figure below the X-T plot

    fig_dxdt = go.Figure()
    dxdt = np.gradient(smoothed_warm_X, smoothed_warm_T)
    fig_dxdt.add_trace(go.Scatter(
        x=smoothed_warm_T,
        y=dxdt,
        mode='markers+lines',
        line=dict(color='red'),
        name='Heating - dX/dT smoothed'
    ))

    dxdt = np.gradient(smoothed_cool_X, smoothed_cool_T)
    fig_dxdt.add_trace(go.Scatter(
        x=smoothed_cool_T,
        y=dxdt,
        mode='markers+lines',
        line=dict(color='blue'),
        name='Cooling - dX/dT smoothed'
    ))

    fig_dxdt.update_layout(
        title={
            'text': f"{experiment['specimen'].unique()[0]} - dX/dT",
            'x': 0.5,  # Center the title
            'xanchor': 'center'
        },
        width=width,
        height=height,
        xaxis_title=f'Temperature (&deg;{temp_unit})',
        yaxis_title='dX/dT',
        showlegend=True,
        paper_bgcolor='white',  # Background color of the entire plot
        plot_bgcolor='white',   # Background color of the plotting area
        font=dict(
            family='Roboto, sans-serif',
            size=18,
            color="Black"
        ),
        xaxis=dict(
            tick0=0,  # Start at 0
            dtick=100,  # Tick every 100 units
            gridcolor='lightgray',  # Color of the grid lines
            gridwidth=1,  # Width of the grid lines
            showline=True,  # Show x-axis line
            linewidth=1,  # Width of the x-axis line
            linecolor='black'  # Color of the x-axis line
        ),
        yaxis=dict(
            # Automatic tick marks on the y-axis
            autorange=True,  # Reversed for easier reading
            tickformat='.1e',  # Scientific notation format
            gridcolor='lightgray',  # Color of the grid lines
            gridwidth=1,  # Width of the grid lines
            showline=True,  # Show y-axis line
            linewidth=1,  # Width of the y-axis line
            linecolor='black'  # Color of the y-axis line
        ),
        shapes=[
            # Frame around the plot
            dict(
                type='rect',
                xref='paper', yref='paper',
                x0=0, y0=0, x1=1, y1=1,
                line=dict(color='black', width=2)
            )
        ],
        margin=dict(
            l=50,  # Adjust left margin
            r=50,  # Adjust right margin
            b=50,  # Adjust bottom margin
            t=80   # Adjust top margin for title
        )
    )

    # add a new figure below the X-T plot
    fig_inv = go.Figure()
    inv_warm_X = [1/X for X in smoothed_warm_X]
    inv_cool_X = [1/X for X in smoothed_cool_X]
    fig_inv.add_trace(go.Scatter(
        x=smoothed_warm_T,
        y=inv_warm_X,
        mode='markers+lines',
        line=dict(color='red'),
        name='Heating - 1/X smoothed'
    ))
    fig_inv.add_trace(go.Scatter(
        x=smoothed_cool_T,
        y=inv_cool_X,
        mode='markers+lines',
        line=dict(color='blue'),
        name='Cooling - 1/X smoothed'
    ))

    fig_inv.update_layout(
        title={
            'text': f"{experiment['specimen'].unique()[0]} - 1/X",
            'x': 0.5,  # Center the title
            'xanchor': 'center'
        },
        width=width,
        height=height,
        xaxis_title=f'Temperature (&deg;{temp_unit})',
        yaxis_title='1/X',
        showlegend=True,
        paper_bgcolor='white',  # Background color of the entire plot
        plot_bgcolor='white',   # Background color of the plotting area
        font=dict(
            family='Roboto, sans-serif',
            size=18,
            color="Black"
        ),
        xaxis=dict(
            tick0=0,  # Start at 0
            dtick=100,  # Tick every 100 units
            gridcolor='lightgray',  # Color of the grid lines
            gridwidth=1,  # Width of the grid lines
            showline=True,  # Show x-axis line
            linewidth=1,  # Width of the x-axis line
            linecolor='black'  # Color of the x-axis line
        ),
        yaxis=dict(
            # Automatic tick marks on the y-axis
            autorange=True,  # Reversed for easier reading
            tickformat='.1e',  # Scientific notation format
            gridcolor='lightgray',  # Color of the grid lines
            gridwidth=1,  # Width of the grid lines
            showline=True,  # Show y-axis line
            linewidth=1,  # Width of the y-axis line
            linecolor='black'  # Color of the y-axis line
        ),
        shapes=[
            # Frame around the plot
            dict(
                type='rect',
                xref='paper', yref='paper',
                x0=0, y0=0, x1=1, y1=1,
                line=dict(color='black', width=2)
            )
        ],
        margin=dict(
            l=50,  # Adjust left margin
            r=50,  # Adjust right margin
            b=50,  # Adjust bottom margin
            t=80   # Adjust top margin for title
        )
    )


    # display the main figure
    fig.show()

    # display and return the other two figures depending on the plot_dxdt and plot_inv
    fig_dxdt.show()
    fig_inv.show()

    return fig, fig_dxdt, fig_inv


def X_T_running_average(temp_list, chi_list, temp_window):
    if not temp_list or not chi_list or temp_window <= 0:
        return temp_list, chi_list, [], []
    
    avg_temps = []
    avg_chis = []
    temp_vars = []
    chi_vars = []
    n = len(temp_list)
    
    for i in range(n):
        # Determine the temperature range for the current point
        temp_center = temp_list[i]
        start_temp = temp_center - temp_window / 2
        end_temp = temp_center + temp_window / 2
        
        # Get the indices within the temperature range
        indices = [j for j, t in enumerate(temp_list) if start_temp <= t <= end_temp]
        
        # Calculate the average temperature and susceptibility for the current window
        if indices:
            temp_range = [temp_list[j] for j in indices]
            chi_range = [chi_list[j] for j in indices]
            avg_temp = sum(temp_range) / len(temp_range)
            avg_chi = sum(chi_range) / len(chi_range)
            temp_var = np.var(temp_range)
            chi_var = np.var(chi_range)
        else:
            avg_temp = temp_center
            avg_chi = chi_list[i]
            temp_var = 0
            chi_var = 0
        
        avg_temps.append(avg_temp)
        avg_chis.append(avg_chi)
        temp_vars.append(temp_var)
        chi_vars.append(chi_var)
    
    return avg_temps, avg_chis, temp_vars, chi_vars

def optimize_X_T_running_average_window(experiment, min_temp_window=0, max_temp_window=50, steps=50):
    warm_T, warm_X, cool_T, cool_X = split_warm_cool(experiment)
    windows = np.linspace(min_temp_window, max_temp_window, steps)
    fig, ax = plt.subplots(figsize=(10, 6))

    # Normalize the colormap
    norm = colors.Normalize(vmin=min_temp_window, vmax=max_temp_window)

    for window in windows:
        _, warm_avg_chis, _, warm_chi_vars = X_T_running_average(warm_T, warm_X, window)
        warm_avg_rms, warm_avg_variance = calculate_avg_variance_and_rms(warm_X, warm_avg_chis, warm_chi_vars)
        _, cool_avg_chis, _, cool_chi_vars = X_T_running_average(cool_T, cool_X, window)  
        cool_avg_rms, cool_avg_variance = calculate_avg_variance_and_rms(cool_X, cool_avg_chis, cool_chi_vars)

        ax.scatter(warm_avg_variance, warm_avg_rms, c=window, cmap='Reds', norm=norm)
        ax.scatter(cool_avg_variance, cool_avg_rms, c=window, cmap='Blues', norm=norm)
        # ax.text(warm_avg_variance, warm_avg_rms, f'{window:.2f}°C', fontsize=12, ha='right')
        # ax.text(cool_avg_variance, cool_avg_rms, f'{window:.2f}°C', fontsize=12, ha='right')
    ax.set_xlabel('Average Variance', fontsize=14)
    ax.set_ylabel('Average RMS', fontsize=14)
    ax.set_title('Average RMS vs Average Variance plot\nfor optimizing the running average window size')
    ax.invert_yaxis()
    # show the colormaps and make sure the range is correct
    warm_cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='Reds', norm=norm))
    warm_cbar.set_label('Warm cycle window size (°C)')
    cool_cbar = plt.colorbar(plt.cm.ScalarMappable(cmap='Blues', norm=norm))
    cool_cbar.set_label('Cool cycle window size (°C)')

    return fig, ax   


def calculate_avg_variance_and_rms(chi_list, avg_chis, chi_vars):
    rms_list = np.sqrt([(chi - avg_chi)**2 for chi, avg_chi in zip(chi_list, avg_chis)])
    total_rms = np.sum(rms_list)
    avg_rms = total_rms / len(rms_list)
    
    total_variance = np.sum(chi_vars)
    avg_variance = total_variance / len(chi_vars)
    
    return avg_rms, avg_variance

def plot_MPMS_AC_X_T(experiment, frequency=None, phase='in', figsize=(6,6)):
    """
    This function plots the AC susceptibility data from the MPMS-X for a given experiment
    
    Parameters
    ----------
    experiment : pandas DataFrame
        The experiment table from the MagIC contribution
    frequency : float
        The frequency of the AC measurement in Hz
    phase : str
        The phase of the AC measurement ('in' or 'out' or 'both')
    """
    assert phase in ['in', 'out', 'both'], 'phase should be either "in" or "out" or "both"'
    assert frequency is None or frequency in experiment['frequency'].unique(), 'frequency should be one of the available frequencies'

    if frequency is None:
        # that means we will plot all frequencies
        meas_freqs = experiment['meas_freq'].unique()
    else:
        meas_freqs = [frequency]

    if phase != 'both':
        
        fig, ax = plt.subplots(figsize=figsize)
        for meas_freq in meas_freqs:
            data = experiment[(experiment['meas_freq']==meas_freq)]
            # plot the data
            if phase == 'in':
                ax.plot(data['meas_temp'], data['susc_chi_mass'], 'o-', label=f'{meas_freq} Hz')
            else:
                ax.plot(data['meas_temp'], data['susc_chi_qdr_mass'], 'o-', label=f'{meas_freq} Hz')

        ax.set_xlabel('Temperature (K)', fontsize=16)
        ax.set_ylabel('$m^3/kg$', fontsize=16)
        ax.set_title('AC Susceptibility '+phase+' phase', fontsize=16)
        ax.legend()
        return fig, ax
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        for meas_freq in meas_freqs:
            data = experiment[(experiment['meas_freq']==meas_freq)]
            # plot the data
            ax1.plot(data['meas_temp'], data['susc_chi_mass'], 'o-', label=f'{meas_freq} Hz')
            ax2.plot(data['meas_temp'], data['susc_chi_qdr_mass'], 'o-', label=f'{meas_freq} Hz')
        ax1.set_xlabel('Temperature (K)', fontsize=16)
        ax1.set_ylabel('$m^3/kg$', fontsize=16)
        ax1.set_title('AC Susceptibility in phase', fontsize=16)
        ax1.legend()
        ax2.set_xlabel('Temperature (K)', fontsize=16)
        ax2.set_ylabel('$m^3/kg$', fontsize=16)
        ax2.set_title('AC Susceptibility out phase', fontsize=16)
        ax2.legend()
        return fig, (ax1, ax2)

def backfield_data_processing(experiment, smooth_frac=0.0, drop_first=False):
    '''
    Function to process the backfield data
    including shifting the magnetic moment to be positive values
              taking the log base 10 of the magnetic field values
    and writing these new fields into the experiment attribute table

    Parameters
    ----------
    experiment : DataFrame
        DataFrame containing the backfield data
    smooth_frac : float
        Fraction of the data to be used for smoothing, value must be between 0 and 1
    drop_first : bool
        Whether to drop the first data point or not
        in some cases you may want to drop the first data point to avoid negative log values
    
    Returns
    -------
    experiment dataframe with the new attributes
    '''
    assert smooth_frac >= 0 and smooth_frac <= 1, 'smooth_frac must be between 0 and 1'
    assert isinstance(drop_first, bool), 'drop_first must be a boolean'
    # check and make sure to force drop first row if the first treat field is in the wrong direction
    if experiment['treat_dc_field'].iloc[0] > 0:
        drop_first = True
    if drop_first:
        experiment = experiment.iloc[1:].reset_index(drop=1)
    
    # to plot the backfield data in the conventional way, we need to shift the magnetization to be positive
    experiment['magn_mass_shift'] = [i - experiment['magn_mass'].min() for i in experiment['magn_mass']]
    # we then calculate the log10 of the treatment fields
    experiment['log_dc_field'] = np.log10(-experiment['treat_dc_field']*1e3)
    # loess smoothing
    spl = lowess(experiment['magn_mass_shift'], experiment['log_dc_field'], frac=smooth_frac)
    experiment['smoothed_magn_mass_shift'] = spl[:, 1]
    experiment['smoothed_log_dc_field'] = spl[:, 0]
    return experiment
    
def plot_backfield_data(experiment, figsize=(5, 12)):
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=figsize)

    # raw data
    ax1.scatter(experiment['treat_dc_field'], experiment['magn_mass'], c='black', marker='o', s=10, label='raw backfield data')
    ax1.plot(experiment['treat_dc_field'], experiment['magn_mass'], c='black')
    ax1.set_xlabel('treatment field (T)', fontsize=14)
    ax1.set_ylabel('magnetization (Am$^2$/kg)', fontsize=14)
    ax1.set_title('raw backfield data')
    ax1.legend()
    # processed data as curve on top of raw scattered data points in log10 space
    ax2.scatter(experiment['log_dc_field'], experiment['magn_mass_shift'], c='grey', marker='o', s=10, label='shifted raw data')
    ax2.plot(experiment['smoothed_log_dc_field'], experiment['smoothed_magn_mass_shift'], c='k', label='smoothed shifted raw data')
    ax2ticks = ax2.get_xticks()
    ax2.set_xticklabels([f'{round(10**i, 1)}' for i in ax2ticks])
    ax2.set_xlabel('treatment field (mT)', fontsize=14)
    ax2.set_ylabel('magnetization (Am$^2$/kg)', fontsize=14)
    ax2.set_title('processed backfield data')
    ax2.legend()
    # coercivity spectrum
    # first show the raw data in the derivative space
    raw_derivatives_y = -np.diff(experiment['magn_mass_shift'])/np.diff(experiment['log_dc_field'])
    # take the middle points of the logB values, and also get rid of the nan values
    raw_derivatives_x = experiment['log_dc_field'].rolling(window=2).mean().dropna()
    # derivatives of smoothed data
    smoothed_derivatives_y = -np.diff(experiment['smoothed_magn_mass_shift'])/np.diff(experiment['smoothed_log_dc_field'])
    smoothed_derivatives_x = experiment['smoothed_log_dc_field'].rolling(window=2).mean().dropna()
    ax3.scatter(raw_derivatives_x, raw_derivatives_y, c='grey', marker='o', s=10, label='raw coercivity spectrum')
    ax3.plot(smoothed_derivatives_x, smoothed_derivatives_y, c='k', label='smoothed coercivity spectrum')
    ax3ticks = ax3.get_xticks()
    ax3.set_xticklabels([f'{round(10**i, 1)}' for i in ax3ticks])
    ax3.set_xlabel('treatment field (mT)', fontsize=14)
    ax3.set_ylabel('dM/dB', fontsize=14)
    ax3.set_title('coercivity spectrum')
    ax3.legend()
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.1, hspace=0.4, wspace=0.2)

    return fig, (ax1, ax2, ax3)

def backfield_unmixing(experiment, n_comps=1, parameters=None):
    '''
    backfield unmixing for a single experiment

    Parameters
    ----------
    experiment : DataFrame
        DataFrame containing the backfield data
        by default we unmix on the smothed data
    n_comps : int
        Number of components to unmix the data into
    params : Pandas DataFrame
        Initial values for the model parameters
        should be constructed as the following columns:
        - amplitude in arbitrary scale
        - center in unit of mT
        - sigma in unit of mT
        - gamma in arbitrary scale
        |amplitude|center|sigma|gamma|
        |---|---|---|---|
        |1.0|100|10|0.0|
        |...|...|...|...|
        the program will automatically go through the rows and extract these inital parameter values
        If the parameters are not given, we will run an automated program to make initial guess
        
    Returns
    '''

    assert n_comps > 0, 'n_component must be greater than 0'
    assert isinstance(n_comps, int), 'n_component must be an integer'
    assert isinstance(parameters, pd.DataFrame), f"Expected a pandas DataFrame, but got {type(variable).__name__}"
    assert n_comps == parameters.shape[0], 'number of components must match the number of rows in the parameters table'
    # re-calculate the derivatives based on the smoothed data columns
    smoothed_derivatives_y = -np.diff(experiment['smoothed_magn_mass_shift'])/np.diff(experiment['smoothed_log_dc_field'])
    smoothed_derivatives_x = experiment['smoothed_log_dc_field'].rolling(window=2).mean().dropna()

    # create the model depending on the number of components specified
    composite_model = None
    params = Parameters()
    for i in range(n_comps):
        prefix = f'g{i+1}_'
        model = SkewedGaussianModel(prefix=prefix)
        
        # Initial parameter guesses
        params.add(f'{prefix}amplitude', value=parameters['amplitude'][i])
        params.add(f'{prefix}center', value=np.log10(parameters['center'][i]))
        params.add(f'{prefix}sigma', value=np.log10(parameters['sigma'][i]))
        params.add(f'{prefix}gamma', value=parameters['gamma'][i])
        
        # now let's set bounds to the parameters to help fitting algorithm converge
        params[f'{prefix}amplitude'].min = 0  # Bounds for amplitude parameters
        params[f'{prefix}amplitude'].max = np.max(smoothed_derivatives_y)
        params[f'{prefix}center'].min = experiment['smoothed_log_dc_field'].min()  # Bounds for center parameters
        params[f'{prefix}center'].max = experiment['smoothed_log_dc_field'].max()  # Bounds for center parameters
        params[f'{prefix}sigma'].min = 0
        params[f'{prefix}sigma'].max = experiment['smoothed_log_dc_field'].max()-experiment['smoothed_log_dc_field'].min()  # Bounds for sigma parameters

        if composite_model is None:
            composite_model = model
        else:
            composite_model += model

    result = composite_model.fit(smoothed_derivatives_y, params, x=smoothed_derivatives_x)

    # now use the resultant fit to update the parameters table and return it
    for i in range(n_comps):
        prefix = f'g{i+1}_'
        parameters.loc[i, 'amplitude'] = result.params[f'{prefix}amplitude'].value
        parameters.loc[i, 'center'] = 10**result.params[f'{prefix}center'].value # convert back to mT
        parameters.loc[i, 'sigma'] = 10**result.params[f'{prefix}sigma'].value # convert back to mT
        parameters.loc[i, 'gamma'] = result.params[f'{prefix}gamma'].value

    return result, parameters

def plot_backfield_unmixing_result(experiment, result, sigma=2, figsize=(8,6)):
    raw_derivatives_y = -np.diff(experiment['magn_mass_shift'])/np.diff(experiment['log_dc_field'])
    raw_derivatives_x = experiment['log_dc_field'].rolling(window=2).mean().dropna()
    smoothed_derivatives_x = experiment['smoothed_log_dc_field'].rolling(window=2).mean().dropna()

    comps = result.eval_components(x=smoothed_derivatives_x)
    dely = result.eval_uncertainty(sigma=sigma)

    fig, ax = plt.subplots(figsize=figsize)
    # first plot the scatter raw dMdB data
    ax.scatter(raw_derivatives_x, raw_derivatives_y, c='grey', marker='o', s=10, label='raw coercivity spectrum')
    # plot the total best fit
    ax.plot(raw_derivatives_x, result.best_fit, '-', label='total spectrum best fit')
    ax.fill_between(smoothed_derivatives_x,
                            result.best_fit-dely,
                            result.best_fit+dely,
                            color="#8A8A8A", 
                            label=f'total {sigma}-$\sigma$ band', alpha=0.5)
    if len(result.components) > 1:
        for i in range(len(result.components)):

            ax.plot(smoothed_derivatives_x, comps[f'g{i+1}_'], c='C0', label=f'component #{i+1}, {sigma}-$\sigma$ band')
            ax.fill_between(smoothed_derivatives_x,
                                    comps[f'g{i+1}_']-result.dely_comps[f'g{i+1}_'],
                                    comps[f'g{i+1}_']+result.dely_comps[f'g{i+1}_'],
                                    color=f'C{i}', alpha=0.5)
    xticks = ax.get_xticks()
    ax.set_xticklabels([f'{int(10**i)}' for i in xticks])
    ax.legend()
    ax.set_title('coercivity unmixing results')
    ax.set_xlabel('treatment field (mT)', fontsize=14)
    ax.set_ylabel('dM/dB', fontsize=14)
    return fig, ax
	
def kDTA_interp(blank_time, blank_temp, sample_time, sample_temp):
	'''
	Function to calculate residual temperature from sample and blank time and temperature data.
	Inputs:
	Outputs:
	'''
	#interpolation
	t_max = max(max(blank_time), max(sample_time))### change to accomodate temp max
	
	t = np.arange(0,t_max + 15,15)
	int_blank = interp1d(blank_time[:-5],blank_temp[:-5], kind = 'linear', bounds_error = False)
	temp_blank = int_blank(t)
	
	int_sample = interp1d(sample_time[:-5],sample_temp[:-5], kind = 'linear', bounds_error = False)
	temp_sample = int_sample(t)
	
	return t, temp_blank, temp_sample