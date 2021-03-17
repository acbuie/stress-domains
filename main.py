"""
Plot a deformation mechanism map of quartz. 
X-Axis: Homologous temperatures, from 0.2 to 1.0.
Y-Axis: Log10(Theta/Mu) (normalized sheer stress) from 0 to -6.
Z-Axis/Contour: Strain rate from 10^-6 to 10^-15.
"""
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
from matplotlib import cm

import numpy as np
from numpy import typing
import pandas as pd

from colour import Color

from typing import List

# Import constants
from constants import *

def dislocation_creep(
    temperature: np.ndarray, n_shear_stress: np.ndarray) -> np.ndarray:
    """
    Constitutive equation for dislocation creep.

    Returns:
    - Strain rate.
    """
    numerator = (
        MU * B * D_L * np.power(n_shear_stress, 3) * 
        np.exp(-H_L / (R * temperature))
    )

    denomenator = K * temperature
    
    strain_rate = numerator/denomenator

    # Note, must use np.exp and np.power instead of math module variants
    # to apply across an ndarray

    return strain_rate

def nabarro_herring_creep(
    temperature: np.ndarray, n_shear_stress: np.ndarray) -> np.ndarray:
    """
    Constitutive equation for Nabarro-Herring creep.

    Returns:
    - Strain rate.
    """
    numerator = (
        A_NH * MU * V * D_L * n_shear_stress * 
        np.exp(-H_L / (R * temperature))
    )
    denomenator = R * D ** 2 * temperature
    
    strain_rate = numerator/denomenator

    return strain_rate

def coble_creep(
    temperature: np.ndarray, n_shear_stress: np.ndarray) -> np.ndarray:
    """
    Constitutive equation for Coble creep.

    Returns:
    - Strain rate.
    """
    numerator = (
        A_C * MU * V * D_G * W * n_shear_stress *
        np.exp(-H_G / (R * temperature))
    )
    denomenator =  R * D ** 3 * temperature
    
    strain_rate = numerator/denomenator

    return strain_rate

def three_array_max(array_list: List[np.ndarray]) -> np.ndarray:
    """
    Finds the maximum value between three np.ndarrays in a list.
    """
    temp = np.maximum(array_list[0], array_list[1])
    all_maxs = np.maximum(temp, array_list[2])

    return all_maxs

def find_intersections() -> np.ndarray:
    """
    Finds the points of intersection between the arrays.

    Will need to work on this, could be tricky.
    """
    pass

def make_plot(X: np.ndarray, Y: np.ndarray, Z: np.ndarray, 
    levels: np.ndarray, color_range: List[str]):
    """
    Create deformation mechanism map.

    X and Y inputs are xx and yy from np.meshgrid(x, y)
    Z input is any value calculated across xx and yy. 
    Contours will be at the levels specified, and colors will be 
    selected from two colors in color_range.
    """

    color1 = Color(color_range[0])
    color2 = Color(color_range[1])
    
    # List of Color objects, need to convert for matplotlib
    color_list = list(color1.range_to(color2, len(levels)))
    
    contour_colors = []
    for color in color_list:
        contour_colors.append(color.hex)

    fig = plt.figure()
    ax1 = fig.add_subplot(111)

    fig.suptitle('Quartz Deformation Mechanism Map')
    
    CS = ax1.contour(X, Y, Z, levels, colors = contour_colors) 
    
    # Labels... hardcoded. May find a way to make this generic
    sr_labels = [
        r'$10^{-15}$', r'$10^{-14}$', r'$10^{-13}$', 
        r'$10^{-12}$', r'$10^{-11}$', r'$10^{-10}$', 
        r'$10^{-9}$', r'$10^{-8}$', r'$10^{-7}$', 
        r'$10^{-6}$']

    for i in range(len(sr_labels)):
        CS.collections[i].set_label(sr_labels[i])

    # Reverse legend order, set location and title
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(reversed(handles), reversed(labels), 
        loc = 'upper right', title = r'$\.\epsilon$')

    # Set log scale and ylim
    ax1.set_yscale('log')
    ax1.set_ylim([10 ** -6, 1])
   
    # Set labels... could be done better in Illustrator...
    ax1.set_ylabel(
        r'Shear Stress $log\left( \frac{\sigma}{\mu}\right)$'
    )
    ax1.set_xlabel(
        r'Homologous Temperature $\left( \frac{T}{T_m} \right)$'
    )

    plt.show()

if __name__ == "__main__":
    # Using np.meshgrid for speed
    x_values = np.linspace(0.2, 1.0, 1000)
    y_values = np.geomspace(1, 10 ** -6, 1000)

    # Specify contour levels
    strain_levels = np.geomspace(10 ** -15, 10 ** -6, 10)

    # Make coordinate values
    # Homologous T (xx), Normalized Shear Stress (y)
    ht_x, n_ss_y = np.meshgrid(x_values, y_values)
    
    # Find strain rates for DC, NH, and CC
    sr_dc = dislocation_creep(ht_x * T_M, n_ss_y)
    sr_nh = nabarro_herring_creep(ht_x * T_M, n_ss_y)
    sr_cc = coble_creep(ht_x * T_M, n_ss_y)

    # Take highest surface (z) values
    max_strain_rates = three_array_max([sr_dc, sr_nh, sr_cc])

    # Set plotting color range
    colors = ['blue', 'yellow']

    # Make the plot
    make_plot(ht_x, n_ss_y, max_strain_rates, strain_levels, colors)
