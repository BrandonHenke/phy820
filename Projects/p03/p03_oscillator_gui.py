import ipywidgets as widgets
from ipywidgets import HBox, VBox, Layout, Tab, Label, Checkbox, Button
from ipywidgets import FloatSlider, IntSlider, Play, Dropdown, HTMLMath 
from ipywidgets import fixed

import numpy as np
import uuid

from time import sleep



def plot_choice_widget(on=True, plot_description=None):
    """
    	Makes a Checkbox to select whether to show a plot.
    """
    return Checkbox(value=on, description=plot_description,
                  disabled=False, indent=False, layout=Layout(width='150px'))


# Widgets for the oscillator parameters (all use FloatSlider, so we made 
#  it a function)
def float_widget(value, min, max, step, description, format):
    """
    	Makes a FloatSlider with the passed parameters and continuous_update
       	set to False.
	"""
    slider_border = Layout(border='solid 1.0px')
    return FloatSlider(value=value,min=min,max=max,step=step,disabled=False,
                       description=description,continuous_update=False,
                       orientation='horizontal',layout=slider_border,
                       readout=True,readout_format=format)


class OscillatorGUI():

	def __init__(self, **kwargs):

		# Widgets for the plot choice (plus a label out front)
		self.plot_choice_w = Label(value='Which plots: ',layout=Layout(width='100px'))

		# checkboxes
		self.q_vs_time_plot_w = plot_choice_widget(True, r'$q$ vs. time')
		self.p_vs_time_plot_w = plot_choice_widget(True, r'$p$ vs. time')
		self.phase_space_plot_w = plot_choice_widget(True, 'phase space')
		self.driving_curve_w = plot_choice_widget(False, 'driving force')
		self.energy_plot_w = plot_choice_widget(False, 'energy')
		self.poincare_plot_w = plot_choice_widget(True, 'Poincare map')



		# controls for the base oscillator class
		self.omega0_w = float_widget(value=1.0, min=0.0, max=2.*np.pi, step=0.001,
                        description=r'natural $\omega_0$:', format='.4f')
		self.beta_w = float_widget(value=0.05, min=0.0, max=2.*np.pi, step=0.001,
                       description=r'damping $\beta$:', format='.4f')
		self.f_ext_w = float_widget(value=0.4, min=0.0, max=2, step=0.001,
                       description=r'$f_\text{ext}$:', format='.4f')
		self.omega_ext_w = float_widget(value=1.4, min=0.0, max=2.*np.pi, step=0.001,
                        description=r'$\omega_\text{ext}$:', format='.4f')
		self.phi_ext_w = float_widget(value=0, min=0.0, max=2*np.pi, step=0.001,
                        description=r'$\phi_\text{ext}$:', format='.4f')




		# Widgets for the initial conditions
		self.q0_w = float_widget(value=0.0, min=0., max=100, step=0.001,
		                        description=r'$q_0$:', format='.3f')
		self.p0_w = float_widget(value=0.0, min=-100., max=100., step=0.001,
		                            description=r'$p_0$:', format='.3f')

		# Widgets for the plotting parameters
		self.t_start_w = float_widget(value=0., min=0., max=20000., step=1.,
		                         description='t start:', format='.1f') 
		self.t_end_w = float_widget(value=5., min=0., max=20000., step=1.,
		                       description='t end:', format='.1f')
		self.delta_t_w = float_widget(value=0.001, min=0.001, max=0.1, step=0.001,
		                         description='delta t:', format='.3f')
		self.plot_start_w = float_widget(value=0., min=0., max=20000., step=1.,
                            description='start plotting:', format='.1f')


		self.t_start_w.observe(self.update_t_end, 'value')
		self.t_end_w.observe(self.update_t_end, 'value')


		self.plot_start_w.observe(self.update_plot_start, 'value')
		self.t_start_w.observe(self.update_plot_start, 'value')
		self.t_end_w.observe(self.update_plot_start, 'value')


		# Widgets for the styling parameters
		self.font_size_w = Dropdown(options=['12', '16', '18', '20', '24'], value='18',
		                       description='Font size:',disabled=False,
		                       continuous_update=False,layout=Layout(width='140px'))


		for key, value in kwargs.items():
			setattr(self, key, value)

        # try:
    	# Set up the interactive_output widget 
		plot_out = widgets.interactive_output(self.plotfunction,
		                          dict(
		                          q_vs_time_plot=self.q_vs_time_plot_w,
		                          p_vs_time_plot=self.p_vs_time_plot_w,
		                          phase_space_plot=self.phase_space_plot_w,
		                          driving_curve = self.driving_curve_w,
		                          energy_plot = self.energy_plot_w,
		                          poincare_plot = self.poincare_plot_w,
		                          omega0=self.omega0_w,
		                          beta=self.beta_w,
		                          f_ext=self.f_ext_w,
                                  omega_ext=self.omega_ext_w,
		                          phi_ext=self.phi_ext_w,
		                          q0=self.q0_w,
		                          p0=self.p0_w,
		                          t_start=self.t_start_w,
		                          t_end=self.t_end_w, 
		                          delta_t=self.delta_t_w,    
		                          plot_start=self.plot_start_w, 
		                          font_size=self.font_size_w)
		                       )
		# except NoneType:
		# 	print("No plot and solver function has been assigned to the GUI.")


		# Some manual layout...
		hbox1 = HBox([self.plot_choice_w, self.q_vs_time_plot_w, self.p_vs_time_plot_w,
		              self.phase_space_plot_w, self.driving_curve_w, self.energy_plot_w, self.poincare_plot_w]) #  choice of plots to show
		hbox2 = HBox([self.q0_w, self.p0_w]) # initial conditions and damping
		hbox2a = HBox([self.omega0_w, self.beta_w]) # oscillator parameters
		hbox3 = HBox([self.f_ext_w, self.omega_ext_w, self.phi_ext_w]) # driving force
		hbox4 = HBox([self.t_start_w, self.t_end_w, self.delta_t_w, self.plot_start_w]) # time, plot ranges
		hbox5 = HBox([self.font_size_w]) # font size
        

        # We'll set up Tabs to organize the controls.  The Tab contents are declared
		#  as tab0, tab1, ... (probably should make this a list?) and the overall Tab
		#  is called tab (so its children are tab0, tab1, ...).
		tab_height = '40px'  # Fixed minimum height for all tabs. Specify another way?
		tab0 = VBox([hbox2, hbox2a, hbox3], layout=Layout(min_height=tab_height))
		tab1 = VBox([hbox1, hbox4], layout=Layout(min_height=tab_height))
		tab2 = VBox([hbox5], layout=Layout(min_height=tab_height))

		tab = Tab(children=[tab0, tab1, tab2])
		tab.set_title(0, 'Physics')
		tab.set_title(1, 'Plotting')
		tab.set_title(2, 'Styling')

		# Release the Kraken!
		self.vbox = VBox([tab, plot_out])


	# Make sure that t_end is at least t_start + 10
	def update_t_end(self, *args):
	    if self.t_end_w.value < self.t_start_w.value:
	        self.t_end_w.value = self.t_start_w.value + 10     


	# Make sure that plot_start is at least t_start and less than t_end
	def update_plot_start(self, *args):
	    if self.plot_start_w.value < self.t_start_w.value:
	        self.plot_start_w.value = self.t_start_w.value
	    if self.plot_start_w.value > self.t_end_w.value:
	        self.plot_start_w.value = self.t_end_w.value
	
