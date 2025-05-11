#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Complex Unit Circle Visualization

An interactive visualization of complex numbers and their relationship to the unit circle.
Shows the connection between complex numbers, trigonometric functions, and polar form.

Author: Matthias Naumann
Version: 1.0.0
Date: May 11, 2025
License: MIT License (see LICENSE file or https://opensource.org/licenses/MIT)
Repository: https://github.com/MaNafromSaar/complex_unit_circle

This software is free to use for educational purposes.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button, TextBox
from matplotlib.patches import Arc
import matplotlib as mpl
import time  # Import time for throttling

class SimpleComplexPlane:
    def __init__(self):
        # Create the plot
        self.fig, self.ax = plt.subplots(figsize=(9, 8))
        
        # Set the window title to include author name
        self.fig.canvas.manager.set_window_title('Complex Unit Circle Visualization - by Matthias Naumann')
        # Adjust the subplot to leave space at bottom and right for controls
        self.fig.subplots_adjust(bottom=0.20, left=0.12, right=0.75, top=0.95)
        
        # Initial view limits
        self.view_scale = 2.0  # Initial scale: -2 to 2
        self.ax.set_xlim(-self.view_scale, self.view_scale)
        self.ax.set_ylim(-self.view_scale, self.view_scale)
        self.ax.grid(True)
        self.ax.set_aspect('equal')
        self.ax.set_title('Complex Number Visualization')
        self.ax.set_xlabel('Real (Cosine Component)')
        self.ax.set_ylabel('Imaginary (Sine Component)')
        
        # Dark mode state
        self.dark_mode = False
        
        # Initialize background for blitting
        self.background = None
        self.need_full_redraw = True
        
        # Draw unit circle
        theta = np.linspace(0, 2*np.pi, 100)
        self.ax.plot(np.cos(theta), np.sin(theta), 'g--', alpha=0.5)
        
        # Mark the axes
        self.ax.axhline(y=0, color='k', linestyle='-', alpha=0.3)
        self.ax.axvline(x=0, color='k', linestyle='-', alpha=0.3)
        
        # Reduce initial operations
        self.special_angle_labels = []
        self.mark_special_angles()
        
        # Initial complex number
        self.real = 1.0
        self.imag = 1.0
        
        # Calculate magnitude and angle
        self.update_polar()
        
        # Draw the vector and point
        self.vector, = self.ax.plot([0, self.real], [0, self.imag], 'r-', lw=2)
        self.point, = self.ax.plot([self.real], [self.imag], 'ro', ms=8)
        
        # Add projections to axes (coordinate projections)
        self.h_projection, = self.ax.plot([0, self.real], [self.imag, self.imag], 'b--', alpha=0.5)
        self.v_projection, = self.ax.plot([self.real, self.real], [0, self.imag], 'b--', alpha=0.5)
        
        # Add unit circle projections (actual sine and cosine)
        self.unit_point, = self.ax.plot([np.cos(self.angle)], [np.sin(self.angle)], 'go', ms=6)
        self.unit_vector, = self.ax.plot([0, np.cos(self.angle)], [0, np.sin(self.angle)], 'g-', lw=1.5)
        
        # Add projections to axes from unit circle point (sine and cosine)
        self.cos_projection, = self.ax.plot([0, np.cos(self.angle)], [0, 0], 'g-', lw=1.5)
        self.sin_projection, = self.ax.plot([0, 0], [0, np.sin(self.angle)], 'g-', lw=1.5)
        
        # Add specific colored dotted lines for sine and cosine
        self.cos_line, = self.ax.plot([np.cos(self.angle), np.cos(self.angle)], [0, np.sin(self.angle)], 
                                     'r:', lw=2, alpha=0.8)
        self.sin_line, = self.ax.plot([0, np.cos(self.angle)], [np.sin(self.angle), np.sin(self.angle)], 
                                     'm:', lw=2, alpha=0.8)
        
        # Add dotted lines connecting the point to its unit circle projection
        self.unit_connection, = self.ax.plot([self.real, np.cos(self.angle)], 
                                            [self.imag, np.sin(self.angle)], 'k:', alpha=0.5)
        
        # Add projection labels
        self.coord_label = self.ax.text(self.real/2, self.imag + 0.1, f"({self.real:.2f}, {self.imag:.2f})", 
                                       ha='center', va='bottom', color='blue', fontsize=8)
        
        self.cos_label = self.ax.text(np.cos(self.angle)/2, -0.1, f"cos(θ)={np.cos(self.angle):.4f}", 
                                     ha='center', va='top', color='green', fontsize=8)
        
        self.sin_label = self.ax.text(-0.2, np.sin(self.angle)/2, f"sin(θ)={np.sin(self.angle):.4f}", 
                                     ha='right', va='center', color='green', fontsize=8)
        
        # Add angle arc
        self.draw_angle_arc()
        
        # Add trigonometric values display directly under the graph
        self.trig_text = self.ax.text(
            0.5, -0.1, '',  
            transform=self.ax.transAxes,
            ha='center', va='top',
            bbox=dict(facecolor='white', alpha=0.9),
            fontsize=9
        )
        
        # Add text displays below the graph - this will now show exact values
        self.info_text = self.ax.text(
            0.5, -0.146, '', 
            transform=self.ax.transAxes,
            ha='center', va='top',
            bbox=dict(facecolor='white', alpha=0.8)
        )

        # Store references to dynamic artists for efficient updates
        self.dynamic_artists = [
            self.vector, self.point, self.h_projection, self.v_projection,
            self.unit_point, self.unit_vector, self.cos_projection, self.sin_projection,
            self.cos_line, self.sin_line, self.unit_connection,
            self.coord_label, self.cos_label, self.sin_label,
            self.trig_text, self.info_text
        ]

        # Add UI elements in a separate method to reduce init time
        self.setup_ui()
        
        # Initial update
        self.update_display()
    
    def setup_ui(self):
        # Add labels for inputs
        label_real = plt.figtext(0.78, 0.75, "Real part:", fontsize=9, ha='left')
        label_imag = plt.figtext(0.78, 0.64, "Imaginary part:", fontsize=9, ha='left')
        label_mag = plt.figtext(0.78, 0.54, "Magnitude:", fontsize=9, ha='left')
        label_angle = plt.figtext(0.78, 0.43, "Angle (degrees):", fontsize=9, ha='left')
        # Add new label for radians input
        label_radians = plt.figtext(0.78, 0.33, "Angle (radians):", fontsize=9, ha='left')
        
        # Add simple text boxes for input - moved right side of graph
        # For each box: [left_position, bottom_position, width, height]
        axbox_real = plt.axes([0.78, 0.70, 0.15, 0.04])     # Real textbox - moved to right side
        axbox_imag = plt.axes([0.78, 0.59, 0.15, 0.04])     # Imaginary textbox - moved to right side
        axbox_mag = plt.axes([0.78, 0.49, 0.15, 0.04])      # Magnitude textbox - moved to right side
        axbox_angle = plt.axes([0.78, 0.38, 0.15, 0.04])    # Angle textbox - moved to right side
        axbox_radians = plt.axes([0.78, 0.28, 0.15, 0.04])  # New radians textbox
        
        # Create text boxes with higher precision
        self.text_real = TextBox(axbox_real, '', initial=f"{self.real:.4f}")
        self.text_imag = TextBox(axbox_imag, '', initial=f"{self.imag:.4f}")
        self.text_mag = TextBox(axbox_mag, '', initial=f"{self.mag:.4f}")
        self.text_angle = TextBox(axbox_angle, '', initial=f"{self.angle_deg:.4f}")
        self.text_radians = TextBox(axbox_radians, '', initial=f"{self.angle:.4f}")  # New radians textbox
        
        # Connect the callbacks
        self.text_real.on_submit(self.submit_real)
        self.text_imag.on_submit(self.submit_imag)
        self.text_mag.on_submit(self.submit_mag)
        self.text_angle.on_submit(self.submit_angle)
        self.text_radians.on_submit(self.submit_radians)
        
        # Add zoom buttons
        label_zoom = plt.figtext(0.78, 0.23, "Zoom:", fontsize=9, ha='left')
        axbtn_zoom_in = plt.axes([0.83, 0.18, 0.05, 0.04])
        axbtn_zoom_out = plt.axes([0.88, 0.18, 0.05, 0.04])
        self.btn_zoom_in = Button(axbtn_zoom_in, '+')
        self.btn_zoom_out = Button(axbtn_zoom_out, '-')
        self.btn_zoom_in.on_clicked(self.zoom_in)
        self.btn_zoom_out.on_clicked(self.zoom_out)
        
        # Add more angle preset buttons with a header
        label_angles = plt.figtext(0.78, 0.16, "Angle Presets:", fontsize=9, ha='left')
        
        # First row of angle buttons
        axbtn_15deg = plt.axes([0.78, 0.11, 0.05, 0.04])
        axbtn_30deg = plt.axes([0.83, 0.11, 0.05, 0.04])
        axbtn_45deg = plt.axes([0.88, 0.11, 0.05, 0.04])
        
        # Second row of angle buttons
        axbtn_60deg = plt.axes([0.78, 0.06, 0.05, 0.04])
        axbtn_90deg = plt.axes([0.83, 0.06, 0.05, 0.04])
        axbtn_120deg = plt.axes([0.88, 0.06, 0.05, 0.04])
        
        # Create buttons for angle presets
        self.btn_15deg = Button(axbtn_15deg, '15°')
        self.btn_30deg = Button(axbtn_30deg, '30°')
        self.btn_45deg = Button(axbtn_45deg, '45°')
        self.btn_60deg = Button(axbtn_60deg, '60°')
        self.btn_90deg = Button(axbtn_90deg, '90°')
        self.btn_120deg = Button(axbtn_120deg, '120°')
        
        # Connect angle preset button callbacks
        self.btn_15deg.on_clicked(self.set_angle_15)
        self.btn_30deg.on_clicked(self.set_angle_30)
        self.btn_45deg.on_clicked(self.set_angle_45)
        self.btn_60deg.on_clicked(self.set_angle_60)
        self.btn_90deg.on_clicked(self.set_angle_90)
        self.btn_120deg.on_clicked(self.set_angle_120)
        
        # Add 15-degree step angle buttons
        axbtn_180deg = plt.axes([0.78, 0.01, 0.05, 0.04])
        axbtn_270deg = plt.axes([0.83, 0.01, 0.05, 0.04])
        axbtn_360deg = plt.axes([0.88, 0.01, 0.05, 0.04])
        
        self.btn_180deg = Button(axbtn_180deg, '180°')
        self.btn_270deg = Button(axbtn_270deg, '270°')
        self.btn_360deg = Button(axbtn_360deg, '360°')
        
        self.btn_180deg.on_clicked(self.set_angle_180)
        self.btn_270deg.on_clicked(self.set_angle_270)
        self.btn_360deg.on_clicked(self.set_angle_360)
        
        # Add dark mode toggle button
        axbtn_dark_mode = plt.axes([0.78, 0.85, 0.15, 0.04])
        self.btn_dark_mode = Button(axbtn_dark_mode, 'Dark Mode')
        self.btn_dark_mode.on_clicked(self.toggle_dark_mode)
        
        # Connect mouse events
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.drag_active = False
        
        # Throttle for update_display calls during dragging
        self.last_update_time = 0
        self.update_throttle = 0.03  # Limit updates to ~30fps

    def mark_special_angles(self):
        """Mark special angles on the unit circle with labels showing degrees, radians, and e^(iθ) form"""
        # Define special angles in degrees - now including all 15-degree increments
        special_angles_deg = [0, 30, 45, 60, 90, 180, 270]
        
        # Clear any existing angle label references
        if not hasattr(self, 'special_angle_labels'):
            self.special_angle_labels = []
        else:
            for label in self.special_angle_labels:
                if label in self.ax.texts:
                    label.remove()
            self.special_angle_labels = []
        
        self.special_angle_points = []  # Store the points for special angles
        self.special_angle_data = {}  # Store the angle data for hover functionality
        
        for angle_deg in special_angles_deg:
            # Convert to radians
            angle_rad = np.radians(angle_deg)
            
            # Calculate point on unit circle
            x = np.cos(angle_rad)
            y = np.sin(angle_rad)
            
            # Draw point
            point, = self.ax.plot(x, y, 'go', ms=6)
            self.special_angle_points.append(point)
            
            # Get radian representation and complex representation
            radian_repr, complex_repr = self.get_angle_representation(angle_deg)
            
            # Create multi-line label with different representations
            label_text = f"{angle_deg}°\n{radian_repr}\n"
            
            # Add simplified complex representation
            complex_parts = complex_repr.split(" = ")
            if len(complex_parts) > 1:
                label_text += complex_parts[1]
            else:
                label_text += complex_repr
            
            # Position the label slightly away from the point
            offset = 0.15
            label_x = x * 1.1
            label_y = y * 1.1
            
            # Adjust text alignment based on quadrant
            ha = 'center'
            va = 'center'
            
            if x > 0.7:  # Right side
                ha = 'left'
            elif x < -0.7:  # Left side
                ha = 'right'
            
            if y > 0.7:  # Top
                va = 'bottom'
            elif y < -0.7:  # Bottom
                va = 'top'
            
            # Create the label but initially make it invisible
            label = self.ax.text(label_x, label_y, label_text, ha=ha, va=va, fontsize=10, 
                         bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'),
                         visible=False)  # Initially invisible
            
            # Store reference to the label
            self.special_angle_labels.append(label)
            
            # Store point coordinates and label for hover functionality
            self.special_angle_data[(x, y)] = label
        
        # Connect hover event handler
        self.hover_cid = self.fig.canvas.mpl_connect('motion_notify_event', self.on_hover)

    def draw_angle_arc(self):
        """Draw an arc to show the angle"""
        radius = 0.3  # Fixed radius for visibility
        theta = np.linspace(0, self.angle, 50)
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        
        # Check if arc already exists
        if hasattr(self, 'angle_arc'):
            self.angle_arc.set_data(x, y)
        else:
            self.angle_arc, = self.ax.plot(x, y, 'b-', lw=1.5)
        
        # Update angle label at midpoint of arc
        mid_angle = self.angle / 2
        label_x = (radius * 1.2) * np.cos(mid_angle)
        label_y = (radius * 1.2) * np.sin(mid_angle)
        
        if hasattr(self, 'angle_label'):
            self.angle_label.remove()
        
        # Format the angle in degrees
        angle_text = f"{self.angle_deg:.2f}°"
        
        # Get special angle representation if available
        radian_repr, complex_repr = self.get_angle_representation(self.angle_deg)
        
        # Add radian representation to angle text only for special angles
        if "π" in radian_repr:
            angle_text += f"\n{radian_repr}"
        
        # Add negative angle representation only for angles between 180° and 360° 
        # to avoid tautologies
        if 180 < self.angle_deg <= 360:
            negative_angle = self.angle_deg - 360
            angle_text += f"\n{negative_angle:.2f}°"
        
        # Add complex representation for unit circle angles
        if abs(self.mag - 1.0) < 0.01:  # If close to unit circle
            angle_text += f"\n{complex_repr.split(' = ')[0]}"
        
        self.angle_label = self.ax.text(label_x, label_y, angle_text, 
                                        fontsize=8, ha='center', va='center',
                                        bbox=dict(facecolor='white', alpha=0.7))

    def update_polar(self):
        """Update magnitude and angle from real and imaginary parts"""
        self.mag = np.sqrt(self.real**2 + self.imag**2)
        self.angle = np.arctan2(self.imag, self.real)
        self.angle_deg = np.degrees(self.angle)
    
    def update_cartesian(self):
        """Update real and imaginary parts from magnitude and angle"""
        self.real = self.mag * np.cos(self.angle)
        self.imag = self.mag * np.sin(self.angle)
    
    def get_exact_form(self, value):
        """Return the exact form of common trigonometric values"""
        # Common angles with exact forms
        special_values = {
            0: "0",
            1: "1",
            -1: "-1",
            0.5: "1/2",
            -0.5: "-1/2",
            np.sqrt(2)/2: "√2/2",
            -np.sqrt(2)/2: "-√2/2",
            np.sqrt(3)/2: "√3/2",
            -np.sqrt(3)/2: "-√3/2",
            1/np.sqrt(2): "1/√2",
            -1/np.sqrt(2): "-1/√2",
            np.sqrt(3)/3: "√3/3",
            -np.sqrt(3)/3: "-√3/3"
        }
        
        # Check for close matches
        for exact_val, exact_form in special_values.items():
            if abs(value - exact_val) < 1e-13: # Tolerance for floating point comparison
                return exact_form
        
        return f"{value:.4f}"
    
    def update_display(self, force_full_redraw=False):
        """Update the visualization"""
        # Check if throttling is needed during drag operations
        current_time = time.time()
        if self.drag_active and (current_time - self.last_update_time < self.update_throttle):
            return
        self.last_update_time = current_time
        
        # Get trigonometric values first for reuse
        cosine, sine, tangent = self.calculate_trig_values()
        
        # Update the data for all dynamic parts
        self.vector.set_data([0, self.real], [0, self.imag])
        self.point.set_data([self.real], [self.imag])
        
        # Update unit circle point and projections
        self.unit_point.set_data([cosine], [sine])
        self.unit_vector.set_data([0, cosine], [0, sine])
        self.cos_projection.set_data([0, cosine], [0, 0])
        self.sin_projection.set_data([0, 0], [0, sine])
        self.unit_connection.set_data([self.real, cosine], [self.imag, sine])
        
        # Update colored dotted lines for sine and cosine
        self.cos_line.set_data([cosine, cosine], [0, sine])
        self.sin_line.set_data([0, cosine], [sine, sine])
        
        # Update coordinate projections
        self.h_projection.set_data([0, self.real], [self.imag, self.imag])
        self.v_projection.set_data([self.real, self.real], [0, self.imag])
        
        # Update projection labels
        self.coord_label.set_position((self.real/2, self.imag + 0.1))
        self.coord_label.set_text(f"({self.real:.2f}, {self.imag:.2f})")
        
        self.cos_label.set_position((cosine/2, -0.1))
        self.cos_label.set_text(f"cos(θ)={cosine:.4f}")
        
        self.sin_label.set_position((-0.2, sine/2))
        self.sin_label.set_text(f"sin(θ)={sine:.4f}")
        
        # Update angle arc
        self.draw_angle_arc()
        
        # Get exact angle representations
        radian_repr, complex_repr = self.get_angle_representation(self.angle_deg)
        
        # Update complex number information text - focusing on exact representations
        info_text = f"Complex number: {self.real:.4f} + {self.imag:.4f}i\n"
        info_text += f"Magnitude: {self.mag:.4f}\n"
        
        # Display the exact angle representation when available
        info_text += f"Angle: {self.angle_deg:.4f}° = {radian_repr}\n"
        
        # Add polar form representation using cos + i*sin format
        cos_exact = self.get_exact_form(cosine)
        sin_exact = self.get_exact_form(sine)
        
        # Use exact forms if available, otherwise use decimal values
        if cos_exact != f"{cosine:.4f}" or sin_exact != f"{sine:.4f}":
            info_text += f"Polar form: {self.mag:.4f}(cos({self.angle_deg:.4f}°) + i·sin({self.angle_deg:.4f}°)) = {self.mag:.4f}({cos_exact} + i·{sin_exact})\n"
        else:
            info_text += f"Polar form: {self.mag:.4f}(cos({self.angle_deg:.4f}°) + i·sin({self.angle_deg:.4f}°))\n"
        
        # Enhanced information - show e^(iθ) form with exact representations for all 15-degree increments
        if abs(self.mag - 1.0) < 0.01:  # If close to unit circle
            # Extract the exact representation from complex_repr
            if " = " in complex_repr:
                exact_expression = complex_repr.split(" = ")[1]
                info_text += f"e^({radian_repr}i) = {exact_expression}"
            else:
                info_text += f"e^({radian_repr}i)"
        
        # Update the info text
        self.info_text.set_text(info_text)
        
        # Update trigonometric values text
        # Format tangent display
        if isinstance(tangent, float) and not np.isinf(tangent):
            tan_text = f"{tangent:.4f}"
        else:
            tan_text = "undefined" if tangent > 0 else "-undefined"
        
        # Create trigonometric values display with exact forms when available
        cos_exact = self.get_exact_form(cosine)
        sin_exact = self.get_exact_form(sine)
        
        trig_text = f"Unit Circle Values:   "
        
        # Show decimal and exact forms when different
        if cos_exact != f"{cosine:.4f}":
            trig_text += f"cos(θ) = {cosine:.4f} = {cos_exact}     "
        else:
            trig_text += f"cos(θ) = {cosine:.4f}     "
            
        if sin_exact != f"{sine:.4f}":
            trig_text += f"sin(θ) = {sine:.4f} = {sin_exact}     "
        else:
            trig_text += f"sin(θ) = {sine:.4f}     "
            
        trig_text += f"tan(θ) = {tan_text}"
        
        # Update the trigonometric values display
        self.trig_text.set_text(trig_text)
        
        # Use blitting for faster updates if possible
        if self.background is not None and not force_full_redraw and not self.need_full_redraw:
            # Restore the background
            self.fig.canvas.restore_region(self.background)
            
            # Redraw only the dynamic artists
            for artist in self.dynamic_artists:
                if artist in self.ax.get_children():
                    self.ax.draw_artist(artist)
            
            # Add the angle_label to dynamic artists if not already there
            if hasattr(self, 'angle_label') and self.angle_label not in self.dynamic_artists:
                self.dynamic_artists.append(self.angle_label)
                self.ax.draw_artist(self.angle_label)
            
            # Blit the updated display
            self.fig.canvas.blit(self.ax.bbox)
            self.fig.canvas.flush_events()
        else:
            # Fall back to standard draw if no background or full redraw needed
            self.fig.canvas.draw_idle()
            # Reset the flag
            self.need_full_redraw = False
    
    def submit_real(self, text):
        """Handle real part changes"""
        try:
            self.real = float(text)
            self.update_polar()
            self.text_mag.set_val(f"{self.mag:.4f}")
            self.text_angle.set_val(f"{self.angle_deg:.4f}")
            self.text_radians.set_val(f"{self.angle:.4f}")
            
            # Adjust scale if needed to fit the point
            self.adjust_scale_to_fit_point()
            
            self.update_display()
        except:
            pass
    
    def submit_imag(self, text):
        """Handle imaginary part changes"""
        try:
            self.imag = float(text)
            self.update_polar()
            self.text_mag.set_val(f"{self.mag:.4f}")
            self.text_angle.set_val(f"{self.angle_deg:.4f}")
            self.text_radians.set_val(f"{self.angle:.4f}")
            
            # Adjust scale if needed to fit the point
            self.adjust_scale_to_fit_point()
            
            self.update_display()
        except:
            pass
    
    def submit_mag(self, text):
        """Handle magnitude changes"""
        try:
            self.mag = float(text)
            self.update_cartesian()
            self.text_real.set_val(f"{self.real:.4f}")
            self.text_imag.set_val(f"{self.imag:.4f}")
            
            # Adjust scale if needed to fit the point
            self.adjust_scale_to_fit_point()
            
            self.update_display()
        except:
            pass
    
    def submit_angle(self, text):
        """Handle angle changes in degrees"""
        try:
            self.angle_deg = float(text)
            self.angle = np.radians(self.angle_deg)
            self.update_cartesian()
            self.text_real.set_val(f"{self.real:.4f}")
            self.text_imag.set_val(f"{self.imag:.4f}")
            self.text_radians.set_val(f"{self.angle:.4f}")
            
            # Adjust scale if needed to fit the point
            self.adjust_scale_to_fit_point()
            
            self.update_display()
        except:
            pass
    
    def submit_radians(self, text):
        """Handle angle changes in radians"""
        try:
            self.angle = float(text)
            self.angle_deg = np.degrees(self.angle)
            self.update_cartesian()
            self.text_real.set_val(f"{self.real:.4f}")
            self.text_imag.set_val(f"{self.imag:.4f}")
            self.text_angle.set_val(f"{self.angle_deg:.4f}")
            
            # Adjust scale if needed to fit the point
            self.adjust_scale_to_fit_point()
            
            self.update_display()
        except:
            pass
    
    def on_click(self, event):
        """Handle mouse clicks"""
        if event.inaxes == self.ax:
            # Check if click is near point for dragging
            if event.xdata is not None and event.ydata is not None:
                dx = event.xdata - self.real
                dy = event.ydata - self.imag
                dist = np.sqrt(dx**2 + dy**2)
                
                if dist < 0.2:  # If near point, start dragging
                    self.drag_active = True
                    self.point.set_color('orange')
                    self.point.set_markersize(10)  # Make point bigger during drag
                    # Simulate an update to show the color change
                    self.update_display()
                else:
                    # Otherwise move point to clicked location
                    self.real = event.xdata
                    self.imag = event.ydata
                    self.update_polar()
                    
                    # Adjust scale if needed to fit the point
                    self.adjust_scale_to_fit_point()
                    
                    # Update text boxes
                    self.text_real.set_val(f"{self.real:.4f}")
                    self.text_imag.set_val(f"{self.imag:.4f}")
                    self.text_mag.set_val(f"{self.mag:.4f}")
                    self.text_angle.set_val(f"{self.angle_deg:.4f}")
                    self.text_radians.set_val(f"{self.angle:.4f}")
                    
                    self.update_display()
    
    def on_motion(self, event):
        """Handle mouse motion for dragging with throttled updates"""
        if self.drag_active and event.inaxes == self.ax:
            if event.xdata is not None and event.ydata is not None:
                # Update the complex number to the dragged point
                self.real = event.xdata
                self.imag = event.ydata
                self.update_polar()
                
                # Only update the display at throttled rate during drag
                self.update_display()
    
    def on_release(self, event):
        """Handle mouse release after dragging"""
        if self.drag_active:
            self.drag_active = False
            self.point.set_color('red')
            self.point.set_markersize(8)  # Reset point size
            
            # Update text boxes
            self.text_real.set_val(f"{self.real:.4f}")
            self.text_imag.set_val(f"{self.imag:.4f}")
            self.text_mag.set_val(f"{self.mag:.4f}")
            self.text_angle.set_val(f"{self.angle_deg:.4f}")
            self.text_radians.set_val(f"{self.angle:.4f}")
            
            # Full update
            self.update_display(force_full_redraw=True)
    
    def zoom_in(self, event):
        """Zoom in by 25% while keeping origin centered"""
        # Reduce scale by 25%
        self.view_scale *= 0.75
        
        # Update the view limits with origin centered
        self.ax.set_xlim(-self.view_scale, self.view_scale)
        self.ax.set_ylim(-self.view_scale, self.view_scale)
        
        # Need a full redraw for this operation
        self.need_full_redraw = True
        self.update_display(force_full_redraw=True)
    
    def zoom_out(self, event):
        """Zoom out by 25% while keeping origin centered"""
        # Increase scale by 25%
        self.view_scale *= 1.25
        
        # Update the view limits with origin centered
        self.ax.set_xlim(-self.view_scale, self.view_scale)
        self.ax.set_ylim(-self.view_scale, self.view_scale)
        
        # Need a full redraw for this operation
        self.need_full_redraw = True
        self.update_display(force_full_redraw=True)
    
    def adjust_scale_to_fit_point(self):
        """Adjust the scale if the current point is outside the visible area"""
        # Calculate the minimum scale needed to show the point
        point_max = max(abs(self.real), abs(self.imag))
        
        # Add 20% margin
        needed_scale = point_max * 1.2
        
        # Only increase scale if needed
        if needed_scale > self.view_scale:
            self.view_scale = needed_scale
            self.ax.set_xlim(-self.view_scale, self.view_scale)
            self.ax.set_ylim(-self.view_scale, self.view_scale)
            # Need a full redraw when changing axis limits
            self.need_full_redraw = True
    
    def calculate_trig_values(self):
        """Calculate trigonometric values from the current angle"""
        # For accuracy, calculate trig values directly from the angle
        sine = np.sin(self.angle)
        cosine = np.cos(self.angle)
        # Handle undefined tangent at 90° and 270°
        if abs(cosine) < 1e-10:
            tangent = float('inf') if sine > 0 else float('-inf')
        else:
            tangent = sine / cosine
        return cosine, sine, tangent
    
    def set_angle_15(self, event):
        """Set angle to 15 degrees"""
        self.angle = np.radians(15)
        self.angle_deg = 15
        self.update_cartesian()
        self.text_real.set_val(f"{self.real:.4f}")
        self.text_imag.set_val(f"{self.imag:.4f}")
        self.text_mag.set_val(f"{self.mag:.4f}")
        self.text_angle.set_val(f"{self.angle_deg:.4f}")
        self.text_radians.set_val(f"{self.angle:.4f}")
        self.update_display()
    
    def set_angle_30(self, event):
        """Set angle to 30 degrees"""
        self.angle = np.radians(30)
        self.angle_deg = 30
        self.update_cartesian()
        self.text_real.set_val(f"{self.real:.4f}")
        self.text_imag.set_val(f"{self.imag:.4f}")
        self.text_mag.set_val(f"{self.mag:.4f}")
        self.text_angle.set_val(f"{self.angle_deg:.4f}")
        self.text_radians.set_val(f"{self.angle:.4f}")
        self.update_display()
    
    def set_angle_45(self, event):
        """Set angle to 45 degrees"""
        self.angle = np.radians(45)
        self.angle_deg = 45
        self.update_cartesian()
        self.text_real.set_val(f"{self.real:.4f}")
        self.text_imag.set_val(f"{self.imag:.4f}")
        self.text_mag.set_val(f"{self.mag:.4f}")
        self.text_angle.set_val(f"{self.angle_deg:.4f}")
        self.text_radians.set_val(f"{self.angle:.4f}")
        self.update_display()
    
    def set_angle_60(self, event):
        """Set angle to 60 degrees"""
        self.angle = np.radians(60)
        self.angle_deg = 60
        self.update_cartesian()
        self.text_real.set_val(f"{self.real:.4f}")
        self.text_imag.set_val(f"{self.imag:.4f}")
        self.text_mag.set_val(f"{self.mag:.4f}")
        self.text_angle.set_val(f"{self.angle_deg:.4f}")
        self.text_radians.set_val(f"{self.angle:.4f}")
        self.update_display()
    
    def set_angle_90(self, event):
        """Set angle to 90 degrees"""
        self.angle = np.radians(90)
        self.angle_deg = 90
        self.update_cartesian()
        self.text_real.set_val(f"{self.real:.4f}")
        self.text_imag.set_val(f"{self.imag:.4f}")
        self.text_mag.set_val(f"{self.mag:.4f}")
        self.text_angle.set_val(f"{self.angle_deg:.4f}")
        self.text_radians.set_val(f"{self.angle:.4f}")
        self.update_display()
    
    def set_angle_120(self, event):
        """Set angle to 120 degrees"""
        self.angle = np.radians(120)
        self.angle_deg = 120
        self.update_cartesian()
        self.text_real.set_val(f"{self.real:.4f}")
        self.text_imag.set_val(f"{self.imag:.4f}")
        self.text_mag.set_val(f"{self.mag:.4f}")
        self.text_angle.set_val(f"{self.angle_deg:.4f}")
        self.text_radians.set_val(f"{self.angle:.4f}")
        self.update_display()
    
    def set_angle_180(self, event):
        """Set angle to 180 degrees"""
        self.angle = np.radians(180)
        self.angle_deg = 180
        self.update_cartesian()
        self.text_real.set_val(f"{self.real:.4f}")
        self.text_imag.set_val(f"{self.imag:.4f}")
        self.text_mag.set_val(f"{self.mag:.4f}")
        self.text_angle.set_val(f"{self.angle_deg:.4f}")
        self.text_radians.set_val(f"{self.angle:.4f}")
        self.update_display()
    
    def set_angle_270(self, event):
        """Set angle to 270 degrees"""
        self.angle = np.radians(270)
        self.angle_deg = 270
        self.update_cartesian()
        self.text_real.set_val(f"{self.real:.4f}")
        self.text_imag.set_val(f"{self.imag:.4f}")
        self.text_mag.set_val(f"{self.mag:.4f}")
        self.text_angle.set_val(f"{self.angle_deg:.4f}")
        self.text_radians.set_val(f"{self.angle:.4f}")
        self.update_display()
    
    def set_angle_360(self, event):
        """Set angle to 360 degrees"""
        self.angle = np.radians(360)
        self.angle_deg = 360
        self.update_cartesian()
        self.text_real.set_val(f"{self.real:.4f}")
        self.text_imag.set_val(f"{self.imag:.4f}")
        self.text_mag.set_val(f"{self.mag:.4f}")
        self.text_angle.set_val(f"{self.angle_deg:.4f}")
        self.text_radians.set_val(f"{self.angle:.4f}")
        self.update_display()
    
    def toggle_dark_mode(self, event):
        """Toggle between dark mode and light mode"""
        self.dark_mode = not self.dark_mode
        
        if self.dark_mode:
            # Dark mode settings
            self.fig.set_facecolor('#2c2c2c')
            self.ax.set_facecolor('#2c2c2c')
            self.ax.grid(True, color='#5b5b5b')
            self.ax.set_title('Complex Number Visualization', color='white')
            self.ax.set_xlabel('Real (Cosine Component)', color='white')
            self.ax.set_ylabel('Imaginary (Sine Component)', color='white')
            self.ax.tick_params(axis='x', colors='white')
            self.ax.tick_params(axis='y', colors='white')
            self.ax.spines['bottom'].set_color('white')
            self.ax.spines['top'].set_color('white')
            self.ax.spines['left'].set_color('white')
            self.ax.spines['right'].set_color('white')
            
            # Update text elements
            self.info_text.set_bbox(dict(facecolor='#444444', alpha=0.8))
            self.info_text.set_color('white')
            self.trig_text.set_bbox(dict(facecolor='#444444', alpha=0.9))
            self.trig_text.set_color('white')
            self.coord_label.set_color('skyblue')
            self.cos_label.set_color('lightgreen')
            self.sin_label.set_color('lightgreen')
            
            # Update angle label
            if hasattr(self, 'angle_label'):
                self.angle_label.set_bbox(dict(facecolor='#444444', alpha=0.7))
                self.angle_label.set_color('white')
            
            # Update special angle labels
            for label in self.special_angle_labels:
                label.set_bbox(dict(facecolor='#444444', alpha=0.7))
                label.set_color('white')
            
            # Update button text
            self.btn_dark_mode.label.set_text('Light Mode')
        else:
            # Light mode settings
            self.fig.set_facecolor('white')
            self.ax.set_facecolor('white')
            self.ax.grid(True, color='lightgray')
            self.ax.set_title('Complex Number Visualization', color='black')
            self.ax.set_xlabel('Real (Cosine Component)', color='black')
            self.ax.set_ylabel('Imaginary (Sine Component)', color='black')
            self.ax.tick_params(axis='x', colors='black')
            self.ax.tick_params(axis='y', colors='black')
            self.ax.spines['bottom'].set_color('black')
            self.ax.spines['top'].set_color('black')
            self.ax.spines['left'].set_color('black')
            self.ax.spines['right'].set_color('black')
            
            # Update text elements
            self.info_text.set_bbox(dict(facecolor='white', alpha=0.8))
            self.info_text.set_color('black')
            self.trig_text.set_bbox(dict(facecolor='white', alpha=0.9))
            self.trig_text.set_color('black')
            self.coord_label.set_color('blue')
            self.cos_label.set_color('green')
            self.sin_label.set_color('green')
            
            # Update angle label
            if hasattr(self, 'angle_label'):
                self.angle_label.set_bbox(dict(facecolor='white', alpha=0.7))
                self.angle_label.set_color('black')
            
            # Update special angle labels
            for label in self.special_angle_labels:
                label.set_bbox(dict(facecolor='white', alpha=0.7))
                label.set_color('black')
            
            # Update button text
            self.btn_dark_mode.label.set_text('Dark Mode')
        
        # Need a full redraw after changing theme
        self.need_full_redraw = True
        self.update_display(force_full_redraw=True)

    def get_angle_representation(self, angle_deg):
        """
        Get string representation of angle in terms of π if it's a special value
        Returns tuple of (radian_repr, complex_repr)
        """
        angle_rad = np.radians(angle_deg)
        
        # Define exact representations for angles in 15-degree increments
        # Format: (numerator, denominator, radian representation, e^(ix) representation with exact values)
        special_angles = {
            0: (0, 1, "0", "e^(0i) = 1"),
            15: (1, 12, "π/12", "e^(iπ/12) = cos(π/12) + i·sin(π/12) = (√6+√2)/4 + i·(√6-√2)/4"),
            30: (1, 6, "π/6", "e^(iπ/6) = cos(π/6) + i·sin(π/6) = √3/2 + i·1/2"),
            45: (1, 4, "π/4", "e^(iπ/4) = cos(π/4) + i·sin(π/4) = 1/√2 + i·1/√2"),
            60: (1, 3, "π/3", "e^(iπ/3) = cos(π/3) + i·sin(π/3) = 1/2 + i·√3/2"),
            75: (5, 12, "5π/12", "e^(i5π/12) = cos(5π/12) + i·sin(5π/12) = (√6-√2)/4 + i·(√6+√2)/4"),
            90: (1, 2, "π/2", "e^(iπ/2) = cos(π/2) + i·sin(π/2) = 0 + i·1 = i"),
            105: (7, 12, "7π/12", "e^(i7π/12) = cos(7π/12) + i·sin(7π/12) = -(√6-√2)/4 + i·(√6+√2)/4"),
            120: (2, 3, "2π/3", "e^(i2π/3) = cos(2π/3) + i·sin(2π/3) = -1/2 + i·√3/2"),
            135: (3, 4, "3π/4", "e^(i3π/4) = cos(3π/4) + i·sin(3π/4) = -1/√2 + i·1/√2"),
            150: (5, 6, "5π/6", "e^(i5π/6) = cos(5π/6) + i·sin(5π/6) = -√3/2 + i·1/2"),
            165: (11, 12, "11π/12", "e^(i11π/12) = cos(11π/12) + i·sin(11π/12) = -(√6+√2)/4 + i·(√6-√2)/4"),
            180: (1, 1, "π", "e^(iπ) = cos(π) + i·sin(π) = -1 + i·0 = -1"),
            195: (13, 12, "13π/12", "e^(i13π/12) = cos(13π/12) + i·sin(13π/12) = -(√6+√2)/4 - i·(√6-√2)/4"),
            210: (7, 6, "7π/6", "e^(i7π/6) = cos(7π/6) + i·sin(7π/6) = -√3/2 - i·1/2"),
            225: (5, 4, "5π/4", "e^(i5π/4) = cos(5π/4) + i·sin(5π/4) = -1/√2 - i·1/√2"),
            240: (4, 3, "4π/3", "e^(i4π/3) = cos(4π/3) + i·sin(4π/3) = -1/2 - i·√3/2"),
            255: (17, 12, "17π/12", "e^(i17π/12) = cos(17π/12) + i·sin(17π/12) = -(√6-√2)/4 - i·(√6+√2)/4"),
            270: (3, 2, "3π/2", "e^(i3π/2) = cos(3π/2) + i·sin(3π/2) = 0 - i·1 = -i"),
            285: (19, 12, "19π/12", "e^(i19π/12) = cos(19π/12) + i·sin(19π/12) = (√6-√2)/4 - i·(√6+√2)/4"),
            300: (5, 3, "5π/3", "e^(i5π/3) = cos(5π/3) + i·sin(5π/3) = 1/2 - i·√3/2"),
            315: (7, 4, "7π/4", "e^(i7π/4) = cos(7π/4) + i·sin(7π/4) = 1/√2 - i·1/√2"),
            330: (11, 6, "11π/6", "e^(i11π/6) = cos(11π/6) + i·sin(11π/6) = √3/2 - i·1/2"),
            345: (23, 12, "23π/12", "e^(i23π/12) = cos(23π/12) + i·sin(23π/12) = (√6+√2)/4 - i·(√6-√2)/4"),
            360: (2, 1, "2π", "e^(i2π) = cos(2π) + i·sin(2π) = 1 + i·0 = 1")
        }
        
        # Normalize angle to 0-360 range
        norm_angle = angle_deg % 360
        
        # Check if the angle is in our special angles dictionary or very close to it
        for special_angle, (num, denom, rad_repr, complex_repr) in special_angles.items():
            if abs(norm_angle - special_angle) < 0.01:
                return rad_repr, complex_repr
        
        # If not a special angle, just return the normalized angle in degrees
        return f"{norm_angle:.4f}", f"e^({norm_angle:.2f}°i)"

    def on_hover(self, event):
        """Handle mouse hover over special angle points"""
        # Skip if not in the axis or if no special angle data
        if not hasattr(self, 'special_angle_data') or event.inaxes != self.ax:
            # Hide all special angle labels when mouse is not over any point
            for label in self.special_angle_labels:
                label.set_visible(False)
            self.fig.canvas.draw_idle()
            return
        
        # Check if mouse is over any special angle point
        hover_found = False
        for (x, y), label in self.special_angle_data.items():
            # Calculate distance from mouse to point
            dist = np.sqrt((event.xdata - x)**2 + (event.ydata - y)**2)
            
            # Show label if mouse is close to point (0.1 units)
            if dist < 0.1:
                label.set_visible(True)
                hover_found = True
            else:
                label.set_visible(False)
        
        # Update the display if any hover state changed
        if hover_found or any(label.get_visible() for label in self.special_angle_labels):
            self.fig.canvas.draw_idle()

# Main program
if __name__ == "__main__":
    plt.ion()  # Turn on interactive mode
    plane = SimpleComplexPlane()
    plt.show(block=True)  # This keeps the window open but non-blocking
