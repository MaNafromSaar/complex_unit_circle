#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Complex Unit Circle Visualization

An interactive visualization of complex numbers and their relationship to the unit circle.
Shows the connection between complex numbers, trigonometric functions, and polar form.

Author: Matthias Naumann
Version: 1.0.0
Date: May 6, 2025
License: MIT License (see LICENSE file or https://opensource.org/licenses/MIT)
Repository: [Your GitHub Repository URL]

This software is free to use for educational purposes.
"""

# An implementation of complex unit circle visualization
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
        
        # Add trigonometric values display directly under the graph
        self.trig_text = self.ax.text(
            0.5, -0.1, '',  
            transform=self.ax.transAxes,
            ha='center', va='top',
            bbox=dict(facecolor='white', alpha=0.9),
            fontsize=9
        )
        
        # Add text displays below the graph
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
        
        # Add angle arc (after dynamic_artists is created)
        self.draw_angle_arc()

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
            if angle_deg > 180 and angle_deg <= 360:
                negative_angle = angle_deg - 360
                label_text = f"{angle_deg}° = {negative_angle}°\n"
                
                # Calculate negative radian representation
                neg_radian_repr = radian_repr
                if radian_repr == "π":
                    neg_radian_repr = "-π"
                elif radian_repr == "2π":
                    neg_radian_repr = "0"
                elif "/" in radian_repr:
                    try:
                        # Extract numerator and handle negative representation
                        parts = radian_repr.split("π")
                        if "/" in parts[0]:
                            num_str, denom_str = parts[0].split("/")
                            num = int(num_str)
                            denom = int(denom_str)
                            # Calculate negative form
                            neg_num = num - (2 * denom)
                            if neg_num == -denom:
                                neg_radian_repr = "-π"
                            else:
                                neg_radian_repr = f"{neg_num}/{denom}π"
                        elif parts[0].strip() == "":
                            neg_radian_repr = "-π"
                        else:
                            num = int(parts[0])
                            neg_radian_repr = f"{num-2}π"
                    except:
                        # If parsing fails, just use the original
                        pass
                
                label_text += f"{radian_repr} = {neg_radian_repr}\n"
            else:
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
        
        # Create arc points
        theta = np.linspace(0, self.angle, 50)
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        
        # Update or create arc line
        if hasattr(self, 'angle_arc'):
            self.angle_arc.set_data(x, y)
        else:
            self.angle_arc, = self.ax.plot(x, y, 'b-', lw=1.5)
            self.dynamic_artists.append(self.angle_arc)  # Add to dynamic artists list for blitting
            
        # Calculate midpoint for angle label positioning
        mid_angle = self.angle / 2
        label_x = (radius * 1.2) * np.cos(mid_angle)
        label_y = (radius * 1.2) * np.sin(mid_angle)
        
        # Format angle text
        angle_text = f"{self.angle_deg:.2f}°"
        
        # Add negative angle representation for angles between 180° and 360°
        if 180 < self.angle_deg <= 360:
            negative_angle = self.angle_deg - 360
            angle_text = f"{self.angle_deg:.2f}° = {negative_angle:.2f}°"
        
        # Get angle representations once to avoid recalculation
        radian_repr, complex_repr = self.get_angle_representation(self.angle_deg)
        
        # Add radian representation to angle text for special angles
        if "π" in radian_repr:
            # For angles between 180° and 360°, show both positive and negative radian forms
            if 180 < self.angle_deg <= 360:
                # Calculate negative radian representation
                neg_radian_repr = radian_repr
                if radian_repr == "π":
                    neg_radian_repr = "-π"
                elif radian_repr == "2π":
                    neg_radian_repr = "0"
                elif "/" in radian_repr:
                    try:
                        # Extract numerator and handle negative representation
                        parts = radian_repr.split("π")
                        if "/" in parts[0]:
                            num_str, denom_str = parts[0].split("/")
                            num = int(num_str)
                            denom = int(denom_str)
                            # Calculate negative form
                            neg_num = num - (2 * denom)
                            if neg_num == -denom:
                                neg_radian_repr = "-π"
                            else:
                                neg_radian_repr = f"{neg_num}/{denom}π"
                        elif parts[0].strip() == "":
                            neg_radian_repr = "-π"
                        else:
                            num = int(parts[0])
                            neg_radian_repr = f"{num-2}π"
                    except:
                        # If parsing fails, just use the original
                        pass
                angle_text += f"\n{radian_repr} = {neg_radian_repr}"
            else:
                angle_text += f"\n{radian_repr}"
        if 180 < self.angle_deg <= 360:
            negative_angle = self.angle_deg - 360
            angle_text = f"{self.angle_deg:.2f}° or {negative_angle:.2f}°"
            if "π" in radian_repr:
                angle_text += f"\n{radian_repr}"
        
        # Add complex representation for unit circle angles
        if abs(self.mag - 1.0) < 0.01:  # If close to unit circle
            complex_parts = complex_repr.split(' = ')
            if len(complex_parts) > 1:
                angle_text += f"\n{complex_parts[0]}"
        
        # Update or create angle label
        if hasattr(self, 'angle_label'):
            # Update existing label instead of removing and recreating
            self.angle_label.set_position((label_x, label_y))
            self.angle_label.set_text(angle_text)
            
            # Update color and background based on dark mode
            label_bg_color = '#444444' if self.dark_mode else 'white'
            label_text_color = 'white' if self.dark_mode else 'black'
            
            self.angle_label.set_bbox(dict(facecolor=label_bg_color, alpha=0.8, boxstyle='round,pad=0.3'))
            self.angle_label.set_color(label_text_color)
        else:
            # Initialize label for the first time
            label_bg_color = '#444444' if self.dark_mode else 'white'
            label_text_color = 'white' if self.dark_mode else 'black'
            
            self.angle_label = self.ax.text(
                label_x, label_y, angle_text, 
                fontsize=9, ha='center', va='center',
                bbox=dict(facecolor=label_bg_color, alpha=0.8, boxstyle='round,pad=0.3'),
                color=label_text_color
            )
            self.dynamic_artists.append(self.angle_label)  # Add to dynamic artists list for blitting

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
        
        # Prevent excessive updates from causing recursion
        if hasattr(self, '_update_in_progress') and self._update_in_progress:
            return
        self._update_in_progress = True
        
        try:
            # Cache background if we don't have it yet (for blitting optimization)
            if self.background is None or force_full_redraw:
                self.fig.canvas.draw()
                self.background = self.fig.canvas.copy_from_bbox(self.ax.bbox)
                force_full_redraw = True  # Need a full redraw on first pass
            
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
            
            # Update angle arc - using optimized calculation
            self.draw_angle_arc()
            
            # Get exact angle representations - using cached values for performance
            radian_repr, complex_repr = self.get_angle_representation(self.angle_deg)
            
            # Update complex number information text - focusing on exact representations
            info_text = f"Complex number: {self.real:.4f} + {self.imag:.4f}i\n"
            info_text += f"Magnitude: {self.mag:.4f}\n"
            
            # Display angle in degrees with appropriate representation
            if 180 < self.angle_deg <= 360:
                negative_angle = self.angle_deg - 360
                info_text += f"Angle: {self.angle_deg:.4f}° = {negative_angle:.4f}°\n"
            else:
                info_text += f"Angle: {self.angle_deg:.4f}°\n"
            
            # Add radian representation on its own line with both positive and negative representation when applicable
            if 180 < self.angle_deg <= 360:
                negative_radians = self.angle - 2*np.pi
                info_text += f"Angle in radians: {self.angle:.4f} = {negative_radians:.4f}"
                if "π" in radian_repr:
                    # Get negative radian representation
                    neg_radian_repr = radian_repr
                    if "π" in radian_repr and "/" in radian_repr:
                        # For fractions, calculate the negative form
                        if radian_repr == "π":
                            neg_radian_repr = "-π"
                        elif radian_repr == "2π":
                            neg_radian_repr = "0"
                        else:
                            try:
                                # Extract numerator and handle negative representation
                                parts = radian_repr.split("π")
                                if "/" in parts[0]:
                                    num_str, denom_str = parts[0].split("/")
                                    num = int(num_str)
                                    denom = int(denom_str)
                                    # Calculate negative form
                                    neg_num = num - (2 * denom)
                                    if neg_num == -denom:
                                        neg_radian_repr = "-π"
                                    else:
                                        neg_radian_repr = f"{neg_num}/{denom}π"
                                elif parts[0].strip() == "":
                                    neg_radian_repr = "-π"
                                else:
                                    num = int(parts[0])
                                    neg_radian_repr = f"{num-2}π"
                            except:
                                # If parsing fails, just use the decimal form
                                pass
                    info_text += f" = {radian_repr} = {neg_radian_repr}"
                info_text += "\n"
            else:
                info_text += f"Angle in radians: {self.angle:.4f}"
                if "π" in radian_repr:
                    info_text += f" = {radian_repr}"
                info_text += "\n"
            
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
                
                # Blit the updated display
                self.fig.canvas.blit(self.ax.bbox)
                self.fig.canvas.flush_events()
            else:
                # Fall back to standard draw if no background or full redraw needed
                self.fig.canvas.draw_idle()
                # Capture a new background for future blitting
                if not self.drag_active:
                    self.background = self.fig.canvas.copy_from_bbox(self.ax.bbox)
                # Reset the flag
                self.need_full_redraw = False
                
        finally:
            # Always reset the in-progress flag when done
            self._update_in_progress = False
    
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
        # Prevent recursive calls
        if hasattr(self, '_submitting_angle') and self._submitting_angle:
            return
        self._submitting_angle = True
        
        try:
            # Parse input angle in degrees - can be any value, negative or positive
            input_angle_deg = float(text)
            
            # Round for precision
            input_angle_deg = round(input_angle_deg, 6)
            
            # Normalize to 0-360 range for internal representation
            self.angle_deg = input_angle_deg % 360
            self.angle = np.radians(self.angle_deg)
            self.update_cartesian()
            
            # Batch update text boxes to avoid multiple redraws
            self.text_real.set_val(f"{self.real:.4f}")
            self.text_imag.set_val(f"{self.imag:.4f}")
            self.text_radians.set_val(f"{self.angle:.4f}")
            
            # Only update angle box if it's significantly different
            # This prevents recursive calls when setting it to a similar value
            displayed_angle = float(self.text_angle.text)
            if abs(displayed_angle - self.angle_deg) > 0.001:
                self.text_angle.set_val(f"{self.angle_deg:.4f}")
            
            # Adjust scale if needed to fit the point
            self.adjust_scale_to_fit_point()
            
            # Update with optimizations for special angles
            self.update_display(force_full_redraw=False)
        except ValueError:
            # Handle invalid input gracefully
            pass
        finally:
            self._submitting_angle = False
    
    def submit_radians(self, text):
        """Handle angle changes in radians"""
        # Prevent recursive calls
        if hasattr(self, '_submitting_radians') and self._submitting_radians:
            return
        self._submitting_radians = True
        
        try:
            self.angle = float(text)
            self.angle_deg = np.degrees(self.angle)
            
            # Round for precision
            self.angle_deg = round(self.angle_deg, 6)
            self.angle = round(self.angle, 6)
            
            self.update_cartesian()
            
            # Batch update text boxes to avoid multiple redraws
            self.text_real.set_val(f"{self.real:.4f}")
            self.text_imag.set_val(f"{self.imag:.4f}")
            
            # Only update angle box if it's significantly different
            displayed_angle = float(self.text_angle.text)
            if abs(displayed_angle - self.angle_deg) > 0.001:
                self.text_angle.set_val(f"{self.angle_deg:.4f}")
            
            # Adjust scale if needed to fit the point
            self.adjust_scale_to_fit_point()
            
            # Update with optimizations for special angles
            self.update_display(force_full_redraw=False)
        except ValueError:
            # Handle invalid input gracefully
            pass
        finally:
            self._submitting_radians = False
    
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
        
        # Only increase scale if needed and if the change is significant
        if needed_scale > self.view_scale and (needed_scale - self.view_scale) > 0.01:
            self.view_scale = needed_scale
            self.ax.set_xlim(-self.view_scale, self.view_scale)
            self.ax.set_ylim(-self.view_scale, self.view_scale)
            # Need a full redraw when changing axis limits
            self.need_full_redraw = True
            # Clear the background to force redrawing
            self.background = None

    def calculate_trig_values(self):
        """Calculate trigonometric values from the current angle"""
        # For critical angles that might cause numerical precision issues, 
        # use exact values instead of floating-point calculations
        
        # Normalize angle to 0-360 range and round to handle floating point issues
        norm_angle = round(self.angle_deg % 360, 6)
        
        # Special handling for common angles that might have precision issues
        if abs(norm_angle - 0) < 0.000001 or abs(norm_angle - 360) < 0.000001:
            return 1.0, 0.0, 0.0
        elif abs(norm_angle - 90) < 0.000001:
            return 0.0, 1.0, float('inf')
        elif abs(norm_angle - 180) < 0.000001:
            return -1.0, 0.0, 0.0
        elif abs(norm_angle - 270) < 0.000001:
            return 0.0, -1.0, float('inf')
            
        # For other angles, calculate using numpy for accuracy
        sine = np.sin(self.angle)
        cosine = np.cos(self.angle)
        
        # Handle very small values that should be zero (floating point cleanup)
        if abs(sine) < 1e-13:
            sine = 0.0
        if abs(cosine) < 1e-13:
            cosine = 0.0
            
        # Handle undefined tangent carefully
        if abs(cosine) < 1e-10:
            tangent = float('inf') if sine > 0 else float('-inf')
        else:
            tangent = sine / cosine
            
        return cosine, sine, tangent
    
    def set_angle(self, angle_deg):
        """Set angle to a specific value in degrees - common method for all angle presets"""
        # Prevent recursive calls
        if hasattr(self, '_setting_angle') and self._setting_angle:
            return
        self._setting_angle = True
        
        try:
            # Round the angle to help with floating-point precision issues
            angle_deg = round(angle_deg, 6)
            
            # Set the new angle
            self.angle_deg = angle_deg % 360  # Normalize to 0-360
            self.angle = np.radians(self.angle_deg)
            
            # Update coordinates with mag=1 (unit circle)
            self.mag = 1.0
            self.update_cartesian()
            
            # Batch update text boxes to avoid multiple redraws
            self.text_real.set_val(f"{self.real:.4f}")
            self.text_imag.set_val(f"{self.imag:.4f}")
            self.text_mag.set_val(f"{self.mag:.4f}")
            self.text_angle.set_val(f"{self.angle_deg:.4f}")
            self.text_radians.set_val(f"{self.angle:.4f}")
            
            # Update display with optimizations for special angles
            self.update_display(force_full_redraw=False)
        finally:
            self._setting_angle = False
    
    def set_angle_15(self, event):
        """Set angle to 15 degrees"""
        self.set_angle(15)
    
    def set_angle_30(self, event):
        """Set angle to 30 degrees"""
        self.set_angle(30)
    
    def set_angle_45(self, event):
        """Set angle to 45 degrees"""
        self.set_angle(45)
    
    def set_angle_60(self, event):
        """Set angle to 60 degrees"""
        self.set_angle(60)
    
    def set_angle_90(self, event):
        """Set angle to 90 degrees"""
        self.set_angle(90)
    
    def set_angle_120(self, event):
        """Set angle to 120 degrees"""
        self.set_angle(120)
    
    def set_angle_180(self, event):
        """Set angle to 180 degrees"""
        self.set_angle(180)
    
    def set_angle_270(self, event):
        """Set angle to 270 degrees"""
        self.set_angle(270)
    
    def set_angle_360(self, event):
        """Set angle to 360 degrees"""
        self.set_angle(360)
    
    def toggle_dark_mode(self, event):
        """Toggle between dark mode and light mode"""
        # Prevent recursive updates
        if hasattr(self, '_toggling_mode') and self._toggling_mode:
            return
        self._toggling_mode = True
        
        try:
            self.dark_mode = not self.dark_mode
            
            # Set colors based on dark/light mode
            bg_color = '#2c2c2c' if self.dark_mode else 'white'
            text_color = 'white' if self.dark_mode else 'black'
            grid_color = '#5b5b5b' if self.dark_mode else 'lightgray'
            box_color = '#444444' if self.dark_mode else 'white'
            
            # Update figure and axis colors
            self.fig.set_facecolor(bg_color)
            self.ax.set_facecolor(bg_color)
            self.ax.grid(True, color=grid_color)
            self.ax.set_title('Complex Number Visualization', color=text_color)
            self.ax.set_xlabel('Real (Cosine Component)', color=text_color)
            self.ax.set_ylabel('Imaginary (Sine Component)', color=text_color)
            
            # Update tick parameters
            self.ax.tick_params(axis='x', colors=text_color)
            self.ax.tick_params(axis='y', colors=text_color)
            
            # Update spines
            for spine in self.ax.spines.values():
                spine.set_color(text_color)
            
            # Update text elements
            self.info_text.set_bbox(dict(facecolor=box_color, alpha=0.8))
            self.info_text.set_color(text_color)
            self.trig_text.set_bbox(dict(facecolor=box_color, alpha=0.9))
            self.trig_text.set_color(text_color)
            
            # Update label colors
            self.coord_label.set_color('skyblue' if self.dark_mode else 'blue')
            self.cos_label.set_color('lightgreen' if self.dark_mode else 'green')
            self.sin_label.set_color('lightgreen' if self.dark_mode else 'green')
            
            # Update angle label if it exists
            if hasattr(self, 'angle_label'):
                self.angle_label.set_bbox(dict(facecolor=box_color, alpha=0.8, boxstyle='round,pad=0.3'))
                self.angle_label.set_color(text_color)
            
            # Update special angle labels
            if hasattr(self, 'special_angle_labels'):
                for label in self.special_angle_labels:
                    label.set_bbox(dict(facecolor=box_color, alpha=0.8, boxstyle='round,pad=0.3'))
                    label.set_color(text_color)
            
            # Update button text
            self.btn_dark_mode.label.set_text('Light Mode' if self.dark_mode else 'Dark Mode')
            
            # Force a full redraw
            self.background = None  # Force regeneration of background
            self.need_full_redraw = True
            self.update_display(force_full_redraw=True)
        finally:
            self._toggling_mode = False

    def get_angle_representation(self, angle_deg):
        """
        Get string representation of angle in terms of π if it's a special value
        Returns tuple of (radian_repr, complex_repr)
        """
        # Normalize angle to 0-360 range
        norm_angle = angle_deg % 360
        
        # Round to 6 decimals to avoid floating point comparison issues
        norm_angle = round(norm_angle, 6)
        
        # Create the special angles dictionary only once
        if not hasattr(self, 'special_angles_cache'):
            # Define exact representations for angles in 15-degree increments
            # Format: (numerator, denominator, radian representation, e^(ix) representation with exact values)
            self.special_angles_cache = {
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
            
            # Add second-level cache for special presets to ensure immediate matches for 
            # values that might have tiny numerical discrepancies
            self.special_angles_rounded_cache = {}
            for angle, data in self.special_angles_cache.items():
                self.special_angles_rounded_cache[round(angle, 6)] = data
        
        # Fast path for exact angles - common case for preset buttons
        if norm_angle in self.special_angles_cache:
            _, _, rad_repr, complex_repr = self.special_angles_cache[norm_angle]
            return rad_repr, complex_repr
        
        # Second fast path for rounded values - handles floating point precision issues
        if norm_angle in self.special_angles_rounded_cache:
            _, _, rad_repr, complex_repr = self.special_angles_rounded_cache[norm_angle]
            return rad_repr, complex_repr
            
        # For preset buttons that might have floating point issues, use an even more precise check
        # This is especially important for 270° and 45° which are causing the recursion error
        for special_angle, (_, _, rad_repr, complex_repr) in self.special_angles_cache.items():
            if abs(norm_angle - special_angle) < 0.000001:  # Much tighter tolerance
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