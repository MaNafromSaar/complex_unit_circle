#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Complex Unit Circle Visualization - Executable Setup File

This script creates a standalone executable for the Complex Unit Circle Visualization application.

Author: Matthias Naumann
Version: 1.0.0
Date: May 11, 2025
License: MIT License
"""

import os
import sys

def create_executable():
    """Create a standalone executable using PyInstaller"""
    try:
        # Import PyInstaller
        import PyInstaller.__main__

        # Get the directory of the current script
        current_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Path to the main script
        main_script = os.path.join(current_dir, "complex_unit_circle.py")
        
        # Create the executable
        PyInstaller.__main__.run([
            main_script,
            '--name=ComplexUnitCircle',
            '--windowed',  # Don't show console window
            '--onefile',   # Create a single executable file
            '--add-data=LICENSE;.',  # Include LICENSE file
            '--add-data=README.md;.',  # Include README file
            '--icon=NONE',  # No icon for now, replace with icon path if available
            '--clean',     # Clean PyInstaller cache
            '--noconfirm'  # Replace existing output without asking
        ])
        
        print("\nSuccess! Executable created in the 'dist' folder.")
        print("Look for ComplexUnitCircle.exe in the dist folder.")
        
    except ImportError:
        print("PyInstaller is not installed. Installing it now...")
        try:
            import pip
            pip.main(['install', 'pyinstaller'])
            print("PyInstaller installed successfully. Run this script again to create the executable.")
        except Exception as e:
            print(f"Error installing PyInstaller: {e}")
            print("Please install PyInstaller manually by running: pip install pyinstaller")
    except Exception as e:
        print(f"Error creating executable: {e}")

if __name__ == "__main__":
    create_executable()
