#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The script generates an animation from png's (e.g. spectra plots).

20240130
@author: lothar
"""

from PIL import Image
import glob


files = input('Enter the name of the png images: ')
imgs = glob.glob(files)
imgs.sort()

frames = []

for i in imgs:
    frames.append(Image.open(i))

frames[0].save('animation.gif', format='GIF',
               append_images=frames[1:],
               save_all=True,
               duration=500, loop=0)  # duration anpassen
