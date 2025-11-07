# -*- coding: utf-8 -*-
"""
Created on Fri Oct 24 09:02:15 2025

@author: sydney
import data strema from phypox
"""

import requests
import time

# Replace with your Phyphox device IP
url = "http://10.143.194.122:8080/get?accX&accY"

while True:
    # Use a short timeout and stream=True
    response = requests.get(url)
    data = response.json()
    x = data['buffer']['accX']['buffer'][0]
    y = data['buffer']['accY']['buffer'][0]
    print(x,y)
