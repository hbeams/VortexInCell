#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
class RK2:
    def __init__(self):
        self._m_RHS = []
    def __init__(self, RHS):
        self._m_RHS = RHS
    def advance(self, time, dt, currentState):
        k = self._m_RHS.compute(time, currentState)*dt
        midState = currentState + k*(1.0/2.0)
        k = self._m_RHS.compute(time + (1.0/2.0)*dt, midState)*dt
        finalState = currentState + k
        finalState.clean()
        return finalState
