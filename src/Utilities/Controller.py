# Copyright (C) 2026  Cody Raskin

import json
from Integrators import *
from Colors import BRIGHT_CYAN, BRIGHT_GREEN, BRIGHT_YELLOW, DIM, colorize

class Controller:
    def __init__(self,integrator,periodicWork=[],statStep=1,tstop=None):
        self.integrator = integrator
        self.periodicWork = periodicWork
        self.statStep = statStep # use this to override how frequently we print to the screen
        self.cycle = 0
        self.time = 0
        self.dt = 0
        self.tstop = tstop

    def checkpointPeriodicWork(self, fileName):
        """Serialize any periodicWork item's own auxiliary state (e.g. running
        accumulators) to a JSON side file. RestartWriter/RestartReader only
        cover NodeList fields and integrator cycle/time/dt -- anything a
        periodicWork callable tracks itself lives outside that C++ state graph
        and must opt in here via a checkpointState() method."""
        state = {}
        for i, work in enumerate(self.periodicWork):
            if hasattr(work, "checkpointState"):
                state[str(i)] = work.checkpointState()
        with open(fileName, "w") as f:
            json.dump(state, f)

    def restorePeriodicWork(self, fileName):
        """Counterpart to checkpointPeriodicWork: restores auxiliary state to
        any periodicWork item that defines restoreState(data), matching by
        position in the periodicWork list. Assumes periodicWork is passed in
        the same order it was at checkpoint time."""
        with open(fileName, "r") as f:
            state = json.load(f)
        for i, work in enumerate(self.periodicWork):
            key = str(i)
            if key in state and hasattr(work, "restoreState"):
                work.restoreState(state[key])
    def Step(self,nsteps=1):
        if self.integrator.Cycle() == 0:
            self.integrator.Initialize()
            if len(self.periodicWork) > 0:
                for work in self.periodicWork:
                    work(0, self.integrator.Time(), self.integrator.dt)
        for i in range(nsteps):
            self.integrator.Step()
            cycle = self.integrator.Cycle()
            time = self.integrator.Time()
            dt = self.integrator.dt
            if self.integrator.Cycle() % self.statStep == 0:
                print(colorize("Cycle:", DIM), colorize("%04d"%cycle, BRIGHT_CYAN),
                    colorize(" Time:", DIM), colorize("%03.3e"%time, BRIGHT_GREEN),
                    colorize(" dt:", DIM), colorize("%03.3e"%dt, BRIGHT_YELLOW))
            if len(self.periodicWork) > 0:
                for work in self.periodicWork:
                    if cycle % work.cycle == 0:
                        work(cycle,time,dt)
            # for package in self.integrator.getPackages():
            #     package.UpdateState()
            self.time = self.integrator.Time()
            self.cycle = self.integrator.Cycle()
            self.dt = self.integrator.dt
            if(self.tstop and self.time >=self.tstop):
                break