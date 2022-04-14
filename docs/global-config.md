# Configuration options

When you first run ViroConstrictor you will be asked to set some global configuration settings.  
These settings are used every time you execute ViroConstrictor and tell the application *how* the analysis, or even the entire application, should behave.

There's two main sections of global configuration that need to be set up:

1. The [computing settings](#setting-up-computing-settings) which ViroConstrictor will use when running the analysis.
2. The [updating-behavior](#updating-settings) of the ViroConstrictor application.

ViroConstrictor will store the information you provide during this first setup so you won't have to provide this information every time you run an analysis.

## Setting up computing settings

You will be asked to provide the information which tells ViroConstrictor *where* the actual analysis should be executed. ViroConstrictor is able to run both on your local (linux) machine as well as High Performace Computing (HPC/grid) infrastructures.

You will therefore be asked in which *computing-mode* you want ViroConstrictor to run.  
The options are <u>local</u> or <u>grid</u> mode.

??? Question "What is local execution mode?"
    Local execution mode means that every step of the ViroConstrictor analysis will be performed on the computer where you turned on said analysis.  
    Use this option when you're running ViroConstrictor on your own computer or if you don't know if your analysis environment has a grid infrastructure.

??? Question "What is grid execution mode?"
    Grid execution mode means that every step of the ViroConstrictor analysis will be performed on other computers within the computing-infrastructure. This allows a large amount of analysis steps to be performed in parallel.

If you choose the grid computing-mode then you'll be asked a follow up question regarding the computing-queue.  
High Performance Computing infrastructures are often set up with task-queues for its users. With these queues it's possible to determine on which (remote) computers the various analysis tasks will be executed and which tasks have priority over other types of tasks.  
You will therefore be asked to provide the name of the computing-queue that you wish to use during analysis with ViroConstrictor.  
If you don't know the computing-queue then please check with your system administrator beforehand.


## Updating settings

ViroConstrictor is able to update itself to newer *minor* and *patch* versions.  
During the initial configuration setup you will be asked your preferences regarding the auto-updating features.

When auto-updating is enabled, ViroConstrictor will update itself to the latest possible version without asking you first. When this setting is enabled then ViroConstrictor will only let you know it *has* updated to the latest possible version, or when it was unable to update itself to the latest version.

If you choose to disable auto-updating during the configuration then you'll get a follow-up question if instead of completely automatic updates you wish to be asked to update to the latest version.  
For most use-cases, this is the preferred option where you can still have ViroConstrictor update itself but you stay in control of the used version.
