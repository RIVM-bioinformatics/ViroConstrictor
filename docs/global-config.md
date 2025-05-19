# Configuration options

When you first run ViroConstrictor you will be asked to set some global configuration settings.  
These settings are used every time you execute ViroConstrictor and dictate how the analysis, or even the entire application, should behave.

There are three main sections of global configuration that need to be set up:

1. The [computing settings](#setting-up-computing-settings) that ViroConstrictor uses when running the analysis.
2. The [reproducibility settings](#setting-up-reproducibility-settings) that ensure the analysis is reproducible.
3. The [updating settings](#updating-settings) for the ViroConstrictor application.

ViroConstrictor will store the information you provide during this initial setup so you won't have to provide it every time you run an analysis.

## Setting up computing settings

You will be asked to provide the information that tells ViroConstrictor *where* the actual analysis should be executed. ViroConstrictor is able to run both on your local (Linux) machine as well as on high-performance computing (HPC/grid) infrastructures.

You will, therefore, be asked in which *computing mode* you want ViroConstrictor to run.  
The options are <u>local</u> or <u>grid</u> mode.

??? Question "What is local execution mode?"
    Local execution mode means that every step of the ViroConstrictor analysis will be performed on the computer where the analysis is initiated.  
    Use this option when you're running ViroConstrictor on your own computer or if you're unsure whether your analysis environment has a grid infrastructure.

??? Question "What is grid execution mode?"
    Grid execution mode means that every step of the ViroConstrictor analysis will be performed on other computers within the computing infrastructure. This allows a large number of analysis steps to be performed in parallel.

If you choose the grid computing mode, you'll be asked a follow-up question regarding the computing queue.  
High-performance computing infrastructures are often set up with task queues for their users. With these queues, it's possible to determine which (remote) computers execute the various analysis tasks and prioritize some tasks over others.  
You will be asked to provide the name of the computing queue that you wish to use during analysis with ViroConstrictor.  
If you don't know the computing queue, please check with your system administrator beforehand.

## Setting up reproducibility settings

ViroConstrictor will attempt to automatically detect whether the use of containers is possible on your system; using containers is the preferred method for ensuring reproducibility of the analysis.  
If ViroConstrictor detects that containers can be used through Apptainer, it will automatically configure this method as the preferred approach. However, you will still be asked to provide a path where the downloaded containers can be stored on your system.  
It is recommended to use a dedicated directory for this purpose. If no path is provided, then ViroConstrictor will create and use the default path, which is the `~/.viroconstrictor/containers` directory.

!!! info "Specifying a path on a shared system"
    If you're using a shared system, as is common in high-performance computing environments, it is recommended to choose a path that is shared between all users.  
    This way, containers can be reused and all users can benefit from them, resulting in the analysis being reproducible for everyone.  
    Additionally, please ensure that the path you provide is accessible by all users and that it is not automatically cleaned up by the system.

## Updating settings

ViroConstrictor is able to update itself to newer *minor* and *patch* versions.  
During the initial configuration setup, you will be asked about your preferences regarding the auto-updating features.

When auto-updating is enabled, ViroConstrictor will update itself to the latest possible version without asking you first. In this mode, ViroConstrictor will only notify you that it *has* updated to the latest version or that it was unable to update itself.

If you choose to disable auto-updating during the configuration, you'll be asked a follow-up question about whether, instead of completely automatic updates, you wish to be prompted before updating to the latest version.  
For most use cases, this is the preferred option because it allows ViroConstrictor to update itself while still giving you control over the version being used.
