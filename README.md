HDX DATA HANDLER (2023)
=======================

HDX Data Handler is an interface package designed to:

	1. Simulate protein morphing from a given start and end state using Linear and RigiMOL.
	2. Estimate SASA and H bonds for each morph state.
	3. Graph SASA, H bonds, and experimental HDX-MS activity.

----------------------
INTERFACE HOUSEKEEPING
----------------------

1. To start, within the folder "hdx-data-handler-repo", find the following .py files:

	- ColorStatesGraph.py
	- HDX_data_handler_main.py
	- InterpolatedHBONDS.py
	- InterpolatedSASA.py
	- UpperSASAlimit_def.py

2. Open each .py file in your IDE of choice.

3. Set filepaths for the indicated variables.

	- see inline commments for formatting
	- file extensions (such as .csv) are indicated where needed
	- if no file extensions are included in the comment, do not include the extension in the path

4. Run HDX_data_handler_main.py

5. Provide all required entries regardless of chosen morphing method to minimise errors.

---------------
IMPORTANT NOTES
---------------

Calculating SASA for Linear and RigiMOL will take a while depending on your machine. 
If calculations are running, the "Calculate SASA" button should remain disabled.
Once the button returns to normal, the calculations should be done.

Morphing via RigiMOL opens PyMOL and shows the two protein states in the window.
It will appear to be unresponsive for five minutes or longer, but it is processing the morph in the background.
Do not close the PyMOL window or terminate the interface program until the RigiMOL morph has concluded.
RigiMOL will have finished processing when a moving animation between the protein states appear in PyMOL.


