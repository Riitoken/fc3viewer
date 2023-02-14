README.TXT

2023-FEB-14 TUE

FC3 and FC3VIEWER 
-----------------
This git repo is the official/original source code for the FC3 API and the FC3VIEWER.

The fc3viewer is a handy windows 10 application for viewing 3D models in FC3 format.
The FC3 API is contained in the files:  fc3.h and fc3.cpp
The fc3viewer demonstrates correct usage of the API.
See the the complete 32 byte header below.

INSTALL
-------
Clone/download the repo/files.  Navigate to your fc3viewer folder and open the fc3viewer.sln with Visual Studio 2015 or later.  Newer Visual Studio versions (2017, 2019, etc.) will force you to convert the project.  The conversion may or may not convert all flags and settings.  One setting in particular is the preprocessor define _CRT_SECURE_NO_WARNINGS.  You want to verify this and define it or you'll want to uncheck the debug build setting for "Treat warnings as errors."

The fc3 source code API itself (fc3.cpp and fc3.h) is OS agnostic and platform agnostic.  It can port anywhere there is a C++ compiler.  So it is developer ready.  You can just drop the files into your current project and it should build clean.

FC3VIEWER USAGE
---------------
You can associate the .fc3 file extension with your fc3viewer.exe.
Doing this will let you click any .fc3 file and launch the viewer.
The viewer responds to:
	mouse-clicks (tilts the model) and 
	keystrokes (opens the file selection dialog)
If you load an fc3 model that does not have an image texture,
the fc3viewer will look for the image "skin.png" and use it if it exists.

WHY FC3?
--------
My motives for FC3 can be characterized as: "Fast Consistent 3D" models, with an emphasis on consistency.

Here are some of the core bullet points that I identified years ago after much frustration with other 3D formats.

I want zero headeaches, zero worry.
I want models that work out of the box (OOB).
I want first time consistency and reliability.
I want no guess work, no ambiguities.
I want as few dependencies as possible.
I want god-like control over the model from the moment it is loaded.
I want the model to load into a consistent predictable well-defined state.
I want the units size of the model to be explicitly stated (ie inches, feet, meters, etc.)
I want no issues with endianess.
I want no burden to know how the model was saved.
I want no burden to worry about how the original model was rotated.
I want no burden to tweak or manipulate the default state.
I want to easily scale the model, or rotate it if I choose.
I want an elegant way to do compression without any external dependencies.
I want the model to just work, with or without an image texture.
I want a vertex format that can model any typical 3D object.
I want a range of scale that permanently eliminates all complaints.
I want a format that can be streamed or wholistically cached.
I want an easily extensible versioning strategy.
I want to do as little as possible to load and render a model.
I want a simple light-weight file header. (see fc3.h)

LICENSE
-------
The license is called the Glory Ware License (GWL) - which is also my invention.
The GWL is a somewhat humorous tongue-in-cheek idea where you pay for a license with a "Thank you."
It could be called the Thank-You Ware License (TYWL).
In this age of social media, it seemed fitting.
The original GWL broiler plate text is a google doc here:
https://docs.google.com/document/d/1xAZ-WAHxwBuu1H-LszPGuMQXLryebBYzBdrHzV2Gom0

FC3 models
----------
There are example fc3 models included.
These models are royalty free and included with the GWL.
The most notable model is the classic OpenGL teapot.
There are 2 versions: 'a' format, and 'b' format.
You will notice a slight visual difference.
The 'a' version is lossy due to 8-bit vertex compression.
The vase and the ufo's were procedurally generated.

