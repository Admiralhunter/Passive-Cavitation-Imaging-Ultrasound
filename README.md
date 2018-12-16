# Passive-Cavitation-Imaging-Ultrasound

Analyze_FrameV3 is a script that was used to analyze the passive imaging (PData file) from ultrasound testing. This analysis allowed my research leader to determine whether certain physical features were happening during the ultrasound. The code was given to me and originally only ran on the CPU and RAM and thus each frame of data took ~5 mins. I was tasked to find how to make this code run on a GPU as the code was easily to vectorize and other optimization's. After all my adjustments and optimizations I was able to get the original data set to run up to 60Hz, thus achieving real time analysis. This was achieved on a very cheap GTX 1080 and did not require the more powerful yet expensive GPU's. The program was used to acquire a new grant for my researcher to allow him to integrate the analysis program to real time analyze while collecting data.

The limitation of the script is that the GPU will run extremely slow if too much of the VRAM is used, so either reduced the frequency, x or Z dimensions. One could also just get a larger GPU VRAM if possible as well as the program scales quite quickly in size.


Frameratetest and the associated exccel file was a way to quantify how the three variables (f,x,z) affected the framerates for my program. My reseacher gave me the task to analyze the code to help him make an informed decision into what parameters we wanted those variables to be. I automated the entire process and had it analyze over 5920 variations of those parameters to understand how they affect the framerate.
