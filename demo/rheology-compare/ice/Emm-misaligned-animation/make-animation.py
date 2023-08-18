# N. M. Rathmann <rathmann@nbi.ku.dk>, 2023

"""
Make animation from scratch
"""

import os, PyPDF2

MAKE_FRAMES = 0
MAKE_PDF    = 0
MAKE_ANI    = 1

#-----------------------

density = 150 
framerate = 16
scale = 750
numfmt = '%03d'

fallframes = 'allframes.pdf'
fout = 'Emm-misaligned'

#-----------------------

if MAKE_FRAMES:
    os.system(r'python3 Emm-misaligned-frames.py')
    
if MAKE_PDF:
    os.system(r'pdflatex %s'%(fallframes[:-4]))
    os.system(r'rm frames/frame-*')
    os.system(r'pdfseparate %s frames/frame-%s.pdf'%(fallframes,numfmt))
    Npages = len(PyPDF2.PdfReader(open(fallframes,'rb')).pages)
    for i in range(1,Npages+1): 
        print('pdftocairo frame %i'%(i))
        os.system('pdftocairo -singlefile -png -r %i frames/frame-%03d.pdf frames/frame-%03d'%(density, i,i))

if MAKE_ANI:
    os.system('rm %s.avi %s.gif'%(fout,fout))
    os.system('ffmpeg -y -f image2 -framerate %i -stream_loop 0 -i frames/frame-%s.png -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -crf 20  -pix_fmt yuv420p %s.avi'%(framerate, numfmt, fout))
    os.system('ffmpeg -i %s.avi -vf "fps=%i,scale=%i:-1:flags=lanczos,split[s0][s1];[s0]palettegen[p];[s1][p]paletteuse" -loop 0 %s.gif'%(fout, framerate, scale, fout))

