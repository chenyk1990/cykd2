from __future__ import absolute_import, print_function
import numpy
import numpy as np
from numpy import dot
import math
import scipy as sp
from bowpy.util.fkutil import nextpow2
from bowpy.util.base import stream2array, array2stream
import sys

def ssa_denoise_recon(st, p, flow, fhigh):
	"""
	SSA method, that de-noises the data given in stream by a rank reduction of the singular values of the
	Hankel matrix, created from the data in st and the sampling interval of the traces, to p.

	:param st:     Stream of data
	:type  st:

	:param dt:     sampling interval
	:type  dt:	   float

	:param p:      number of singular values used to reconstuct the data
	:type  p:	   int
	
	:param flow:   min  freq. in the data in Hz
	:type  flow:   float

	:param fhigh:  max  freq. in the data in Hz
	:type  fhigh:  float


	Example
	st = stream
	dt = st[0].stats.delta
	p = 4
	flow = 1
	fhigh = 250

	st_ssa = ssa_denoise_recon(st, dt, p, flow, fhigh)
	"""
	st_tmp = st.copy()
	
	data = stream2array(st_tmp)

	dt = st_tmp[0].stats.delta

	data_ssa = fx_ssa(data,dt,p,flow,fhigh)
	
	st_ssa = array2stream(data_ssa, st_tmp)
	
	return st_ssa

def ssa(d,nw,p,ssa_flag):
	"""
	SSA: 1D Singular Spectrum Analysis for snr enhancement

	  dp,sing,R = ssa(d,nw,p,ssa_flag);

	  IN   d:   1D time series (column)
	       nw:  view used to make the Hankel matrix
	       p:   number of singular values used to reconstuct the data
	       ssa_flag = 0 do not compute R

	  OUT  dp:  predicted (clean) data
	       R:   matrix consisting of the data predicted with
	            the first eof (R[:,0]), the second eof (R[:,1]) etc
	       sing: singular values of the Hankel matrix

	  Example:
		from math import pi
		import numpy as np
		from numpy import cos
		import matplotlib.pyplot as plt
		from bowpy.filter.ssa import ssa
		import scipy.io as sio
		
		rand =  sio.loadmat("../../mtz_ssa/randomnumbers.mat")
		r = rand['r']
		d = (cos(2*pi*0.01*np.linspace(1,200,200)) + 0.5*r[:,0])
		dp, sing, R = ssa(d,100,2,0)
		
		
		plt.plot(d/d.max())
		plt.plot(dp/dp.max()+3)
		plt.ion()
		plt.draw()
		plt.show()
		plt.ioff()

	  Based on: 

	  M.D.Sacchi, 2009, FX SSA, CSEG Annual Convention, Abstracts,392-395.
	                    http://www.geoconvention.org/2009abstracts/194.pdf

	  Copyright (C) 2008, Signal Analysis and Imaging Group.
	  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
	  Author: M.D.Sacchi
	  Translated to Python by: S. Schneider, 2016



	  This program is free software: you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published
	  by the Free Software Foundation, either version 3 of the License, or
	  any later version.

	  This program is distributed in the hope that it will be useful,
	  but WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	  GNU General Public License for more details: http://www.gnu.org/licenses/

	"""



	# Check for Data type of variables.
	if not type(d) == numpy.ndarray:
		print( "Wrong input type of d, must be numpy.ndarray" )
		raise TypeError

	nt = d.size
	N = int(nt-nw+1)
	l = np.arange(0,nw,1)
	R = np.zeros((nt,p))

	# Make Hankel Matrix.
	M = np.zeros((N-1,N)).astype('complex')
	Mp = np.zeros((N-1,N)).astype('complex')

	for k in range(N):
		M[:,k] = d[k+l]

	# Eigenimage decomposition

	U,S,V = sp.linalg.svd(M)
	


	# Reconstruct with one oscillatory component at the time.
	if not ssa_flag == 0:
		for k in range(p):
			u = np.zeros((N-1,2)).astype('complex')
			u[:,0] = U[:,k]
			Mp = dot( dot(u, u.conj().transpose()), M )
			R[:,k] = average_anti_diag(Mp)
		dp = sum(d)

	else:
		
		for k in range(p):
			u = np.zeros((N-1,2)).astype('complex')
			u[:,0] = U[:,k]
			Mp = Mp + dot( dot(u, u.conj().transpose()), M )

		R = None
		dp = average_anti_diag(Mp)
		

	sing = S

	return(dp,sing,R)

def fx_ssa(data,dt,p,flow,fhigh):
	"""
	FX_SSA: Singular Spectrum Analysis in the fx domain for snr enhancement
	
	
	 [data_f] = fx_ssa(data,dt,p,flow,fhigh);
	
	  IN   data:      data (traces are columns)
	       dt:     sampling interval
	       p:      number of singular values used to reconstuct the data
	       flow:   min  freq. in the data in Hz
	       fhigh:  max  freq. in the data in Hz
	
	
	  OUT  data_f:  filtered data
	
	  Example:
	
	        d = linear_events;
	        [df] = fx_ssa(d,0.004,4,1,120);
	        wigb([d,df]);
	
	  Based on:
	
	  M.D.Sacchi, 2009, FX SSA, CSEG Annual Convention, Abstracts,392-395.
	                    http://www.geoconvention.org/2009abstracts/194.pdf
	
	  Copyright (C) 2008, Signal Analysis and Imaging Group.
	  For more information: http://www-geo.phys.ualberta.ca/saig/SeismicLab
	  Author: M.D.Sacchi
	  Translated to Python by: S. Schneider 2016
	
	  This program is free software: you can redistribute it and/or modify
	  it under the terms of the GNU General Public License as published
	  by the Free Software Foundation, either version 3 of the License, or
	  any later version.
	
	  This program is distributed in the hope that it will be useful,
	  but WITHOUT ANY WARRANTY; without even the implied warranty of
	  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	  GNU General Public License for more details: http://www.gnu.org/licenses/
	
	"""
	nt, ntraces = data.shape
	nf = 2 * 2 ** nextpow2(nt)

	# First and last samples of the DFT.

	ilow = int(math.floor(flow*dt*nf)+1)
	if ilow < 1:
		ilow = 1

	ihigh = int(math.floor(fhigh*dt*nf)+1)
	if ihigh > math.floor(nf/2)+1:
		ihigh = int(math.floor(nf/2)+1)
	
	data_FX = np.fft.fft(data, nf, axis=0)
	data_FX_f = np.zeros(data_FX.shape).astype('complex')
	
	nw = int(math.floor(ntraces/2))

	print("		Loop through frequencies \n")
	i=1
	rend = len(range(ilow,ihigh+1))
	for k in range(ilow,ihigh+1):
		
		prcnt = 100*i/rend
		print("		%i %% done" % (prcnt), end="\r")
		sys.stdout.flush()

		tmp = data_FX[k-1,:].transpose()
		
		for j in range(10):
			tmp_out = ssa(tmp,nw,p,0)[0]
			tmp = tmp_out

		data_FX_f[k-1,:] = tmp_out
		i+=1

	for k in range(nf/2+2, nf):
		data_FX_f[k-1,:] = data_FX_f[nf-k+1,:].conj()
		
	data_f = np.fft.ifft(data_FX_f, axis=0)
	data_f = data_f[0:nt,:].real
	
	return data_f

def average_anti_diag(A):
	"""
	Given a Hankel matrix A,  this program retrieves
	the signal that was used to make the Hankel matrix
	by averaging along the antidiagonals of A.

	M.D.Sacchi
	2008
	SAIG - Physics - UofA
	msacchi@ualberta.ca


	In    A: A hankel matrix

	Out   s: signal (column vector)
	"""

	"""
	MATLAB
	[m,n] = size(A);
	N = m+n-1;

	 s = zeros(N,1);

	 for i = 1 : N

	  a = max(1,i-m+1);
	  b = min(n,i);

	   for k = a : b
	    s(i,1) = s(i,1) + A(i-k+1,k);
	   end

	 s(i,1) = s(i,1)/(b-a+1);

	 end;
	"""

	m,n = A.shape

	N = m+n-1

	s = np.zeros(N).astype('complex')

	for i in range(N):
		a = max(1,(i+1)-m+1)
		b = min(n,(i+1))
		
		if a == b:
			k = a
			s[i] = s[i] + A[i-k+1,k-1]
		else:
			for k in range(a,b+1):
				s[i] = s[i] + A[i-k+1,k-1]

		s[i]= s[i]/(b-a+1)
		
	return(s)











