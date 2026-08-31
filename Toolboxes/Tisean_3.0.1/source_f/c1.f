c===========================================================================
c
c   This file is part of TISEAN
c 
c   Copyright (c) 1998-2007 Rainer Hegger, Holger Kantz, Thomas Schreiber
c 
c   TISEAN is free software; you can redistribute it and/or modify
c   it under the terms of the GNU General Public License as published by
c   the Free Software Foundation; either version 2 of the License, or
c   (at your option) any later version.
c
c   TISEAN is distributed in the hope that it will be useful,
c   but WITHOUT ANY WARRANTY; without even the implied warranty of
c   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c   GNU General Public License for more details.
c
c   You should have received a copy of the GNU General Public License
c   along with TISEAN; if not, write to the Free Software
c   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
c
c===========================================================================
c   information dimension, fixed mass
c   see  H. Kantz, T. Schreiber, Nonlinear Time Series Analysis, Cambridge
c      University Press (1997)
c   author T. Schreiber (1999)
c===========================================================================

      parameter(nx=100000,mx=10)
      dimension x(nx,mx), icol(mx)
c
c   hctsa 2.0: filename buffers raised from character*72 to character*512.
c   TISEAN declared every filename character*72 and gfortran truncates a longer
c   name silently, so the program writes to a path the caller never asked for
c   and the caller then looks for a file that does not exist. A Windows
c   tempname is already 72 characters. Only the filename variables are widened:
c   aline is a data-line buffer that is sliced with a hardcoded :72, and the
c   option-list buffers in arguments.f are the same, so both keep their size.
c   Every routine that receives these is character*(*) and adapts on its own.
c
      character*512 file, fout
      data kmax/100/, res/2./
      data iverb/1/
      external rand

      call whatido("fixed mass approach to d1 estimation",iverb)
      id=imust("d")
      mfrom=imust("m")
      mto=imust("M")
      ntmin=imust("t")
      ncmin=imust("n")
      res=fcan("#",res)
      r=rand(sqrt(abs(fcan("I",0.0))))
      kmax=ican("K",kmax)
      resl=log(2.)/res
      nmax=ican("l",nx)
      nexcl=ican("x",0)
      call columns(mc,mx,icol)
      mcmax=max(1,mc)
      isout=igetout(fout,iverb)
      if(fout.eq." ") isout=1

      call nthstring(1,file)
      call xreadfile(nmax,mcmax,nx,x,nexcl,icol,file,iverb)
c
c   hctsa 2.0 guard. d1 computes ncomp=nmax-(mt-1)*id and then
c   log(real(ncomp-2*nmin-1)). With a time separation (-t) large enough that
c   2*ntmin+1 reaches the number of points, that argument is zero or negative:
c   the logarithm returns NaN, the neighbour count k becomes meaningless and
c   the eps-doubling sweep in d1 never terminates. Measured: -d1 -m1 -M3 -t40
c   -n50 on 50 points segfaults, and hangs when built with bounds checking.
c   hctsa's NL_TISEAN_c1 applies no minimum length of its own.
c
      mtmax=(mto-1)/mcmax+1
      if(nmax-(mtmax-1)*id-2*ntmin-1.le.0) then
         write(istderr(),'(a)')
     .      "c1: series too short for the embedding and -t requested"
         write(istderr(),*) "  points=",nmax," -m=",mto," -d=",id,
     .      " -t=",ntmin
         stop
      endif
      if(file.eq."-") file="stdin"
      if(isout.eq.1) call addsuff(fout,file,"_c1")
      call outfile(fout,iunit,iverb)
      do 10 m=mfrom,mto
         write(iunit,'(4h#m= ,i5)') m
         pr=0.
C          Ben Fulcher, 2015-03-06 requires a fix to integer do loop for some compilers
         do 20 pl=log(1./(nmax-(m-1)*id)),0.,resl
            pln=pl
            call d1(nmax,mcmax,nx,x,id,m,ncmin,pr,pln,rln,ntmin,kmax)
            if(pln.eq.pr) goto 20
            it=it+1
            pr=pln
            write(iunit,*)   exp(rln), exp(pln)
 20         continue
         write(iunit,'()')
 10      write(iunit,'()')
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "-d# -m# -M# -t# -n# "//
     .   "[-## -K# -o outfile -I# -l# -x# -c#,# -V# -h] file")
      call popt("d","delay")
      call popt("m","minimal total embedding dimension")
      call popt("M","maximal total embedding dimension")
      call popt("t","minimal time separation")
      call popt("n","minimal number of center points")
      call popt("#","resolution, values per octave (2)")
      call popt("K","maximal number of neighbours (100)")
      call popt("I","seed for random numbers")
      call popt("l","number of values to be read (all)")
      call popt("x","number of values to be skipped (0)")
      call popt("c","column(s) to be read (1 or file,#)")
      call pout("file_c1")
      call pall()
      stop
      end
