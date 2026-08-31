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
c   c2t.f
c   Takens' estimator from c2
c   author T. Schreiber (1998)
c===========================================================================

      parameter(meps=1000)
      dimension e(meps), c(meps), lw(meps)
      double precision a, b
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
      character*72 aline
      data iverb/1/

      call whatido("Takens' estimator from correlation sum data",iverb)
      isout=igetout(fout,iverb)
      if(nstrings().eq.0) call usage()
      call nthstring(1,file)
      call infile(file,iunit,iverb)
      if(isout.eq.1) call addsuff(fout,file,"_t")
      call outfile(fout,iunit2,iverb)
 1    read(iunit,'(a)',end=999) aline
 4    if(aline(1:1).ne."#") goto 1
      if(aline(1:1).eq."#") 
     .   read(aline(index(aline,"m=")+2:72),'(i20)',err=1) m
      me=0
 2    read(iunit,'(a)') aline
      if(aline(1:72).eq." ") goto 3
      read(aline,*,err=999,end=999) ee, cc
      if(cc.le.0.) goto 3
      me=me+1
      e(me)=log(ee)
      c(me)=log(cc)
      goto 2
 3    write(iunit2,'(4h#m= ,i5)') m
      call indexx(me,e,lw)
      call index2sort(me,e,lw)
      call index2sort(me,c,lw)
      cint=0
      do 10 i=2,me
         b=(e(i)*c(i-1)-e(i-1)*c(i))/(e(i)-e(i-1))
         a=(c(i)-c(i-1))/(e(i)-e(i-1))
         if(a.ne.0) then
            cint=cint+(exp(b)/a)*(exp(a*e(i))-exp(a*e(i-1)))
         else
            cint=cint+exp(b)*(e(i)-e(i-1))
         endif
 10      write(iunit2,*) exp(e(i)), exp(c(i))/cint
      write(iunit2,'()') 
      write(iunit2,'()') 
      goto 4
 999  stop
      end

      subroutine usage()
c usage message

      call whatineed(
     .   "[-o outfile -V# -h] file")
      call pout("file_t")
      call pall()
      stop
      end
