      subroutine tsreset()
c   Close the fixed TISEAN i/o units between in-process runs.
c
c   readfile.f's infile/outfile open units 10 and 11 (istdio_temp.f) and never
c   close them. In a standalone program that is harmless -- process exit closes
c   them. In a MEX file the second call re-opens a unit that is still connected
c   and positioned at end-of-file, and gfortran raises
c
c      Fortran runtime error: Sequential READ or WRITE not allowed after EOF
c
c   which calls exit() from inside libgfortran and takes MATLAB with it.
c   Measured: c2t call 1 matched the standalone executable byte for byte and
c   call 2 died. With this routine, three consecutive calls are identical.
c
c   The args common block is NOT reset here: whatido() calls argdel(0) at the
c   start of every program, and touching it outside a run would read an argv
c   the gateway has already released.
c   c1 seeds with rand(sqrt(abs(fcan("I",0.0)))), and with no -I on the command
c   line that argument is 0.0 -- which SLATEC's RAND treats as "advance the
c   sequence", not "reseed". Its state (IX1,IX0) lives in a SAVE block that
c   DATA sets to 0,0 exactly once, when the image is loaded. A fresh process
c   therefore starts at 0,0 and a second in-process call does not, so the point
c   permutation in d1 differs and c1 returns different numbers. Measured: call
c   1 matched the executable byte for byte and call 2 did not.
c
c   RAND(1.0) takes the R>0 branch with MOD(1.,1.)=0, which sets IX1=IX0=0 and
c   returns without advancing -- precisely the freshly-loaded state.
      logical op
      real rand, rdummy
      external rand
c   libgfortran buffers unit 0 and flushes it when the image exits, which in a
c   MEX file is long after the gateway has restored fd 2 -- the text then
c   escapes into MATLAB's command window instead of the captured string.
c   Flushing here, before the gateway restores the descriptor, keeps it.
      call flush(0)
      rdummy = rand(1.0)
      inquire(unit=10,opened=op)
      if(op) close(10)
      inquire(unit=11,opened=op)
      if(op) close(11)
      end

      subroutine tsdie(s)
c   Print a message and unwind to the MEX gateway. The generated MEX copies of
c   the Fortran sources call this where the originals say `stop "message"`,
c   so the diagnostic survives the conversion.
      character*(*) s
      write(0,'(a)') s
      call tsfatal(1)
      end
