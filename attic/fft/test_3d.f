/* parallel FFT/remap/pack functions - Feb 1998

   Steve Plimpton, MS 1111, Dept 9221, Sandia National Labs  (505) 845-7873
   sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level directory of the distribution.
*/

c test 3d tranpose and 3d FFT
c set precisions in test_3d.h and fft_3d.h

      program test_3d
      include "mpif.h"

#include "test_3d.h"

c common blocks insure all arrays are allocated on 8-byte boundaries

c      parameter (maxsize=524288+16)
c      parameter (maxsize=262144+16)
      parameter (maxsize=2097152+16)
      COMPLEX_DATA fft(maxsize)
      REAL_DATA xpose(maxsize),bufxpose(maxsize)

      common /blk1/ fft
      common /blk2/ xpose
      common /blk3/ bufxpose

      integer in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
      integer out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
      integer istatus(mpi_status_size)
      real*8 time

      call mpi_init(ierror)
      call mpi_comm_rank(mpi_comm_world,node,ierror)
      call mpi_comm_size(mpi_comm_world,nprocs,ierror)

c get input values

 100  call prompt(node,nx,ny,nz,iteration,
     $     iwhich,idecomp_in,idecomp_out,ipermute,
     $     nqty,memory,iscale,idebug,icheck)

c assign input and output sections of matrix to procs

      call decompose(node,nprocs,nx,ny,nz,idecomp_in,
     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi)

c      write (6,*) 'AAA',node,
c     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi

      call decompose(node,nprocs,nx,ny,nz,idecomp_out,
     $     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi)

c      write (6,*) 'BBB',node,
c     $     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi

c check if enough memory for input and output

      isize = (in_ihi-in_ilo+1) * (in_jhi-in_jlo+1) *
     $     (in_khi-in_klo+1) * nqty
      call mpi_allreduce(isize,jsize,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (jsize.ge.maxsize) then
        if (node.eq.0)
     $       write (6,*) 'ERROR: Insufficient input space'
        call exit(0)
      endif

      isize = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
     $     (out_khi-out_klo+1) * nqty
      call mpi_allreduce(isize,jsize,1,mpi_integer,mpi_max,
     $     mpi_comm_world,ierror)

      if (isize.ge.maxsize) then
        if (node.eq.0)
     $       write (6,*) 'ERROR: Insufficient output space'
        call exit(0)
      endif

c initialize matrix

      if (iwhich.eq.0) then
        call xpose_init(xpose,nqty,nx,ny,nz,
     $       in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi)
      else
        call fft_init(fft,nx,ny,nz,
     $       in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi)
      endif

c store value = 0.0 in 1 extra memory location for bounds checking

      num1 = (in_ihi-in_ilo+1) * (in_jhi-in_jlo+1) *
     $     (in_khi-in_klo+1)
      num2 = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
     $     (out_khi-out_klo+1)
      nummax = max(num1,num2)*nqty

      if (iwhich.eq.0) then
        xpose(nummax+1) = 0.0
        if (memory.eq.0) then
          bufxpose(nummax+1) = FLOAT(0.0)
        else
          bufxpose(1) = FLOAT(0.0)
        endif
      else
        fft(nummax+1) = COMPLEX(FLOAT(0.0),FLOAT(0.0))
      endif

c write out array values

      if (idebug.eq.1) then
        num = (in_ihi-in_ilo+1) * (in_jhi-in_jlo+1) *
     $       (in_khi-in_klo+1)
        if (node.gt.0) call mpi_recv(iflag,0,mpi_integer,
     $       node-1,0,mpi_comm_world,istatus,ierror)
        if (iwhich.eq.0) then
          do m = 1,num*nqty
            write (6,*) 'BEFORE',node,m,xpose(m)
          enddo
        else
          do m = 1,num
            write (6,*) 'BEFORE',node,m,fft(m)
          enddo
        endif
        if (node.lt.nprocs-1) call mpi_send(iflag,0,mpi_integer,
     $       node+1,0,mpi_comm_world,ierror)
      endif

c perform the tranpose or FFT

      if (iwhich.eq.0) then
        call xpose_perform(xpose,bufxpose,iteration,
     $       in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $       out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
     $       nqty,ipermute,memory,time)
      else
        call fft_perform(fft,iteration,nx,ny,nz,
     $       in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $       out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
     $       ipermute,iscale,time)
      endif

c write out array values
c if doing non-multiple-of-2 # of transposes,
c   compute # of values in transposed state

      if (idebug.eq.1) then
        if (mod(iteration,2).eq.1) then
          num = (out_ihi-out_ilo+1) * (out_jhi-out_jlo+1) *
     $         (out_khi-out_klo+1)
        else
          num = (in_ihi-in_ilo+1) * (in_jhi-in_jlo+1) *
     $         (in_khi-in_klo+1)
        endif
        if (node.gt.0) call mpi_recv(iflag,0,mpi_integer,
     $       node-1,0,mpi_comm_world,istatus,ierror)
        if (iwhich.eq.0) then
          do m = 1,num*nqty
            write (6,*) 'AFTER',node,m,xpose(m)
          enddo
        else
          do m = 1,num
            write (6,*) 'AFTER',node,m,fft(m)
          enddo
        endif
        if (node.lt.nprocs-1) call mpi_send(iflag,0,mpi_integer,
     $       node+1,0,mpi_comm_world,ierror)
      endif

c if requested check results for accuracy

      if (icheck.eq.1) then
        if (iwhich.eq.0) then
          call xpose_check(xpose,nx,ny,nz,
     $         iteration,node,nqty,ipermute,
     $         in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $         out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi)
        else
          call fft_check(fft,nx,ny,nz,
     $         iteration,node,iscale,
     $         in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi)
        endif
      endif

c check for memory overflow

      iflag = 0
      if (iwhich.eq.0) then
        if (xpose(nummax+1).ne.0.0) iflag = 1
      else
      endif
      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) then
        if (node.eq.0)
     $       write (6,*) 'ERROR: Overwritten data',jflag
        call exit(0)
      endif

      iflag = 0
      if (iwhich.eq.0) then
        if (memory.eq.0.and.bufxpose(nummax+1).ne.0.0) iflag = 1
        if (memory.eq.1.and.bufxpose(1).ne.0.0) iflag = 1
      else
      endif
      call mpi_allreduce(iflag,jflag,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (jflag.gt.0) then
        if (node.eq.0)
     $       write (6,*) 'ERROR: Overwritten buffer',jflag
        call exit(0)
      endif
      
c print timing

      if (node.eq.0) then
        if (iwhich.eq.0) then
          write (6,*) 'Time for',iteration,' transposes of size',
     $         nx,ny,nz,' =',time,' on',nprocs,' procs'
        else
          write (6,*) 'Time for',iteration,' FFTs of size',
     $         nx,ny,nz,' =',time,' on',nprocs,' procs'
        endif
      endif

c start over

 900  call mpi_barrier(mpi_comm_world,ierror)
      if (node.eq.0) then
        write (6,*) 'Continue (0) No (1) Yes'
        read (5,*) iyes
      endif

      call mpi_bcast(iyes,1,mpi_integer,0,mpi_comm_world,ierror)
      if (iyes.eq.1) goto 100

      end


c -----------------------------------------------------------------------
c Subroutines
c -----------------------------------------------------------------------

c -----------------------------------------------------------------------
c prompt for user inputs

      subroutine prompt(node,nx,ny,nz,iteration,
     $     iwhich,idecomp_in,idecomp_out,ipermute,
     $     nqty,memory,iscale,idebug,icheck)

      character*80 title
 900  format (a)

      include "mpif.h"

      if (node.eq.0) then
        write (6,*) 'Title'
        read (5,900) title
        write (6,900) title(1:length(title))
        write (6,*) 'Nx Ny Nz'
        read (5,*) nx,ny,nz
      endif

      call mpi_bcast(nx,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(ny,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(nz,1,mpi_integer,0,mpi_comm_world,ierror)

c exit if 0 size

      if (nx.eq.0.or.ny.eq.0.or.nz.eq.0) call exit(0)

      if (node.eq.0) then
        write (6,*) 'Iterations'
        read (5,*) iteration
        write (6,*) '(0) Transpose, (1) FFT'
        read (5,*) iwhich
        write (6,*) 'Input  decomp by (0) xyz (blocks),'
        write (6,*) '                 (1) xy, (2) yz (3) xz,'
        write (6,*) '                 (4) x, (5) y, (6) z'
        read (5,*) idecomp_in
        write (6,*) 'Output decomp by (0) xyz (blocks),'
        write (6,*) '                 (1) xy, (2) yz (3) xz,'
        write (6,*) '                 (4) x, (5) y, (6) z'
        read (5,*) idecomp_out
        write (6,*) '(0) No permute (1) one permute (2) two permutes'
        read (5,*) ipermute

        if (iwhich.eq.0) then
          write (6,*) '# of quantities per matrix element'
          read (5,*) nqty
          write (6,*) '(0) Provide memory (1) No buffer memory'
          read (5,*) memory
        else
          nqty = 1
          write (6,*) '(0) No scaling (1) Yes scaling'
          read (5,*) iscale
        endif

        write (6,*) 'Printout results (0) No (1) Yes'
        read (5,*) idebug
        write (6,*) 'Check results (0) No (1) Yes'
        read (5,*) icheck
      endif

      call mpi_bcast(iteration,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(iwhich,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(idecomp_in,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(idecomp_out,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(ipermute,1,mpi_integer,0,mpi_comm_world,ierror)

      if (iwhich.eq.0) then
        call mpi_bcast(nqty,1,mpi_integer,0,mpi_comm_world,ierror)
        call mpi_bcast(memory,1,mpi_integer,0,mpi_comm_world,ierror)
      else
        call mpi_bcast(nqty,1,mpi_integer,0,mpi_comm_world,ierror)
        call mpi_bcast(iscale,1,mpi_integer,0,mpi_comm_world,ierror)
      endif

      call mpi_bcast(idebug,1,mpi_integer,0,mpi_comm_world,ierror)
      call mpi_bcast(icheck,1,mpi_integer,0,mpi_comm_world,ierror)

      return
      end

c -----------------------------------------------------------------------
c partition array for all possibilities
c nx = global # in fast dimension
c ny = global # in mid dimension
c nz = global # in slow dimension

c idecomp = 0 = xyz (blocks)
c           1 = xy, 2 = yz, 3 = xz
c           4 = x, 5 = y, 6 = z

      subroutine decompose(node,nprocs,nx,ny,nz,
     $     idecomp,ilo,ihi,jlo,jhi,klo,khi)

      if (idecomp.eq.0) then
        call proc2grid3d(nprocs,nx,ny,nz,npx,npy,npz)
        node_x = mod(node,npx)
        node_y = mod(node/npx,npy)
        node_z = node/(npx*npy)
        ilo = node_x*nx/npx + 1
        ihi = (node_x+1)*nx/npx
        jlo = node_y*ny/npy + 1
        jhi = (node_y+1)*ny/npy
        klo = node_z*nz/npz + 1
        khi = (node_z+1)*nz/npz
      endif

      if (idecomp.eq.1) then
        call proc2grid2d(nprocs,nx,ny,npx,npy)
        node_x = mod(node,npx)
        node_y = node/npx
        ilo = node_x*nx/npx + 1
        ihi = (node_x+1)*nx/npx
        jlo = node_y*ny/npy + 1
        jhi = (node_y+1)*ny/npy
        klo = 1
        khi = nz
      endif

      if (idecomp.eq.2) then
        call proc2grid2d(nprocs,ny,nz,npy,npz)
        node_y = mod(node,npy)
        node_z = node/npy
        ilo = 1
        ihi = nx
        jlo = node_y*ny/npy + 1
        jhi = (node_y+1)*ny/npy
        klo = node_z*nz/npz + 1
        khi = (node_z+1)*nz/npz
      endif

      if (idecomp.eq.3) then
        call proc2grid2d(nprocs,nx,nz,npx,npz)
        node_x = mod(node,npx)
        node_z = node/npx
        ilo = node_x*nx/npx + 1
        ihi = (node_x+1)*nx/npx
        jlo = 1
        jhi = ny
        klo = node_z*nz/npz + 1
        khi = (node_z+1)*nz/npz
      endif

      if (idecomp.eq.4) then
        ilo = node*nx/nprocs + 1
        ihi = (node+1)*nx/nprocs
        jlo = 1
        jhi = ny
        klo = 1
        khi = nz
      endif

      if (idecomp.eq.5) then
        ilo = 1
        ihi = nx
        jlo = node*ny/nprocs + 1
        jhi = (node+1)*ny/nprocs
        klo = 1
        khi = nz
      endif

      if (idecomp.eq.6) then
        ilo = 1
        ihi = nx
        jlo = 1
        jhi = ny
        klo = node*nz/nprocs + 1
        khi = (node+1)*nz/nprocs
      endif

      return
      end

c -----------------------------------------------------------------------
c matrix values are ordered by fast to mid to slow
c nqty values for each matrix element

      subroutine xpose_init(data,nqty,nx,ny,nz,
     $     ilo,ihi,jlo,jhi,klo,khi)
      REAL_DATA data(*)

      m = 0
      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            n = (k-1)*ny*nx + (j-1)*nx + i
            do mm = 1,nqty
              m = m + 1
              data(m) = FLOAT(n)
            enddo
          enddo
        enddo
      enddo

      return
      end


c -----------------------------------------------------------------------
c matrix values are ordered by fast to mid to slow
c complex datum is +/- of matrix location

      subroutine fft_init(data,nx,ny,nz,ilo,ihi,jlo,jhi,klo,khi)
      COMPLEX_DATA data(*)

      m = 0
      do k = klo,khi
        do j = jlo,jhi
          do i = ilo,ihi
            n = (k-1)*ny*nx + (j-1)*nx + i
            m = m + 1
            data(m) = COMPLEX(FLOAT(n),FLOAT(-n))
          enddo
        enddo
      enddo
      
      return
      end

c -----------------------------------------------------------------------
c perform remaps

      subroutine xpose_perform(data,buf,iteration,
     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
     $     nqty,ipermute,memory,time)
      REAL_DATA data(*),buf(*)
      real*8 time,time1,time2
      real*8 plan1,plan2
      integer in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
      integer out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi

      include "mpif.h"

      call remap_3d_create_plan(mpi_comm_world,
     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
     $     nqty,ipermute,memory,PREC,plan1)

      if (ipermute.eq.0) then
        call remap_3d_create_plan(mpi_comm_world,
     $       out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
     $       in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $       nqty,0,memory,PREC,plan2)
      else if (ipermute.eq.1) then
        call remap_3d_create_plan(mpi_comm_world,
     $       out_jlo,out_jhi,out_klo,out_khi,out_ilo,out_ihi,
     $       in_jlo,in_jhi,in_klo,in_khi,in_ilo,in_ihi,
     $       nqty,2,memory,PREC,plan2)
      else if (ipermute.eq.2) then
        call remap_3d_create_plan(mpi_comm_world,
     $       out_klo,out_khi,out_ilo,out_ihi,out_jlo,out_jhi,
     $       in_klo,in_khi,in_ilo,in_ihi,in_jlo,in_jhi,
     $       nqty,1,memory,PREC,plan2)
      endif

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()
      
      do i = 1,iteration
        if (mod(i,2).eq.1) then
          call remap_3d(data,data,buf,plan1)
        else
          call remap_3d(data,data,buf,plan2)
        endif
      enddo
      
      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      
      time = time2 - time1

      call remap_3d_destroy_plan(plan1)
      call remap_3d_destroy_plan(plan2)
      
      return
      end


c -----------------------------------------------------------------------
c perform FFTs

      subroutine fft_perform(data,iteration,nx,ny,nz,
     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
     $     ipermute,iscale,time)
      COMPLEX_DATA data(*)
      integer iteration,nx,ny,nz
      integer in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
      integer out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi
      real*8 time,time1,time2
      real*8 plan1,plan2

      include "mpif.h"

      call fft_3d_create_plan(mpi_comm_world,nx,ny,nz,
     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
     $     iscale,ipermute,nbuf,plan1)

      if (ipermute.eq.0) then
        call fft_3d_create_plan(mpi_comm_world,nx,ny,nz,
     $       out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
     $       in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $       iscale,0,nbuf,plan2)
      else if (ipermute.eq.1) then
        call fft_3d_create_plan(mpi_comm_world,ny,nz,nx,
     $       out_jlo,out_jhi,out_klo,out_khi,out_ilo,out_ihi,
     $       in_jlo,in_jhi,in_klo,in_khi,in_ilo,in_ihi,
     $       iscale,2,nbuf,plan2)
      else if (ipermute.eq.2) then
        call fft_3d_create_plan(mpi_comm_world,nz,nx,ny,
     $       out_klo,out_khi,out_ilo,out_ihi,out_jlo,out_jhi,
     $       in_klo,in_khi,in_ilo,in_ihi,in_jlo,in_jhi,
     $       iscale,1,nbuf,plan2)
      endif

c LAMMPS-like calls

c      call fft_3d_create_plan(mpi_comm_world,nx,ny,nz,
c     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
c     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
c     $     1,0,nbuf,plan1)
c
c      call fft_3d_create_plan(mpi_comm_world,nx,ny,nz,
c     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
c     $     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi,
c     $     1,0,nbuf,plan2)

      call mpi_barrier(mpi_comm_world,ierror)
      time1 = mpi_wtime()
      
      do i = 1,iteration

        if (mod(i,2).eq.1) then
          call fft_3d(data,data,1,plan1)
        else
          call fft_3d(data,data,-1,plan2)
        endif

c LAMMPS-like calls

c        call fft_3d(data,data,1,plan1)
c        call fft_3d(data,data,-1,plan2)
c        call fft_3d(data,data,-1,plan2)
c        call fft_3d(data,data,-1,plan2)

      enddo
      
      call mpi_barrier(mpi_comm_world,ierror)
      time2 = mpi_wtime()
      
      time = time2 - time1

      call fft_3d_destroy_plan(plan1)
      call fft_3d_destroy_plan(plan2)
      
      return
      end


c -----------------------------------------------------------------------
c check remap results

      subroutine xpose_check(data,nx,ny,nz,iteration,
     $     node,nqty,ipermute,
     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi,
     $     out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi)
      REAL_DATA data(*)
      integer in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi
      integer out_ilo,out_ihi,out_jlo,out_jhi,out_klo,out_khi

      include "mpif.h"

      nerrors = 0

c data should be back in input layout, with no permutation

      if (mod(iteration,2).eq.0) then

        m = 0
        do k = in_klo,in_khi
          do j = in_jlo,in_jhi
            do i = in_ilo,in_ihi
              n = (k-1)*ny*nx + (j-1)*nx + i
              do mm = 1,nqty
                m = m + 1
                if (data(m).ne.REAL(n)) then
                  if (nerrors.eq.0)
     $                 write (6,*) 'Bad Value:',node,n,data(m)
                  nerrors = nerrors + 1
                endif
              enddo
            enddo
          enddo
        enddo

c data is in output layout, not permuted, so i is inner loop

      else if (ipermute.eq.0) then

        m = 0
        do k = out_klo,out_khi
          do j = out_jlo,out_jhi
            do i = out_ilo,out_ihi
              n = (k-1)*ny*nx + (j-1)*nx + i
              do mm = 1,nqty
                m = m + 1
                if (data(m).ne.REAL(n)) then
                  if (nerrors.eq.0)
     $                 write (6,*) 'Bad Value:',node,n,data(m)
                  nerrors = nerrors + 1
                endif
              enddo
            enddo
          enddo
        enddo

c data is in output layout, but permuted once, so loop order is i,k,j

      else if (ipermute.eq.1) then

        m = 0
        do i = out_ilo,out_ihi
          do k = out_klo,out_khi
            do j = out_jlo,out_jhi
              n = (k-1)*ny*nx + (j-1)*nx + i
              do mm = 1,nqty
                m = m + 1
                if (data(m).ne.REAL(n)) then
                  if (nerrors.eq.0)
     $                 write (6,*) 'Bad Value:',node,n,data(m)
                  nerrors = nerrors + 1
                endif
              enddo
            enddo
          enddo
        enddo

c data is in output layout, but permuted twice, so loop order is j,i,k

      else if (ipermute.eq.2) then

        m = 0
        do j = out_jlo,out_jhi
          do i = out_ilo,out_ihi
            do k = out_klo,out_khi
              n = (k-1)*ny*nx + (j-1)*nx + i
              do mm = 1,nqty
                m = m + 1
                if (data(m).ne.REAL(n)) then
                  if (nerrors.eq.0)
     $                 write (6,*) 'Bad Value:',node,n,data(m)
                  nerrors = nerrors + 1
                endif
              enddo
            enddo
          enddo
        enddo

      endif

      call mpi_allreduce(nerrors,ncounts,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (node.eq.0) then
        if (ncounts.eq.0) write (6,*) 'Total # of errors = NONE'
        if (ncounts.ne.0) write (6,*) 'Total # of errors =',ncounts
      endif

      return
      end


c -----------------------------------------------------------------------
c check FFT results

      subroutine fft_check(data,nx,ny,nz,
     $     iteration,node,iscale,
     $     in_ilo,in_ihi,in_jlo,in_jhi,in_klo,in_khi)
      COMPLEX_DATA data(*)

      include "mpif.h"

      nerrors = 0

      if (mod(iteration,2).eq.1) then
        if (node.eq.0)
     $       write (6,*) 'Cannot check oneway FFT for errors'
        return
      endif

c data should be back in input layout, with no permutation
c due to roundoff, accuracy check is to nearest int

      m = 0
      do k = in_klo,in_khi
        do j = in_jlo,in_jhi
          do i = in_ilo,in_ihi
            n = (k-1)*ny*nx + (j-1)*nx + i
            if (iscale.eq.0) n = n * (nx*ny*nz) ** (iteration/2)
            m = m + 1
            if (nint(REAL(data(m))).ne.n.or.
     $           nint(IMAG(data(m))).ne.-n) then
              if (nerrors.eq.0)
     $             write (6,*) 'Bad Value:',node,n,data(m)
              nerrors = nerrors + 1
            endif
          enddo
        enddo
      enddo

      call mpi_allreduce(nerrors,ncounts,1,mpi_integer,mpi_sum,
     $     mpi_comm_world,ierror)
      if (node.eq.0) then
        if (ncounts.eq.0) write (6,*) 'Total # of errors = NONE'
        if (ncounts.ne.0) write (6,*) 'Total # of errors =',ncounts
      endif

      return
      end

c ----------------------------------------------------------------------
c assign nprocs to 3-d nx,ny,nz grid so as to minimize surface area
c return 3-d px,py,pz grid of procs

      subroutine proc2grid3d(nprocs,nx,ny,nz,px,py,pz)
      integer nprocs
      integer nx,ny,nz
      integer px,py,pz

      integer surf,boxx,boxy,boxz
      integer bestsurf,bestboxx,bestboxy,bestboxz
      integer ipx,ipy,ipz

      bestsurf = 2 * (nx*ny + ny*nz + nz*nx)
      bestboxx = 0
      bestboxy = 0
      bestboxz = 0

c loop thru all possible factorizations of nprocs
c surf = surface area of largest proc sub-domain
c innermost if test minimizes surface area and surface/volume ratio

      ipx = 1
      do while (ipx.le.nprocs)
        if (mod(nprocs,ipx).eq.0) then
          nremain = nprocs/ipx
          ipy = 1
          do while (ipy.le.nremain)
            if (mod(nremain,ipy).eq.0) then
              ipz = nremain/ipy
              boxx = nx/ipx
              if (mod(nx,ipx).ne.0) boxx = boxx + 1
              boxy = ny/ipy
              if (mod(ny,ipy).ne.0) boxy = boxy + 1
              boxz = nz/ipz
              if (mod(nz,ipz).ne.0) boxz = boxz + 1
              surf = boxx*boxy + boxy*boxz + boxz*boxx
              if (surf.lt.bestsurf.or.(surf.eq.bestsurf.and.
     $             boxx*boxy*boxz.gt.bestboxx*bestboxy*bestboxz)) then
                bestsurf = surf
                bestboxx = boxx
                bestboxy = boxy
                bestboxz = boxz
                px = ipx
                py = ipy
                pz = ipz
              endif
            endif
            ipy = ipy + 1
          enddo
        endif
        ipx = ipx + 1
      enddo

      if (px*py*pz.ne.nprocs) then
        write (6,*) 'Bad result in proc2grid3d'
        call exit(0)
      endif

      return
      end


c ----------------------------------------------------------------------
c assign nprocs to 2d nx,ny grid so as to minimize surface area
c return 2d px,py grid of procs

      subroutine proc2grid2d(nprocs,nx,ny,px,py)
      integer nprocs
      integer nx,ny
      integer px,py

      integer surf,boxx,boxy
      integer bestsurf,bestboxx,bestboxy
      integer ipx,ipy

      bestsurf = 2 * (nx + ny)
      bestboxx = 0
      bestboxy = 0

c loop thru all possible factorizations of nprocs
c surf = surface area of largest proc sub-domain
c innermost if test minimizes surface area and surface/volume ratio

      ipx = 1
      do while (ipx.le.nprocs)
        if (mod(nprocs,ipx).eq.0) then
          ipy = nprocs/ipx
          boxx = nx/ipx
          if (mod(nx,ipx).ne.0) boxx = boxx + 1
          boxy = ny/ipy
          if (mod(ny,ipy).ne.0) boxy = boxy + 1
          surf = boxx + boxy
          if (surf.lt.bestsurf.or.(surf.eq.bestsurf.and.
     $         boxx*boxy.gt.bestboxx*bestboxy)) then
            bestsurf = surf
            bestboxx = boxx
            bestboxy = boxy
            px = ipx
            py = ipy
          endif
        endif
        ipx = ipx + 1
      enddo

      if (px*py.ne.nprocs) then
        write (6,*) 'Bad result in proc2grid2d'
        call exit(0)
      endif

      return
      end

c ----------------------------------------------------------------------
c returns the actual length of str
c backtracks from end of string skipping over spaces
c break if test into 2 to avoid compiler evaluation of str(0:0)

      integer function length(str)
      character*(*) str
      
      n = len(str)
 10   if (n.gt.0) then
        if (str(n:n).eq.' ') then
          n = n - 1
          goto 10
        endif
      endif
      length = n

      return
      end
