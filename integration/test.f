      program main
      double precision a,b,result,abserr,resabs,resasc
      double precision book1,book3,book11,book16,book455
      double precision alpha
      double precision alist(1000),blist(1000),rlist(1000)
      double precision elist(1000)
      integer iord(1000)
      integer inf
      common /ALPHA/alpha
      external book1,book3,book11,book16,book455

      call gsl_ieee_env_setup

      a = 0.0
      b = 1.0

c      alpha = 2.6
c      print *,'alpha = ',alpha
c      call dqk15(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk21(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk31(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk41(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk51(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk61(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc

c      alpha = -0.9
c      print *,'alpha = ',alpha
c      call dqk15(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk21(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk31(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk41(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk51(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc
c      call dqk61(book1,a,b,result,abserr,resabs,resasc)
c      write(6,1) result, abserr, resabs, resasc


c     alpha = 1.3
      a = 0.3
      b = 2.71
c     print *,'OSCILLATINg alpha = ',alpha
c     call dqk15(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk21(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk31(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk41(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk51(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc
c     call dqk61(book3,a,b,result,abserr,resabs,resasc)
c     write(6,1) result, abserr, resabs, resasc



c      ier = 0
c      neval = 0
c      alpha = 2.6
c      epsabs = 1.0e-1
c      epsrel = 0.0
c      a = 0.0
c      b = 1.0
c      print *,'QNG book1'
c      print *,'alpha = ',alpha
c      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier
c      epsabs = 0.0
c      epsrel = 1e-9
c      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier
c      epsabs = 0.0
c      epsrel = 1e-13
c      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier
c      alpha = -0.9
c      epsabs = 0.0
c      epsrel = 1e-3
c      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier

c     alpha = 1.3
c      a = 0.3
c      b = 2.71
c      epsabs = 0.0
c      epsrel = 1e-12
c      call dqng(book3,a,b,epsabs,epsrel,result,abserr,neval,ier)
c      write(6,2) result, abserr, neval, ier

c      alpha = 2.6
c      a = 0.0
c      b = 1.0
c      epsabs = 0.0
c      epsrel = 1d-10
c      key = 1
c      limit = 1000
c      print *, 'DQAGE'
c      call dqage(book1,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do 10 i=1,10
c         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c 10   continue

c      alpha = 2.6
c      a = 0.0
c      b = 1.0
c      epsabs = 1d-14
c      epsrel = 0
c      key = 2
c      limit = 1000
c      print *, 'DQAGE'
c      call dqage(book1,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do 11 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c 11   continue

c     alpha = 1.3
c     a = 0.3
c     b = 2.71
c     epsabs = 1d-14
c     epsrel = 0
c     key = 3
c     limit = 1000
c     print *, 'DQAGE oscill'
c     call dqage(book3,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 12 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c12   continue

c     alpha = 2.0
c     a = -1.0
c     b =  1.0
c     epsabs = 1d-14
c     epsrel = 0
c     key = 5
c     limit = 1000
c     print *, 'DQAGE sing'
c     call dqage(book16,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 13 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c13   continue


c     alpha = 1.0
c     a = -1.0
c     b =  1.0
c     epsabs = 1d-7
c     epsrel = 0
c     key = 6
c     limit = 3
c     print *, 'DQAGE sing'
c     call dqage(book16,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 14 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c14   continue


c     alpha = 2.6
c     a = 0.0
c     b = 1.0
c     epsabs = 0.0
c     epsrel = 1d-10
c     limit = 1000
c     print *, 'DQAGSE'
c     call dqagse(book1,a,b,epsabs,epsrel,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 15 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c15   continue


c     alpha = 2.6
c     a = 0.0
c     b = 1.0
c     epsabs = 1d-14
c     epsrel = 0.0
c     limit = 1000
c     print *, 'DQAGSE abs'
c     call dqagse(book1,a,b,epsabs,epsrel,LIMIT,result,abserr,
c    $     neval,ier,alist,blist,rlist,elist,iord,last)
c     write(6,3) result, abserr, neval, ier, last
c     do 16 i=1,10
c        write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c16   continue

c      alpha = 2.0
c      a = 1.0
c      b = 1000.0
c      epsabs = 1d-7
c      epsrel = 0.0
c      limit = 1000
c      print *, 'DQAGSE abs'
c      call dqagse(book11,a,b,epsabs,epsrel,LIMIT,result,abserr,
c     $     neval,ier,alist,blist,rlist,elist,iord,last)
c      write(6,3) result, abserr, neval, ier, last
c      do i=1,10
c         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
c      enddo

      alpha = 2.0
      a = 0.0
      inf = 1
      epsabs = 0.0
      epsrel = 1.0d-3
      limit = 1000
      print *, 'DQAGSI abs'
      call dqagie(book455,a,inf,epsabs,epsrel,limit,result,abserr,
     $     neval,ier,alist,blist,rlist,elist,iord,last)
      write(6,3) result, abserr, neval, ier, last
      do i=1,20
         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
      enddo


 1    format(
     $     '-----------------',/
     $     'double exp_result =',1pe25.18, ';',/
     $     'double exp_abserr =',1pe25.18, ';',/
     $     'double exp_resabs =',1pe25.18, ';',/
     $     'double exp_resasc =',1pe25.18, ';')
 2    format(
     $     '-----------------',/
     $     'double exp_result =',1pe25.18, ';',/
     $     'double exp_abserr =',1pe25.18, ';',/
     $     'double exp_neval  =',I8, ';',/
     $     'double exp_ier    =',I8, ';')
 3    format(
     $     '-----------------',/
     $     'double exp_result =',1pe25.18, ';',/
     $     'double exp_abserr =',1pe25.18, ';',/
     $     'int    exp_neval  =',I8, ';',/
     $     'int    exp_ier    =',I8, ';',/
     $     'int    last       =',I8, ';')
 4    format('i=',i4,' a=',1pe25.18,' b=',1pe25.18,
     $     ' r=',1pe25.18,' e=',1pe25.18, ' iord=',i4)
      end




      double precision function book1(x)
      double precision alpha,x
      common /ALPHA/alpha
      book1=x**alpha*log(1.0/x)
      end

      double precision function book3(x)
      double precision alpha,x
      common /ALPHA/alpha
      book3=cos((2**alpha)*sin(x))
      end

      double precision function book11(x)
      double precision alpha,x
      common /ALPHA/alpha
      book11=(log(1/x))**(alpha-1.0)
      end

      double precision function book16(x)
      double precision alpha,x
      common /ALPHA/alpha
      book16=(x**(alpha-1.0))/((1.0+10.0*x)**2.0)
      end

      double precision function book455(x)
      double precision alpha,x
      common /ALPHA/alpha
      book455=log(x)/(1.0 + 100.0*x*x)
      end





