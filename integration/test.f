      program main
      double precision a,b,result,abserr,resabs,resasc
      double precision book1,book3,book11,book16
      double precision alpha
      double precision alist(1000),blist(1000),rlist(1000)
      double precision elist(1000)
      integer iord(1000)
      common /ALPHA/alpha
      external book1,book3,book11,book16

      call gsl_ieee_env_setup

      a = 0.0
      b = 1.0
      alpha = 2.6
      print *,'alpha = ',alpha
      call dqk15(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk21(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk31(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk41(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk51(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk61(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc

      alpha = -0.9
      print *,'alpha = ',alpha
      call dqk15(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk21(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk31(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk41(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk51(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk61(book1,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc


      alpha = 1.3
      a = 0.3
      b = 2.71
      print *,'OSCILLATINg alpha = ',alpha
      call dqk15(book3,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk21(book3,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk31(book3,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk41(book3,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk51(book3,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc
      call dqk61(book3,a,b,result,abserr,resabs,resasc)
      write(6,1) result, abserr, resabs, resasc



      ier = 0
      neval = 0
      alpha = 2.6
      epsabs = 1.0e-1
      epsrel = 0.0
      print *,'QNG book1'
      print *,'alpha = ',alpha
      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
      write(6,2) result, abserr, neval, ier
      epsabs = 0.0
      epsrel = 1e-9
      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
      write(6,2) result, abserr, neval, ier
      epsabs = 0.0
      epsrel = 1e-13
      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
      write(6,2) result, abserr, neval, ier
      alpha = -0.9
      epsabs = 0.0
      epsrel = 1e-3
      call dqng(book1,a,b,epsabs,epsrel,result,abserr,neval,ier)
      write(6,2) result, abserr, neval, ier

      alpha = 1.3
      a = 0.3
      b = 2.71
      epsabs = 0.0
      epsrel = 1e-12
      call dqng(book3,a,b,epsabs,epsrel,result,abserr,neval,ier)
      write(6,2) result, abserr, neval, ier

      alpha = 2.6
      a = 0.0
      b = 1.0
      epsabs = 0.0
      epsrel = 1d-10
      key = 1
      limit = 1000
      print *, 'DQAGE'
      call dqage(book1,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
     $     neval,ier,alist,blist,rlist,elist,iord,last)
      write(6,3) result, abserr, neval, ier, last
      do 10 i=1,10
         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
 10   continue

      alpha = 2.6
      a = 0.0
      b = 1.0
      epsabs = 1d-14
      epsrel = 0
      key = 2
      limit = 1000
      print *, 'DQAGE'
      call dqage(book1,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
     $     neval,ier,alist,blist,rlist,elist,iord,last)
      write(6,3) result, abserr, neval, ier, last
      do 11 i=1,10
         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
 11   continue

      alpha = 1.3
      a = 0.3
      b = 2.71
      epsabs = 1d-14
      epsrel = 0
      key = 3
      limit = 1000
      print *, 'DQAGE oscill'
      call dqage(book3,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
     $     neval,ier,alist,blist,rlist,elist,iord,last)
      write(6,3) result, abserr, neval, ier, last
      do 12 i=1,10
         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
 12   continue

      alpha = 2.0
      a = -1.0
      b =  1.0
      epsabs = 1d-14
      epsrel = 0
      key = 5
      limit = 1000
      print *, 'DQAGE sing'
      call dqage(book16,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
     $     neval,ier,alist,blist,rlist,elist,iord,last)
      write(6,3) result, abserr, neval, ier, last
      do 13 i=1,10
         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
 13   continue


      alpha = 1.0
      a = -1.0
      b =  1.0
      epsabs = 1d-7
      epsrel = 0
      key = 6
      limit = 3
      print *, 'DQAGE sing'
      call dqage(book16,a,b,epsabs,epsrel,KEY,LIMIT,result,abserr,
     $     neval,ier,alist,blist,rlist,elist,iord,last)
      write(6,3) result, abserr, neval, ier, last
      do 14 i=1,10
         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
 14   continue


      alpha = 2.6
      a = 0.0
      b = 1.0
      epsabs = 0.0
      epsrel = 1d-10
      limit = 1000
      print *, 'DQAGSE'
      call dqagse(book1,a,b,epsabs,epsrel,LIMIT,result,abserr,
     $     neval,ier,alist,blist,rlist,elist,iord,last)
      write(6,3) result, abserr, neval, ier, last
      do 15 i=1,10
         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
 15   continue


      alpha = 2.6
      a = 0.0
      b = 1.0
      epsabs = 1d-14
      epsrel = 0.0
      limit = 1000
      print *, 'DQAGSE abs'
      call dqagse(book1,a,b,epsabs,epsrel,LIMIT,result,abserr,
     $     neval,ier,alist,blist,rlist,elist,iord,last)
      write(6,3) result, abserr, neval, ier, last
      do 16 i=1,10
         write(6,4) i,alist(i),blist(i),rlist(i),elist(i),iord(i)
 16   continue

      alpha = 2.0
      a = 1.0
      b = 1000.0
      epsabs = 1d-7
      epsrel = 0.0
      limit = 1000
      print *, 'DQAGSE abs'
      call dqagse(book11,a,b,epsabs,epsrel,LIMIT,result,abserr,
     $     neval,ier,alist,blist,rlist,elist,iord,last)
      write(6,3) result, abserr, neval, ier, last
      do i=1,10
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





