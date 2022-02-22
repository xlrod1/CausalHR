      INTEGER FUNCTION atpos(v, n, x)

C     Determines the position of the scalar x, within the vector v

      INTEGER n
      DOUBLE PRECISION v(n), x

      INTEGER i

      IF (x .LT. v(1)) THEN
         atpos = 0
         RETURN
      END IF

      IF (x .GT. v(n)) THEN
         atpos = n
         RETURN
      END IF

      DO 100, i=1, n

         IF (x-v(i) .GE. 0.0d0) THEN
            atpos = i
         END IF

 100  CONTINUE

      RETURN
      END
      SUBROUTINE bsmoth(n, z, bopt, m, zz, bsm, b, ks, kflag,
     *                   endl, endr)

C     Computes the smoothed bandwidths at zz, using optimal bandwidths
C     at z, and smoothing bandwidth b

      INTEGER n, m, ks, kflag
      DOUBLE PRECISION  z(n), bopt(n), zz(m), bsm(m), b, endl, endr

      INTEGER ilo, ihi, i, j
      DOUBLE PRECISION z0, sum1, sum2, ker, kernel, q, u, ZERO, ONE
      PARAMETER (ZERO=0.0d0, ONE=1.0d0)

      EXTERNAL ibnds, kernel

      DO 200, j=1, m
         z0 = zz(j)

C --- compute for "active" points only

         CALL ibnds(z, n, z0, b, ilo, ihi)

         sum1 = ZERO
         sum2 = ZERO

         DO 100 i=ilo, ihi

C --- u in [-1, 1], for any i in [ilo, ihi] (this is what ibnds does)

            u = (z0-z(i)) / b
         
            IF ( (kflag.EQ.0) .OR.
     *           ((endl+b.LE.z0) .AND. (z0.LE.endr-b)) ) THEN


C --- INTERIOR point; no boundary correction

               ker = kernel(ONE, u, ks)

            ELSE IF ((endl.LE.z0) .AND. (z0.LT.endl+b)) THEN

C     LEFT boundary

               q = (z0-endl) / b

               ker = kernel(q, u, ks)

            ELSE

C     RIGHT boundary

               q = (endr-z0) / b

               ker = kernel(q, -u, ks)

            END IF

            sum1 = sum1 + ker*bopt(i)
            sum2 = sum2 + ker

 100     CONTINUE

         bsm(j) = sum1 / sum2

 200  CONTINUE

      RETURN
      END
      SUBROUTINE func(n, ks, x, delta, z, b, endl, endr,
     *                q, y, bb, vv, bpilot, kflag, w, naes, nw)

      INTEGER n, delta(n), ks, kflag
      DOUBLE PRECISION q, z, y, negy, bpilot, b , zz
      DOUBLE PRECISION endl, endr, newz, fz, k, x(n), bb, vv,
     *                 hazden, kernel, surfct, w(n), naes(n),
     *                 nw
      EXTERNAL hazden, kernel, surfct
    
      newz = z - b*y

      fz = hazden(n, ks, x, delta, newz, bpilot, endl, endr, 
     *                kflag, naes, nw)

      negy = y

      IF ( (endr-b.LT.z) .AND. (z.LE.endr) ) THEN
         negy = -y
      END IF

      k = kernel(q, negy, ks)

      bb = fz * k

      zz = surfct(x, delta, n, newz, w, nw)

      vv = k*k*fz / zz

      RETURN
      END
      DOUBLE PRECISION FUNCTION gets(x, n, x0)

C Computes the value of the survival function at x0
C Uses BINARY search

C     x --> table (matrix) computed using Kaplan-Meier
C           1st column: survival times
C           2nd column: survival function values for the corresp. time
C     n --> number of unique survival times
C     x0 --> survival time for which the survival value is desired

      INTEGER n
      DOUBLE PRECISION x(20000,2), x0

      INTEGER ilo, ihi, ihf

      IF (x0 .LT. x(1,1)) THEN
         gets = 1.0d0
      ELSE IF (x0 .GE. x(n,1)) THEN
         gets = x(n,2)
      ELSE

C     binary search (halving at each step)

         ilo = 1
         ihi = n

 100     CONTINUE

         IF (ihi-ilo .EQ. 1) THEN
            gets = x(ilo,2)
            RETURN
         END IF

         ihf = (ilo+ihi) / 2

         IF (x(ihf,1) .LT. x0) THEN
            ilo = ihf
            GO TO 100
         ELSE IF (x(ihf,1) .GT. x0) THEN
            ihi = ihf
            GOTO 100
         ELSE
            gets = x(ihf,2)
         END IF

      END IF

      RETURN
      END
      SUBROUTINE glmin (n, x, delta, ks, z, gridz, bw, gridb, endl,
     *                  endr, bpilot, imsemn, globlb, glmse, kflag,
     *                  w, naes, nw)

C Computes optimal GLOBAL bandwidth, by minimizing the 
C Integrated Mean Squared Error (IMSE).

C                         Algorithm

C For each bandwidth, compute IMSE.  Take as OPTIMAL bw the one yielding
C the smallest IMSE.  If optimal bw is the first, i.e. bw(1), it means
C that no minimum was reached.  In this case, set optimal global bw to
C the largest, i.e. bw(gridb)

      INTEGER n, delta(n), ks, gridz, gridb, kflag
      DOUBLE PRECISION x(n), z(gridz), bw(gridb), endl, endr, bpilot,
     *     imsemn, globlb, glmse(gridb), w(n), naes(n), nw

      INTEGER i, j
      DOUBLE PRECISION imse, mse, bias, var, ZERO, HUGE, hpilot(1000)
      PARAMETER (ZERO=0.0d0, HUGE=1.0d30)
      COMMON /hazpil/ hpilot

      imsemn = HUGE
      globlb = bw(gridb)

      DO 100, j=1, gridb
         imse = ZERO
         
         DO 50, i=1, gridz

            CALL msemse(n, ks, z(i), endl, endr, x, delta,
     *                  bw(j), mse, bias, var, bpilot, hpilot(i),
     *                  kflag, w, naes, nw)
            
            imse = imse + mse
 50      CONTINUE
         
         IF ( (imse.GT.ZERO) .AND. (imse.LT.imsemn) ) THEN
            imsemn = imse
            globlb = bw(j)
         END IF

         glmse(j) = imse
 100  CONTINUE
      
      RETURN
      END
      DOUBLE PRECISION FUNCTION hazden(n, ks, x, delta, z, b,
     *                                  endl, endr, kflag, naes,
     *                                 nw)       

      INTEGER n, ks, delta(n), kflag
      DOUBLE PRECISION x(n), z, b, endl, endr, naes(n)

      INTEGER i, ilo, ihi
      DOUBLE PRECISION q, u, ker, kernel, ZERO, ONE, nw
      PARAMETER (ZERO=0.0d0, ONE=1.0d0)
      EXTERNAL kernel, ibnds

C compute indices of "active" points

      CALL ibnds(x, n, z, b, ilo, ihi)

      hazden = ZERO

      DO 100, i=ilo, ihi

         IF (delta(i) .EQ. 1) THEN
            u = (z-x(i)) / b

C --- u in [-1, 1], for any i in [ilo, ihi] (this is what ibnds does)

            IF ( (kflag.EQ.0) .OR. 
     *           ((endl+b .LE.z) .AND. (z.LE.endr-b)) ) THEN

C --- INTERIOR point; no correction

               ker = kernel(ONE, u, ks)

             ELSE IF ( (endl.LE.z) .AND. (z.LT.endl+b) ) THEN

C --- LEFT boundary correction

               q = (z-endl) / b

               ker = kernel(q, u, ks)
            ELSE

C --- RIGHT boundary correction

               IF (kflag .EQ. 1) THEN
C --- LEFT only boundary correction; like INTERIOR
                  ker = kernel(ONE, u, ks)
               ELSE
                  q = (endr-z) / b
                  
                  IF ( -q .LE. u) THEN
                     ker = kernel(q, -u, ks)
                  ELSE
                     GOTO 100
                  END IF

               END IF

            END IF

            hazden = hazden + ker*naes(i)

         END IF

 100  CONTINUE

      hazden = hazden / b

C --- do not allow negative hazards

      IF ( hazden .LT. ZERO) hazden = ZERO

      RETURN
      END
      SUBROUTINE ibnds(x, n, z, b, ilo, ihi)

C Ensures that all components of x, with indices ilo...ihi,
C are within (z-b, z+b) ("active points", i.e. yielding weights
C in the kernel smoothing)

C     x --> vector of observation points
C     n --> length of x
C     z --> grid point for which "active" observation points are computed
C     b --> bandwidth; "active" points are those in [z-b, z+b]
C     ilo <-- index of the first "active" observation point
C     ihi <-- index of the last "active" observation point

      INTEGER n, ilo, ihi
      DOUBLE PRECISION x(n), z, b

      INTEGER i
      DOUBLE PRECISION what

      what = z - b

      DO 100, i=1, n
      
         IF (what .LT. x(i)) THEN
            ilo = i
            GO TO 150
         END IF

 100  CONTINUE

      ilo = n + 1

 150  what = z + b

      IF (what .GE. x(n)) THEN
         ihi = n
         RETURN
      END IF

      DO 200, i=n, ilo, -1

         IF (what .GT. x(i)) THEN
            ihi = i
            RETURN
         END IF

 200  CONTINUE

      ihi = 0

      RETURN
      END
      SUBROUTINE intgrl(n, ks, x, delta, z, b, endl, endr,
     *                    q, r, s, valueb, valuev, bpilot, 
     *                    kflag, w, naes, nw)
      INTEGER n, ks, delta(n), kflag
      DOUBLE PRECISION x(n), z, b, endl, endr, q, r, s, valueb, valuev,
     *                 bpilot,nw
      
      INTEGER j, JMAX
      DOUBLE PRECISION  oldb, oldv, EPSI, HUGE, w(n), naes(n)
      PARAMETER (JMAX=6, EPSI=1.0d-3, HUGE=1.0d30)
      INTRINSIC ABS

      oldb = -HUGE
      oldv = -HUGE

      DO 10, j=1, JMAX

         CALL try(n, ks, x, delta, z, b, endl, endr, q, r, s,
     *            valueb, valuev, j, bpilot, kflag, w,
     *            naes, nw)

         IF (ABS(valueb-oldb) .LE. EPSI*ABS(oldb) .AND.
     *       ABS(valuev-oldv) .LE. EPSI*ABS(oldv) ) RETURN

         oldb = valueb
         oldv = valuev
10    CONTINUE

      RETURN
      END
      SUBROUTINE kapmei(times, delta, n, x, count)

C Computes Kaplan-Meier estimates

      INTEGER n, delta(n), count
      DOUBLE PRECISION times(n), x(20000,2)

      INTEGER i, j, equals, lsteql, atrisk, events
      DOUBLE PRECISION prob, ONE
      PARAMETER (ONE=1.0d0)

      INTRINSIC DBLE

      atrisk = n
      count = 0
      lsteql = 0
      prob = ONE

C Loop over all observations

      i = 1

 10   IF (i .GE. n) GO TO 200

C     find how many equal times

      equals = 0
      events = delta(i)

      DO 100, j=i+1, n
         IF (times(j) .EQ. times(i)) THEN
            equals = equals + 1
            events = events + delta(j)
         ELSE
            GO TO 150
         END IF
 100  CONTINUE

 150  CONTINUE

      count = count + 1

      atrisk = atrisk - lsteql
      lsteql = equals + 1

      x(count,1) = times(i)
      x(count,2) = prob*(ONE-DBLE(events)/DBLE(atrisk))
      prob = x(count,2)

      i = i + lsteql

      GO TO 10

 200  CONTINUE

      RETURN
      END
      DOUBLE PRECISION FUNCTION kernel(q, x, ks)

C Computes the value of the modified kernel K(q,x)
C Implements Table.1 from pg. 63

      INTEGER ks
      DOUBLE PRECISION q, x

      DOUBLE PRECISION v1, v2, HALF, ONE, TWO
      PARAMETER (ONE=1.0d0, TWO=2.0d0, HALF=ONE/TWO)

      IF (ks .EQ. 0) THEN

C --- Rectangle kernel

         IF (q .EQ. ONE) THEN

C ------ Interior

            v1 = ONE
            v2 = HALF
         ELSE
            v1 = TWO / (ONE+q)**3
            v2 = TWO*(ONE-q+q*q) + 3.0d0*(ONE-q)*x
         END IF

      ELSE if (ks .EQ. 1) THEN
         
C --- Epanechnikov

         IF (q .EQ. ONE) THEN
            v1 = 0.75d0
            v2 = ONE - x*x
         ELSE
            v1 = 12.0d0*(x+ONE)/(ONE+q)**4
            v2 = HALF*(3.0d0*q*q-q-q+ONE) + x*(ONE-q-q)
         END IF

      ELSE IF (ks .EQ. 2) THEN

C --- Biquadratic

         IF (q .EQ. ONE) THEN
            v1 = 15.0d0 * (ONE-x*x)
            v2 = (ONE-x*x) / 16.0d0
         ELSE
            v1 = 60.0d0*(x+ONE)*(x+ONE)*(q-x)/(ONE+q)**6
            v2 = TWO*q*q-q-q+ONE + x*(TWO-3.0d0*q)
         END IF

      ELSE IF (ks .EQ. 3) THEN

C --- Triquadratic
         
         IF (q .EQ. ONE) THEN
            v1 = 35.0d0*(ONE-x*x)
            v2 = (ONE-x*x)*(ONE-x*x) / 32.0d0
         ELSE
            v1 = 280.0d0*(ONE+x)**3*(q-x)*(q-x)/(ONE+q)**8
            v2 = HALF*(5.0d0*q*q-6.0d0*q+3.0d0) + x*(3.0d0-4.0d0*q)
         END IF

      END IF

      kernel = v1 * v2

      RETURN
      END
      SUBROUTINE knncen(times, status, n, z, nz, k, bw)

c     k-th nearest neighbor for survival data with censorship

c     z  --> grid points where bandwidths are calculated
c     nz --> number of grid points

c     For each grid point z0 in z, computes the bandwidth bw,
c     so that in each [z0-bw, z0+bw] there are k neighbors
c     (uncensored observations, i.e. status=1)

      INTEGER n, status(n), nz, k
      DOUBLE PRECISION times(n), z(nz), bw(nz)

      INTEGER i, j, iv, ilo, ihi, ipos, count, atpos
      DOUBLE PRECISION z0, tcopy(20000), td(20000)
      EXTERNAL atpos, sorter
      INTRINSIC ABS, MIN, MAX

C     "Clean" the survival times vector,
C     i.e. eliminate all censored observations

      count = 0

      DO 100, i=1, n
         
         IF (status(i) .EQ. 1) THEN

C --- UNCENSORED observation

            count = count + 1
            tcopy(count) = times(i)
         END IF

 100  CONTINUE

C     Compute bandwidth for each grid point z, so that in each vicinity
C     of z, there are k survival times (uncensored)

      DO 200, i=1, nz
         z0 = z(i)

         ipos = atpos(tcopy, count, z0)

         ilo = MAX(ipos-k, 1)
         ihi = MIN(ipos+k, count)

         iv = 0
         DO 150, j=ilo, ihi
            iv = iv + 1
            td(iv) = ABS(tcopy(j)-z0)
 150     CONTINUE

         CALL sorter(td, iv)

         bw(i) = td(k)
 200  CONTINUE

      RETURN
      END
      SUBROUTINE knnhad (n, x, delta, ks, bwchoi,
     *       gridz, z, m, zz, bpilot, endl, endr, bsmo, kflag, fzz,
     *       kmin, kmax, bopt, bopt1, kimse, w, naes, nw) 

C     This procedure provides hazard function estimates
C     with boundary modifications and local/global bandwidth choice,
C     for censored and uncensored case.

C Nearest methods algorithms to compute the bandwidth to be used with
C modified kernel polynomials
C
C      PARAMETERS:
C
C      n --> number of observations
C      x --> vector of survival times
C      delta --> censoring indicator for hazard estimation
C                   0 - censored observation
C                   1 - uncensored observation
C      ks --> boundary kernel type; see (5.3) pg. 70
C                   0 - rectangle
C                   1 - Epanechnikov
C                   2 - biquadratic
C                   3 - triquadratic
C      bwchoi --> method to determine the bandwidth
C                   1 - nearest neighbors, eliminating censored observations
C                   2 - nearest neighbors, with censored observations
C      gridz --> number of points in the minimization grid z(gridz)
C      z --> the minimization grid to determine optimal bw.
C            gridz equidistant points between startz and endz;
C            Number of points gridz influences computing time strongly.
C            Usually: gridz < m
C      m --> number of points in output (estimation) grid
C      zz <-- output (estimation) grid, computed as m equidistant points
C             between startz and endz at which curve is estimated
C      bpilot --> initial bandwidth which is used to  estimate bias and var
C      endl --> assumed left endpoint of function to be estimated.
C               Boundary kernels are used in [endl,endl+b).
C      endr --> assumed right endpoint of function to be estimated.
C               Boundary kernels are used in (endr-b,endr],
C               b being the bandwidth.
C      bsmo --> bandwidth for smoothing of local bandwidths,
C               Small value recommended for local bandwidth choice.
C               Specification not necessary for global bandwidth choice.
C      kflag --> 0 - unmodified kernels (Epanechnikov)
C                1 - boundary corrected kernels
C      fzz <--  function estimate at zz(i), i=1, m
C      kmin <-> minimum number of neighbors.
C               OUTPUT: optimum number of neighbors
C      kmax --> maximum number of neighbors.
C      bopt <-- chosen bandwidth at each z(j), j=1, gridz, for local choice.

c               Note that final bandwidth is used on grid zz(m)
c               and is obtained from bopt(gridz) by smoothing
c               with bandwidth bsmo.

C      bopt1 <-- smoothed bandwidth used at output point zz(i)
C      kimse <-- IMSE at each k
C
      INTEGER n, delta(n), gridb, gridz, ks, kflag, kmin, kmax,
     *        bwchoi, m, k
      DOUBLE PRECISION z(gridz), zz(m), fzz(m), bopt(gridz), x(n),
     *     endr, endl, bsmo, bopt1(m), kimse(1), bpilot,
     *     hpilot(1000), w(n), naes(n),nw

      INTEGER   i
      DOUBLE PRECISION hazden
      COMMON /hazpil/ hpilot

      EXTERNAL knnmin, olafmn, bsmoth, hazden

C --- Compute the hazard estimates for bw=bpilot.
C     They will be used in msemse

      DO 10, i=1, gridz
         hpilot(i) = hazden(n, ks, x, delta, z(i), bpilot,
     *                   endl, endr, kflag, naes,nw)
 10   CONTINUE

      IF (bwchoi .EQ. 1) THEN

C     --- simple nearest neighbor approach

         CALL knnmin(x, delta, n, z, gridz, ks, endl, endr,       
     *               bpilot, bopt, kmin, kmax, kimse, kflag, w,
     *               naes,nw)

      ELSE IF (bwchoi .EQ. 2) THEN
         
C     --- modified nearest neighbor approach

         CALL olafmn(x, delta, n, z, gridz, ks, endl,endr,       
     *                bpilot, bopt, kmin, kmax, kimse, kflag,
     *                w, naes,nw)

      ELSE

C     for future algorithms ...

      END IF

C       --- smooth the local bandwidths: BOPT ==> BOPT1(i)

      CALL bsmoth(gridz, z, bopt, m, zz, bopt1, bsmo, ks, kflag,
     *             endl, endr)     

C     compute the estimates, using optimal bandwidths

      DO 40 i=1, m

         fzz(i) = hazden(n, ks, x, delta, zz(i), bopt1(i),
     *                    endl, endr, kflag, naes,nw)

 40   CONTINUE

      RETURN
      END
      SUBROUTINE knnmin(x, delta, n, z, gridz, ks, endl,
     *                  endr, bpilot, bopt, kmin, kmax, kimse,
     *                  kflag, w, naes,nw)

C     Computes the bandwidth at each grid point z
C     First it finds optimum number of neighbors, minimizing the IMSE

      INTEGER n, delta(n), gridz, ks, kmin, kmax, kflag
      DOUBLE PRECISION x(n), z(gridz), bopt(gridz), bpilot,
     *                 endl, endr, kimse(1),nw

      INTEGER k, i, kopt
      DOUBLE PRECISION imse, imsemn, bias, var, mse, bwi, zi, 
     *     hpilot(1000), ZERO, HUGE, w(n), naes(n)
      PARAMETER (ZERO=0.0d0, HUGE=1.0d5)
      COMMON /hazpil/ hpilot

      EXTERNAL knncen, msemse 

      IF (kmin .EQ. kmax) THEN
         CALL knncen (x, delta, n, z, gridz, kmin, bopt)

         RETURN
      END IF

      imsemn = HUGE

C     For the moment, the maximum number of neighbors to be considered
C     is half of the non-censored observations

      DO 100, k=kmin, kmax

C     compute the bandwidths bopt for k neighbors

         CALL knncen (x, delta, n, z, gridz, k, bopt)

         imse = ZERO

C     compute MSE at each gridpoint

         DO 50, i=1, gridz
            zi = z(i)
            bwi = bopt(i)

            CALL msemse(n, ks, zi, endl, endr, x, delta,
     *                  bwi, mse, bias, var, bpilot, hpilot(i),
     *                  kflag, w, naes, nw)

            imse = imse + mse
 50      CONTINUE

         IF (imse .LT. imsemn) THEN
            kopt = k
            imsemn = imse
         END IF

         kimse(k-kmin+1) = imse
 100  CONTINUE

C     kopt is returned in kmin 

      kmin = kopt

C     compute the bandwidths for kopt

      CALL knncen (x, delta, n, z, gridz, kopt, bopt)

      RETURN
      END
      SUBROUTINE loclmn (n, x, delta, ks, z, gridz, bw, gridb, bopt,
     *     endl, endr, bpilot, msemin, biasmn, varmin, kflag, 
     *     w, naes, nw)

C Computes optimal LOCAL bandwidths at each gridpoint of z.

C                         Algorithm

C For each gridpoint in z, computes MSE for each bandwidth in the grid.
C The OPTIMAL bandwidth is the one yielding the smallest MSE.
C NOTE: If no minimum was found for MSE, then set bopt to largest bw

      INTEGER   n, delta(n), gridb, gridz, ks, kflag
      DOUBLE PRECISION z(gridz), x(n), bw(gridb)
      DOUBLE PRECISION bopt(gridz), endl, endr, bpilot, mse, bias, var
      DOUBLE PRECISION msemin(gridz), biasmn(gridz), varmin(gridz),
     *                 amin, hpilot(1000), ZERO, HUGE, w(n),
     *                 naes(n),nw
      PARAMETER (ZERO=0.0d0, HUGE=1.0d30)
      COMMON /hazpil/ hpilot

      INTEGER i, j
      EXTERNAL msemse

      DO 20, i=1, gridz

         amin = HUGE
         bopt(i) = bw(gridb)

         DO 10, j=1, gridb

            CALL msemse(n, ks, z(i), endl, endr, x, delta,
     *                  bw(j), mse, bias, var, bpilot, hpilot(i),
     *                  kflag, w, naes,nw)
            
            IF ( (mse.GT.ZERO) .AND. (mse .LT. amin) ) THEN
               amin = mse
               bopt(i) = bw(j)
               biasmn(i) = bias
               varmin(i) = var
            END IF

 10      CONTINUE

         IF (amin .EQ. HUGE) THEN
C     --- all mse were 0.  Set bopt to the largest bw
            msemin(i) = ZERO
         END IF

         msemin(i) = amin
 20   CONTINUE

      RETURN
      END
      SUBROUTINE locolf(x, nobs, xgrid, ngrid, n, k, bw)

C     Computes array of bandwidths at each grid point in xgrid

C     x --> matrix generated by a Kaplan-Meier survival estimation
C           1st column: unique survival times
C           2nd column: corresponding survival function values
C     nobs --> initial number of observations (censored and uncensored)
C     xgrid --> vector of survival times, where the bandwidths are
C               to be computed
C     ngrid --> number of grid points
C     n --> number of unique survival times
C     k --> number of neighbors
C     bw <-- vector at bandwidths, computed at each grid point

      INTEGER nobs, n, ngrid, k
      DOUBLE PRECISION x(20000,2), xgrid(ngrid), bw(ngrid)

      INTEGER i
      DOUBLE PRECISION oneolf
      EXTERNAL oneolf

      DO 100, i=1, ngrid
         bw(i) = oneolf(x, n, xgrid(i), nobs, k)
 100  CONTINUE

      RETURN
      END
      SUBROUTINE msemse(n, ks, z, endl, endr, x, delta, b, mse,       
     *                  bias, var, bpilot, fz, kflag, w, naes,nw)

C Computes MSE at z, for bw=b, using (2.4) from Mueller, pg. 64

      INTEGER n, delta(n), ks, kflag
      DOUBLE PRECISION z, endl, endr, x(n), b, mse, bias, var,
     *                 bpilot, fz,w(n),naes(n), nw

      DOUBLE PRECISION r, s, q, valueb, valuev, hazden, ONE
      PARAMETER (ONE=1.0d0)

      IF ((kflag.EQ.0) .OR. (endl+b.LE.z .AND. z.LE.endr-b)) THEN

C     --- INTERIOR point; no correction

         q = ONE
         r = -ONE
         s = ONE
      ELSE IF (endl.LE.z .AND. z.LT.endl+b) THEN

C     --- LEFT boundary correction

         q = (z-endl) / b
         r = -ONE
         s = q
      ELSE

C     --- RIGHT boundary correction

         IF (kflag .EQ. 1) THEN

C--- actually, only LEFT boundary correction; this is like INTERIOR

            q = ONE
            r = -ONE
            s = ONE
         ELSE

C --- indeed, RIGHT boundary correction
            
            q = (endr-z) / b
            r = -q
            s = ONE
         END IF

      END IF

      CALL intgrl(n, ks, x, delta, z, b, endl,endr, q, r, s,      
     *               valueb, valuev, bpilot, kflag, w, naes,
     *           nw)
C     nw=SUM(w)
      bias = valueb - fz
      var = DBLE(valuev/nw/b)
      mse = bias*bias + var

      RETURN
      END
      SUBROUTINE newhad(n, x, delta, ks, local, z, gridz,
     *                    zz, m, bpilot, bw, gridb, endl, endr,
     *                    bsmo, kflag, fzz, bopt, bopt1, msemin,
     *                    biasmn, varmin, imsemn, globlb,
     *                    glmse,w,naes,nw)

C     This procedure provides hazard function estimates
C     with boundary modifications and local/global
c     bandwidth choice, for censored and uncensored case.

C        INPUT DATA NEED NOT BE ORDERED

C        ** VERSION JULY 92
C        ** MODIFIED DEC 93
C        ** NO RESPONSIBILITY IS ASSUMED FOR CORRECTNESS OF CODE
C        ** COPYRIGHT H.G.MUELLER & J.L.WANG, DIVISION OF
C        ** STATISTICS, UNIVERSITY OF CALIFORNIA, DAVIS, CA 95616 USA
C        ** PROCEDURE DESCRIBED IN:
C        **  MUELLER,H.G.,WANG,J.L.: HAZARD RATE ESTIMATION UNDER RANDOM
C        **  CENSORING WITH VARYING KERNELS AND BANDWIDTHS,
C        **  BIOMETRICS 50, 61-76, 1994.
C

C ***** Modified: Dan M. Serachitopol, Nov 1997

C      PARAMETERS:
C
C      n --> number of observations
C      x --> vector of observations (survival times)
C      delta --> censoring vector for the observations
C                   0 -  censored observation
C                   1 -  uncensored observation
C      ks --> boundary kernel type; see (5.3) pg. 70
C                   0 - rectangle
C                   1 - Epanechnikov
C                   2 - biquadratic
C                   3 - triquadratic
C      local -->  local or global bandwidth choice by minimizing
c                 direct convolution type estimates of mse/imse.
c                   0 - global minimizer
c                   1 - local minimizer
C      z --> the minimization grid
C      gridz --> number of points for the minimization grid.
C                NOTE: gridz influences computing time strongly.
C                      Usually gridz < m
C      zz --> estimation grid grid, where the hazards are estimated
C      m --> number of points in the estimation grid
C      bpilot --> initial bandwidth which is used to  estimate bias/var
C      bw    --> the bandwidths grid, scanned in the MSE minimization
C      gridb --> number of bandwidths in the grid
C                NOTE: If gridb=1, startb is used as a global optimal
C                      bandwidth to compute the hazard estimates
C      endl, endr --> bounds of the  function to be estimated.
C                     Corrected boundary kernels are used in:
C                     LEFT  boundary: [endl, endl+b]
C                     RIGHT boundary: [endr-b, endr]
C      bsmo --> bandwidth for smoothing the local optimal bandwidths.
C               Small value recommended for local bandwidth choice.
C               Specification not necessary for global bandwidth choice.
C      kflag --> 0 - no boundary correction
C                1 - LEFT boundary correction, ONLY
C                2 - LEFT and RIGHT boundary correction

C              OUTPUT

C      fzz <-- hazard estimates at zz(1:m)
C      bopt <-- optimal local bandwidth at z(1:gridz)
C      bopt1 <-- bandwidth used to compute estimated hazards.
C                It is obtained from bopt, smoothing with the bandwidth bsmo
C      msemin <-- minimum MSE at each z, for local choice
C      biasmn <-- minimum bias at each z, for local choice
C      varmin  <-- minimum variance at each z, for local choice
C      imsemn <-- minimum IMSE,  for global and local choice
C      globlb  <-- optimal global bandwidth, resulting from MSE minimization
C      glmse <-- IMSE at bandwidth grid, for global choice
C      b <-- bandwidth grid for MSE minimization, for global choice
C
      INTEGER i, n, delta(n), gridb, gridz, ks, local, m, kflag
      DOUBLE PRECISION z(gridz), zz(m), fzz(m), bopt(gridz), x(n), endr,
     *     endl, bsmo, bopt1(m), bw(gridb), msemin(gridz), hpilot(1000),
     *     biasmn(gridz), varmin(gridz), bpilot, imsemn, globlb,
     *     glmse(gridb), hazden, w(n),naes(n), nw
      COMMON /hazpil/ hpilot

      EXTERNAL loclmn, glmin, bsmoth, hazden

C      call dblepr('naes', 4, naes(1), 1)
      
      IF (gridb .EQ. 1) THEN
         globlb = bw(1)
         GOTO 90
      END IF
      
     
C --- Compute the hazard estimates for bw=bpilot.
C     They will be used in msemse

      DO 10, i=1, gridz
         hpilot(i) = hazden(n, ks, x, delta, z(i), bpilot,
     *                   endl, endr, kflag, naes,nw)
 10   CONTINUE

      IF (local .EQ. 1) THEN

         CALL loclmn (n, x, delta, ks, z, gridz, bw, gridb,
     *         bopt, endl, endr, bpilot, msemin, biasmn, varmin,
     *         kflag, w, naes,nw)

C --- compute IMSE

         imsemn = 0.0d0
         DO 30, i = 1, gridz
            imsemn = imsemn + msemin(i)
 30      CONTINUE

C --- smooth the optimal local bandwidths

         CALL bsmoth(gridz, z, bopt, m, zz, bopt1, bsmo, ks, kflag,
     *                endl, endr)

      ELSE
    
         CALL glmin (n, x, delta, ks, z, gridz, bw, gridb,
     *         endl, endr, bpilot, imsemn, globlb, glmse, kflag,
     *         w, naes,nw)  
     
      END IF
      
 90   DO 40, i=1, m
         IF (gridb.EQ.1 .OR. local.EQ.0) THEN
            fzz(i) = hazden(n, ks, x, delta, zz(i), globlb,
     *                       endl, endr, kflag, naes,nw)
         ELSE
            fzz(i) = hazden(n, ks, x, delta, zz(i), bopt1(i),
     *                   endl, endr, kflag, naes,nw)
         END IF
 40   CONTINUE

      RETURN
      END

      SUBROUTINE olafbw(times, delta, n, z, gridz, k, bopt)

C Computes local bandwidths at each grid point in z

      INTEGER n, delta(n), gridz, k
      DOUBLE PRECISION times(n), z(gridz), bopt(gridz)

      INTEGER count
      DOUBLE PRECISION x(20000,2)

      EXTERNAL kapmei, locolf

C     Call Kaplan-Meier to compute the survival function

      CALL kapmei(times, delta, n, x, count)

C     now compute the bandwidths

      CALL locolf(x, n, z, gridz, count, k, bopt)

      RETURN
      END
      SUBROUTINE olafmn(x, delta, n, z, gridz, ks, endl, endr, bpilot,
     *                   bopt, kmin, kmax, kimse, kflag, w, naes, nw)

C     Computes the bandwidth at each grid point z
C     First it finds optimum number of neighbors, minimizing the IMSE

      INTEGER n, delta(n), gridz, ks, kmin, kmax, kflag
      DOUBLE PRECISION x(n), z(gridz), bopt(gridz), bpilot,
     *                 endl, endr, kimse(1)

      INTEGER k, i, kopt
      DOUBLE PRECISION imse, imsemn, bias, var, mse, bwi, zi, 
     *                 hpilot(1000), ZERO, HUGE, w(n),
     *                 naes(n), nw      
      PARAMETER (ZERO=0.0d0, HUGE=1.0d5)
      COMMON /hazpil/ hpilot

      EXTERNAL olafbw, msemse 

      IF (kmin .EQ. kmax) THEN
         CALL olafbw (x, delta, n, z, gridz, kmin, bopt)

         RETURN
      END IF

      imsemn = HUGE

C     For the moment, the maximum number of neighbors to be considered
C     is half of the non-censored observations

      DO 100, k=kmin, kmax

C     compute the bandwidths bopt for k neighbors

         CALL olafbw (x, delta, n, z, gridz, k, bopt)

         imse = ZERO

C     compute MSE at each gridpoint

         DO 50, i=1, gridz
            zi = z(i)
            bwi = bopt(i)

            CALL msemse(n, ks, zi, endl, endr, x, delta,
     *                  bwi, mse, bias, var, bpilot, hpilot(i), 
     *                  kflag, w, naes, nw)

            imse = imse + mse
 50      CONTINUE

         IF (imse .LT. imsemn) THEN
            kopt = k
            imsemn = imse
         END IF

         kimse(k-kmin+1) = imse
 100  CONTINUE

C     kopt is returned in kmin 

      kmin = kopt

C     compute the bandwidths for kopt

      CALL olafbw (x, delta, n, z, gridz, kopt, bopt)

      RETURN
      END
      DOUBLE PRECISION FUNCTION oneolf(x, n, x0, nobs, k)

C Computes bandwidth at grid point x0, using the survival matrix x

C     x --> survival matrix computed using Kaplan-Meier
C           first column: distinct survival times
C           second column: survival values
      DOUBLE PRECISION x(20000,2)

C     n --> unique survival times
      INTEGER n

C     x0 --> grid point for which a nearest neighbor bandwidth is to be
C            calculated
      DOUBLE PRECISION x0

C     nobs --> number of observations (censored and uncensored)
      INTEGER nobs

C     k --> number of neighbors within the bandwidth
      INTEGER k

      INTEGER i, ix0, ilo, ihi, count, atpos
      DOUBLE PRECISION dx(200000), const, bw, bw0, r, r0, gets, ds       
      DOUBLE PRECISION EPSI, ONE, ONEMOR, ONELES
      PARAMETER (EPSI=1.0D-5, ONE=1.0D0, ONEMOR=ONE+EPSI,
     &           ONELES=ONE-EPSI)
      
      EXTERNAL atpos, gets, sorter
      INTRINSIC min, max, abs

      ix0 = atpos(x(1,1), n, x0)

      ilo = max(1, ix0-k)
      ihi = min(n, ix0+k)

      count = 0

C     compute the distances from x0 to nearest neighbors

      DO 100, i=ilo, ihi
         count = count + 1
         dx(count) = abs(x(i,1)-x0)
 100  CONTINUE

      CALL sorter(dx, count)

C     compute the largest distance so that: S(x0-r)-S(x0+r-0) <= (k-1)/n

      const = ONEMOR * (k-1) / nobs

      bw = -99.99d0

      DO 200, i=1, count
         r = dx(i)

         ds = gets(x, n, x0-r) - gets(x, n, x0+r)

         IF (ds .GT. const) THEN
            GO TO 300
         ELSE
            bw = r
         END IF

 200  CONTINUE

 300  CONTINUE

C     Decide from bw, bw+, r-

      bw0 = ONEMOR * bw

      ds = gets(x, n, x0-bw0) - gets(x, n, x0+bw0)

      IF (ds .GT. const) THEN
         oneolf = bw
         RETURN
      END IF

      r0 = ONELES * r

      ds = gets(x, n, x0-r0) - gets(x, n, x0+r0)

      IF (ds .GT. const) THEN
         oneolf = bw0
      ELSE
         oneolf = r0
      END IF

      RETURN
      END
      SUBROUTINE sorter(v, n)

C     Sort the vector v, increasingly

C     length of v
      INTEGER n

C     vector to be sorted
      DOUBLE PRECISION v(n)

      LOGICAL qdone
      INTEGER i
      DOUBLE PRECISION temp

      IF (n .EQ. 1) RETURN

 100  qdone = .TRUE.

      DO 200, i=1, n-1

         IF (v(i) .GT. v(i+1)) THEN
            qdone = .FALSE.
            temp = v(i)
            v(i) = v(i+1)
            v(i+1) = temp
         END IF

 200  CONTINUE

      IF (.NOT. qdone) GO TO 100

      RETURN
      END
       DOUBLE PRECISION FUNCTION surfct(x, delta, n, xx, w, nw)

C Computes the empirical survival function of the UNCENSORED observations
C In the original code from Mueller, he uses ALL observations

       INTEGER n, delta(n)
       DOUBLE PRECISION x(n), xx, w(n), nw

       INTEGER index, i

       INTRINSIC DBLE

       index = 0
       
       DO 10, i=1, n

          IF ( (x(i).LE.xx) .AND. (delta(i).EQ.1) ) THEN
             index = index + DBLE(w(i))
          END IF

10      CONTINUE
C       nw=SUM(w)
        surfct = 1.0d0 - index/DBLE(nw+1)
        
        RETURN
        END
      SUBROUTINE try(n, ks, x, delta, z, b, endl, endr, q, r, s,
     *               valueb, valuev, iterat, bpilot, kflag, w, naes, nw)       

      INTEGER n, ks, delta(n), iterat, kflag
      DOUBLE PRECISION x(n), z, b, endl, endr, q, r, s, valueb, valuev,
     *                 bpilot, w(n), naes(n), nw

      INTEGER i, it
      DOUBLE PRECISION sumb, sumv, del, xx, br, bs, bxx, vr, vs, vxx,
     *                 tnm, ZERO, HALF
      PARAMETER (ZERO=0.0d0, HALF=0.5d0)
      INTRINSIC DBLE

      IF (iterat .EQ. 1) THEN

         CALL func(n, ks, x, delta, z, b, endl, endr,
     *             q, r, br, vr, bpilot, kflag, w, naes, nw)

         CALL func(n, ks, x, delta, z, b, endl, endr,
     *             q, s, bs, vs, bpilot, kflag, w, naes, nw)

         valueb = HALF*(s-r)*(br+bs)
         valuev = HALF*(s-r)*(vr+vs)
      ELSE
         it = 2**(iterat-2)
         tnm = DBLE(it)
         del = (s-r) / tnm
         xx = r + HALF*del

         sumb = ZERO
         sumv = ZERO

         DO 10, i=1, it

            CALL func(n, ks, x, delta, z, b, endl, endr,
     *                q, xx, bxx, vxx, bpilot, kflag, w, naes, nw)

            sumb = sumb + bxx
            sumv = sumv + vxx

            xx = xx + del
10       CONTINUE

         valueb = HALF*(valueb+(s-r)*sumb/tnm)
         valuev = HALF*(valuev+(s-r)*sumv/tnm)
      END IF

      RETURN
      END
