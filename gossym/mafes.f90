!MS$IF DEFINED (WINDOWS)
! note that windows is defined as a preprocessor directive
! program with windowing capabilities 

	integer function WinMain(hInstance, hPrevInstance, &
   	                         lpszCmdLine, nCmdShow)

    !MS$ATTRIBUTES STDCALL, ALIAS:'_WinMain@16'::WinMain

    USE MSFLIB
	USE MSFWIN

	INTEGER, intent(IN) :: hInstance, hPrevInstance
	INTEGER, intent(IN) :: lpszCmdLine
	INTEGER, intent(IN) :: nCmdShow
	INTEGER retval
	LOGICAL retlog

! FIND GUICS WINDOW

  

	iWndMain = FindWindow(NULL,"GUICS"C)
	IF(iWndMain.EQ.0)THEN		 
	   retval = MessageBox(NULL,"Interface is needed!"C, &
	                        'Error in GOSSYM'C,IOR(MB_OK,MB_ICONSTOP))
       WinMain = -1
	   RETURN
	ENDIF

	CALL ModelMain()

! *** 1028 is	code for END-OF-WORK message
! *** INFO1 = 0 and INFO2 = 0 may be some data to hand to interface

	retlog = PostMessage(iWndMain,1028,0,0)
	WinMain = 0

	END

    SUBROUTINE ModelMain()

!MS$ELSE
      PROGRAM MAFES
!MS$ENDIF
    USE MSFLIB
    USE MSFWIN

!*****************************************************************************
!
!                            MAFES Version 2004
!     
!           *********************  WARNING  **********************
!Temperatures must be between 13 and 37oC--will not extrapolate beyond limits                     
!           *********************  WARNING  **********************
!      
!             This is the application version of MAFES for GUICS
!
!*****************************************************************************

 	use common_block

    CHARACTER XFILE*150,DUMBO*1


1020 FORMAT(/// '    THE FOLLOWING MESSAGES ARE FOR PROFILE ',A15)

	iDanoflg = 0
	ipdays = 0
	co2 = 0
	runmode = 'gossym'
	opsys = ' '
	errfle = 'runtime.err'
	binfle = 'runtime.g00'
	xfile = 'run.dat'
	insfle = ''
	fngfle = ''

    OPEN(1,FILE=xfile,STATUS='OLD',err=200)
    READ(1,"(a)",END=200) actwth
    READ(1,"(a)",END=200) furwth
    READ(1,"(a)",END=200) profile
    READ(1,"(a)",END=200) irrfrt
    READ(1,"(a)",END=200) solhyd
    READ(1,"(a)",END=200) solimpd
    READ(1,"(a)",END=200) intsol
    READ(1,"(a)",END=200) variety
    READ(1,"(a)",END=200) pgrhrb
    READ(1,"(a)",END=200) pmafil
    READ(1,"(a)",END=200) plotfle
    READ(1,"(a)",END=200) listfle
    READ(1,"(a)",END=200) sumryfle

    loc0 = index(plotfle,'.g01')
    plotfle2 = plotfle(1:loc0) // 'g02'
    OPEN(24,FILE='plot.tmp', status='unknown')
!	loc0 = (index(profile,'\', BACK=.true.)+1)
!	profle = profile(loc0:)
!	loc1 = index(actwth,'\', BACK=.true.)+1
!	actwea = actwth(loc1:)
!	loc2 = index(furwth,'\', BACK=.true.)+1
!	furwea = furwth(loc2:)
!	loc3 = index(irrfrt,'\', BACK=.true.)+1
!	irrfle = irrfrt(loc3:)
!	loc4 = index(solhyd,'\', BACK=.true.)+1
!	hydfle = solhyd(loc4:)
!	loc5 = index(intsol,'\', BACK=.true.)+1
!	intfle = intsol(loc5:)
!	loc6 = index(pgrhrb,'\', BACK=.true.)+1
!	pgrfle = pgrhrb(loc6:)
!	loc7 = index(pmafil,'\', BACK=.true.)+1
!	pmafle = pmafil(loc7:)

	profle = profile((index(profile,'\', BACK=.true.)+1):)
	actwea = actwth((index(actwth,'\', BACK=.true.)+1):)
	furwea = furwth((index(furwth,'\', BACK=.true.)+1):)
	irrfle = irrfrt((index(irrfrt,'\', BACK=.true.)+1):)
	hydfle = solhyd((index(solhyd,'\', BACK=.true.)+1):)
	intfle = intsol((index(intsol,'\', BACK=.true.)+1):)
	pgrfle = pgrhrb((index(pgrhrb,'\', BACK=.true.)+1):)
	pmafle = pmafil((index(pmafil,'\', BACK=.true.)+1):)

	OPEN(21,FILE=binfle,FORM='UNFORMATTED',STATUS='unknown')
    OPEN(23,FILE=ERRFLE,STATUS='UNKNOWN')
    open(77,file=plotfle,status='unknown') 
    open(88,file=plotfle2,status='unknown') 
!MS$IF DEFINED (ALONE)
     Write(23,*) 'not a windows program'
!MS$ENDIF

      
140 CONTINUE
    READ(23,*,END=160,ERR=160) DUMBO
    GO TO 140
160 CONTINUE
    backspace(23)
    WRITE(23,1020) profle
    CALL GOSSYM
    CALL GOSLOOP
    CALL PrintOut
200 CONTINUE
	close(1)
	close(21,status='delete')
	close(22)
	close(23)
 	close(75)
 	close(77)
 	close(88)
    close(LSTNG)


!MS$IF DEFINED(WINDOWS)
	return
!MS$ENDIF
	end


      SUBROUTINE MEMPLANT

      use common_block

      PMAFIL='IAMDONE'
      ABEND=.TRUE.
      RETURN
      END


      SUBROUTINE PHLPLT(PROFLE)
      CHARACTER*13 PROFLE
      PROFLE=PROFLE
      RETURN
      END


      FUNCTION CKUSER(I)
      INTEGER*2 CKUSER,I
      JDUM=I
      CKUSER=0
      RETURN
      END


      SUBROUTINE GOSLOOP

      use common_block

  300 CONTINUE
        IF(SEND.OR.ABEND) THEN
          RETURN
        ELSE
          CALL DAILY
        ENDIF
      GO TO 300
      END


      SUBROUTINE GOSSYM

      use common_block

      CALL GBLOCK
      IF(.NOT.ABEND) CALL PROFILES
      IF(.NOT.ABEND) CALL WEATHER
      IF(.NOT.ABEND) CALL SOILHYDR
      IF(.NOT.ABEND) CALL SOILIMPD
      IF(.NOT.ABEND) CALL INITSOIL
      IF(.NOT.ABEND) CALL AGINPUTS
!      IF(.NOT.ABEND) CALL INSECTCDE
!      IF(.NOT.ABEND) CALL FUNGICDE
      IF(.NOT.ABEND) CALL PGRHRBCDE
      IF(.NOT.ABEND) CALL PLANTMAPS
      IF(.NOT.ABEND) CALL INITIALIZE
      RETURN
      END


      SUBROUTINE DAILY

      use common_block

      INTEGER*2 IDUM,CKUSER
 
      IDAY=IDAY+1
      CALL CLYMAT
      KDAY=DAYNUM-EMERGE+1
      IF(KDAY.LT.0) KDAY=0
      dae = kday
 	  if((iniFertDate.lt.emerge).and.(kday.eq.0)) then
         CALL SOIL
	  else
         CALL SOIL
         IF(DAYNUM.GE.EMERGE) THEN    
            IF((PIXDAY(1).GT.0).AND.(DAYNUM.GE.PIXDAY(1))) CALL PIX
            IF((DEFBGN.GT.0).AND.(DAYNUM.GE.DEFBGN)) call defoliat
            CALL PNET
            CALL GROWTH
            CALL PLTMAPS
            CALL ABCISE
            DO 320 I=1,15
               IPLTNO=I
               IF(DAYNUM.EQ.MSADTE(I)) CALL MDSEA(IPLTNO)
  320       CONTINUE
            CALL MATBAL
         ENDIF
	  endif
      CALL OUTPUT
      IDUM = 0
      IF(CKUSER(IDUM).GT.0) THEN
        ABEND=.TRUE.
      ELSEIF(DAYNUM.GE.JDSTPS) THEN
        SEND=.TRUE.
      ELSEIF(DAYNUM.GE.JDAYLW) THEN
        SEND=.TRUE.
      ENDIF
      RETURN
      END


      SUBROUTINE PROFILES

      use common_block

      INTEGER*2 ISKWTH,IROWSP,NACRES
      CHARACTER EMERG*10,STRSIM*10,STPSIM*10,PRFNAM*13

 1000 FORMAT('1',28X,'MAFES-GOSSYM (TM)', /	&
             29X,'COPYRIGHT (C)  2003' / &
             12X,'Mississippi Agriculture & Forestry Experiment Station' /	&
             16X,'Contact:  Dr. K.R. Reddy and David W. Brand' /	&
             32X,'(662) 325-9462'/	&
             29X,'Version - ',A8 /)
 1020 FORMAT(A80)
 1040 FORMAT(// 29X,'Profile Name: ',A13 / 7X,'Description: ',A51)
 1060 FORMAT(/	&
          4X,'Variety............ ',A17,2X,'Latitude........... ',I11 / &
          4X,'Start Simulation... ',I3,'/',I2,'/',I4,	&
          8X,'Stop Simulation.... ',I3,'/',I2,'/',I4/	&
          4X,'Emergence.......... ',I3,'/',I2,'/',I4,	&
          8X,'Season Length...... ',I11 /				&
          4X,'Row Spacing (in)... ',f11.2,8X,'Plants Per Row-Ft.. ',F11.2/ &
          4X,'Skip Width (in).... ',f11.2,8X,'Plants Per Acre.... ',I11)	
 1080 FORMAT(/	&
          4X,'Actual Weather..... ',A15,6X,'Future Weather..... ',A15 / &
          4X,'Predicted Weather.. ',A15,6X,'Cultural Inputs.... ',A15 / &
          4X,'Soil Hydrology..... ',A15,6X,'Initial Soil....... ',A15 / &
          4X,'PGR & Herbicides... ',A15,6X,'Insecticides....... ',A15 / &
          4X,'Fungicides......... ',A15,6X,'Plant Maps......... ',A15/)

      OPEN(11,FILE=profile,FORM=FILFRM,STATUS='OLD',ERR=570)
      READ(11,*,ERR=580,END=600) PRFNAM,PDESCP,rundate,EMERG,STRSIM,STPSIM, &
                                 LATUDE,PLTPFT,ISKWTH,IROWSP,NACRES, &
					             iDanoflg,(ioutfg(i),i=1,6), CO2
! Added CO2 as an input from profile, 6/1/06, SK
! See "initialize" for CO2 correction routines

      CLOSE(11)
	  ivarty = 0
      LSTNG=7
      OPEN(LSTNG,FILE=listfle,STATUS='UNKNOWN')
      WRITE(LSTNG,1000) VERSION
      WRITE(LSTNG,1040) profle,PDESCP
!cdt      
      OPEN(99,file='debug.dat')

      CALL CALTOJULAN(EMERG,IDUM1,IDUM2,IYEAR,EMERGE)
      CALL CALTOJULAN(STRSIM,IDUM1,IDUM2,IDUM3,JDSTRS)
      CALL CALTOJULAN(STPSIM,IDUM1,IDUM2,IDUM3,JDSTPS)

! *** Routine to vary the start and end of simulation   kit  6/20/01

   	  if(ipdays.ne.0) then

! ***	reset emergence, start and stop dates

	     emerge = emerge + ipdays
!	     jdstrs = emerge - 1 
	     jdstrs = jdstrs + ipdays
	     jdstps = jdstps + ipdays
 	  endif

! ***	end of planting date routine

      IF((JDSTRS.LE.0).OR.(JDSTPS.LE.0)) THEN
         PRINTBUF=' error on reading date in profile'
         CALL WRITERR
         ABEND=.TRUE.
         RETURN
      ELSE
         IF((JDSTRS.GT.EMERGE).AND.(EMERGE.GT.0)) THEN
           JDSTRS=EMERGE
           STRSIM=EMERG
         ELSEIF(EMERGE.LE.0) THEN
           EMERGE = 366
         ENDIF
         IVARTY = IVARTY+1
         SKIPWD = ISKWTH*2.54
         ROWSP = IROWSP*2.54
         ROWSAVE = ROWSP
         IF(SKIPWD.GT.1.0) ROWSP = (ROWSP+SKIPWD)/2.0
         POPPLT = PLTPFT/ROWSP*1327709.
         IF(.NOT.ABEND) CALL VARIETYS
         J0=POPPLT+.5
         J1=JDSTPS-JDSTRS+1

! ***	get the emergence,start and stop dates

         CALL JULANTOCAL(MO,DAZE,IYEAR,emerge)
         IMoemg = MO
         IDyemg= DAZE
	     iyremg= iyear
         CALL JULANTOCAL(MO,DAZE,IYEAR,jdstrs)
         IMostr = MO
         IDystr = DAZE
         iyrstr = iyear
         CALL JULANTOCAL(MO,DAZE,IYEAR,jdstps)
         IMostp = MO
         IDystp = DAZE
         iyrstp = iyear

         WRITE(LSTNG,1060) VARITY(IVARTY),LATUDE,imostr,idystr,iyrstr, &
	                       imostp,idystp,iyrstp,imoemg,idyemg,iyremg, &
                           J1,ROWSP/2.54,PLTPFT,SKIPWD/2.54,J0
        WRITE(LSTNG,1080) actwea,furwea,furwea,irrfle,hydfle,intfle, &
                          pgrfle,insfle,fngfle,pmafle
      ENDIF
      RETURN
  570 CONTINUE
      PRINTBUF=' error on open of profile file'
      CALL WRITERR
      ABEND=.TRUE.
      RETURN
  580 CONTINUE
      PRINTBUF=' error on read in profile'
      CALL WRITERR
      ABEND=.TRUE.
      CLOSE(11)
      RETURN
  600 CONTINUE
      PRINTBUF=' end-of-file while reading profile'
      CALL WRITERR
      ABEND=.TRUE.
      CLOSE(11)
      RETURN
      END


      SUBROUTINE WEATHER
!                      W E A T H E R
!  ************************************************************
!  *   SUBROUTINE TO READ IN WEATHER, FERTILIZATION, AND      *
!  *       IRRIGATION                                         *
!  ************************************************************
 
      use common_block
 
      DIMENSION C13(7)
      CHARACTER A00*5,A01*1,A02*6

 1000 FORMAT(4X,'Last Actual Weather ',I2.2,'/',I2.2,'/',I4.4, &
            13X,'Last Future Weather ',I2.2,'/',I2.2)
 1020 FORMAT(A1)
 1040 FORMAT(' Error, missing data in actual weather on day ',A5)
 1060 FORMAT(' Error, missing data in final weather on day ',I5)
 1080 FORMAT(' Error, final weather is missing day ',I5)
 1100 FORMAT(A80)
 1120 FORMAT(' Error Opening ',A6,' Weather File.')
 1140 FORMAT(' Error Reading ',A6,' Weather File on day ',A5)

      OPEN(23,FILE=ERRFLE,STATUS='UNKNOWN')
  100 CONTINUE
      READ(23,*,END=120,ERR=100) A01
      GO TO 100
  120 CONTINUE
      BACKSPACE(23)

! *** READ FUTURE WEATHER FILE

      OPEN(11,FILE=furwth,FORM='FORMATTED',STATUS='OLD',ERR=220)

      C13(4)=0
      C13(7)=1
  200 CONTINUE
      READ(11,1020,END=240,ERR=230) A01
      IF(A01.EQ.'#') THEN
         GO TO 200
      ELSE
         BACKSPACE(11)
         READ(11,*,END=240,ERR=230) C13(7),A00,(C13(I),I=1,3), &
                                    C13(5),C13(6)
         J=C13(7)
         IF(C13(6).LE.0.) C13(6)=88.
         IF(J.GT.0) THEN
            DO 210 IC=1,7
               CLIMAT(J,IC) = C13(IC)
  210       CONTINUE
         ENDIF
      ENDIF
      GO TO 200
  220 CONTINUE
      A02 = 'FUTURE'
      WRITE(PRINTBUF,1120) A02
      WRITE(23,1100) PRINTBUF
      ABEND=.TRUE.
      CLOSE(23)
      RETURN
  230 CONTINUE
      A02 = 'FUTURE'
      WRITE(PRINTBUF,1140) A02,A00
      WRITE(23,1100) PRINTBUF
  240 CONTINUE
      LDAYFW = C13(7)
      C13(7) = 1
      JDAYLW = LDAYFW
      CLOSE(11)

! *** READ IN CURRENT WEATHER FILE ON TOP OF FUTURE WEATHER

      OPEN(11,FILE=actwth,FORM='FORMATTED',STATUS='OLD',ERR=330)

	  iFirstday = 0

  300 CONTINUE
      READ(11,1020,END=360,ERR=340) A01
      IF(A01.EQ.'#') THEN
         GO TO 300
      ELSE
         BACKSPACE(11)
         READ(11,*,END=360,ERR=340) C13(7),A00,(C13(I),I=1,3), &
                                    C13(5),C13(6)

! *** variable to store the date of first available actual weather data

	 	 iFirstday = iFirstday + 1
		 if(iFirstday.eq.1) iFDayAW = C13(7)
         IF((C13(7).GE.iFDayAW).AND.(C13(7).LE.JDSTPS)) THEN

!          IF((C13(7).GE.JDSTRS).AND.(C13(7).LE.JDSTPS)) THEN
! *** kit      06/24/2003

             DO 310 I=1,6
                IF(ABS(ABS(C13(I))-6999.0).LT..01) THEN
                   WRITE(PRINTBUF,1040) A00
                   WRITE(23,1100) PRINTBUF
                   GO TO 300
                ELSEIF(ABS(ABS(C13(I))-9111.0).LT..01) THEN
                   WRITE(PRINTBUF,1040) A00
                   WRITE(23,1100) PRINTBUF
                   GO TO 300
                ENDIF
  310        CONTINUE

!            J=C13(7)-JDSTRS+1

             J=C13(7)
             IF(C13(6).LE.0.) C13(6)=88.
             DO 320 IC=1,7
                CLIMAT(J,IC) = C13(IC)
  320        CONTINUE
         ENDIF
      ENDIF
      GO TO 300
  330 CONTINUE
      A02 = 'ACTUAL'
      WRITE(PRINTBUF,1120) A02
      WRITE(23,1100) PRINTBUF
      ABEND=.TRUE.
      CLOSE(23)
      RETURN
  340 CONTINUE
      A02 = 'ACTUAL'
      WRITE(PRINTBUF,1140) A02,A00
      WRITE(23,1100) PRINTBUF
  360 CONTINUE
      CLOSE(11)
      LDAYAW = C13(7)
      IF(LDAYAW.GT.JDAYLW) JDAYLW=LDAYAW
      CALL JULANTOCAL(I0,I1,IYEAR,LDAYAW)
      CALL JULANTOCAL(I2,I3,IYEAR,LDAYFW)
      WRITE(lstng,1000) I0,I1,IYEAR,I2,I3

      DO 400 K=JDSTRS,JDSTPS
        DO 380 J=1,6
          IF((ABS(CLIMAT(K,J))-6999..LT..01).AND. &
             (ABS(CLIMAT(K,J))-6999..GT.-.01)) THEN
            WRITE(PRINTBUF,1060) CLIMAT(K,7)
            WRITE(23,1100) PRINTBUF
            ABEND=.TRUE.
          ELSEIF((ABS(CLIMAT(K,J))-9111..LT..01).AND. &
                 (ABS(CLIMAT(K,J))-9111..GT.-.01)) THEN
            WRITE(PRINTBUF,1060) CLIMAT(K,7)
            WRITE(23,1100) PRINTBUF
            ABEND=.TRUE.
          ENDIF
  380   CONTINUE
        IF(CLIMAT(K,7).LT.0.01) THEN
          WRITE(PRINTBUF,1080) K
          WRITE(23,1100) PRINTBUF
          ABEND=.TRUE.
        ENDIF
  400 CONTINUE
      CLOSE(23)
      RETURN
      END


      SUBROUTINE SOILHYDR

      use common_block

      INTEGER*4 POINTR
      CHARACTER FILNAM*13,DSCRIP*51,XDUM1*3,XDUM2*3

 1220 FORMAT(//	 &
             7X,'NUMBER OF SOIL LAYERS ... ',I1, &
            10X,'WATER TABLE DEPTH (IN) .. ',F4.0, //  &
            16X,'LAYER NO.    LAYER DEPTH    % SAND    % CLAY' / &
            16X,'                 (IN)')
 1240 FORMAT(18X,I3,10X,F6.0,7X,I4,6X,I4)
 1260 FORMAT(/ 8X, 'Rainfall Runoff: ',A3,18X,'Irrigation Runoff: ',A3)

      OPEN(11,FILE=solhyd,FORM=FILFRM,STATUS='OLD',ERR=200)
      READ(11,*,ERR=200,END=300) FILNAM,DSCRIP,LYRSOL,DIFF0C,THTA0C, &
                BETAC,THTASC,FCINIC,THTARC,AIRDRC,BDC,DIFF0W,THTA0W, &
                BETAW,THTASW,FCININW,THTARW,AIRDRW,BDW,TD,THETAI, &
                BDSLOP,BDRATO,PSISFC,WATTBL,POINTR
      DO 120 I=1,LYRSOL
         READ(11,*,ERR=200,END=300) LDEPTH(I),DIFF0(I),THETA0(I), &
                BETA(I),THETAS(I),FCININ(I),THETAR(I),AIRDR(I),BD(I), &
                IPSAND(I),IPCLAY(I),POINTR
  120 CONTINUE
      READ(11,*,ERR=150,END=150) IFGRAIN,IFGIRR
  150 CONTINUE
      CLOSE(11)
      IF(WATTBL.LT.1) WATTBL=200
      IF(PSISFC.GT.0.) PSISFC=-PSISFC
      WRITE(lstng,1220) LYRSOL,WATTBL/2.54
      DO 170 I=1,LYRSOL
          WRITE(lstng,1240) I,LDEPTH(I)/2.54,IPSAND(I),IPCLAY(I)
  170 CONTINUE
      XDUM1='ON '
      XDUM2='ON '
      IF(IFGRAIN.GT.0) XDUM1='OFF'
      IF(IFGIRR.GT.0) XDUM2='OFF'
      WRITE(lstng,1260) XDUM1,XDUM2
      RETURN
  200 CONTINUE
      PRINTBUF=' error on read in soilhydr'
      CALL WRITERR
      CLOSE(11)
      ABEND = .TRUE.
      RETURN
  300 CONTINUE
      PRINTBUF=' end-of-file while reading soilhydr'
      CALL WRITERR
      CLOSE(11)
      ABEND = .TRUE.
      RETURN
      END


      SUBROUTINE SOILIMPD

      use common_block

      OPEN(11,FILE=solimpd,FORM=FILFRM,STATUS='OLD',ERR=200)
      READ(11,*,ERR=200,END=300) SNAME
      READ(11,*,ERR=200,END=300) NCURVE
      DO 180 I=1,NCURVE
        READ(11,*,ERR=200,END=300) INRIM,GH2OC(I)
        DO 160 J=1,INRIM
          READ(11,*,ERR=200,END=300) TSTBD(I,J),TSTIMP(I,J)
  160   CONTINUE
  180 CONTINUE
      CLOSE(11)
      RETURN
  200 CONTINUE
      PRINTBUF=' error on read in soilimpd'
      CALL WRITERR
      ABEND = .TRUE.
      CLOSE(11)
      RETURN
  300 CONTINUE
      PRINTBUF=' end-of-file while reading soilimpd'
      CALL WRITERR
      ABEND = .TRUE.
      CLOSE(11)
      RETURN
      END


      SUBROUTINE INITSOIL

      use common_block

      INTEGER*2 ISOILW(12)
      CHARACTER FILNAM*13,DSCRIP*51

 1200 FORMAT(//	&
      '  DEPTH     RESIDUAL N (LB/ACRE)     ORGANIC     H2O CONTENT' / &
      '  (IN)      NITRATE    AMMONIA     MATTER (%)   % FIELD CAP'/ &
      '   1- 6',2(4X,F7.1),8X,F6.2,8X,F7.0 / &
      '   7-12',2(4X,F7.1),8X,F6.2,8X,F7.0 / &
      '  13-18',2(4X,F7.1),8X,F6.2,8X,F7.0 / &
      '  19-24',2(4X,F7.1),8X,F6.2,8X,F7.0 / &
      '  25-30',2(4X,F7.1),8X,F6.2,8X,F7.0 / &
      '  31-36',2(4X,F7.1),8X,F6.2,8X,F7.0)

      OPEN(11,FILE=intsol,FORM=FILFRM,STATUS='OLD',ERR=200)
      DO 60 I=1,12
        ISOILW(I)=100
   60 CONTINUE
      READ(11,*,ERR=200,END=80) FILNAM,DSCRIP,(RNNH4(I),RNNO3(I), &
                OMA(I),ISOILW(I),I=1,12)
   80 CONTINUE
      DO 100 I=1,12
        H2OINT(I) = ISOILW(I)
  100 CONTINUE
      WRITE(lstng,1200) (RNNO3(I),RNNH4(I),OMA(I),H2OINT(I),I=1,6)
      CLOSE(11)
      RETURN
  200 CONTINUE
      PRINTBUF=' error on read in initsoil'
      CALL WRITERR
      ABEND = .TRUE.
      CLOSE(11)
      RETURN
      END


      SUBROUTINE AGINPUTS

      use common_block

      REAL CSTIRR,AMTAMM,AMTNIT,AMTURA,CSTFRT
      INTEGER*2 MTHFRT,ISDHRZ,ISDDPH
      CHARACTER FILNAM*13,DSCRIP*51,CDATE*10,DBDATE*8
 
 1000 FORMAT( 'TABLE:', / &
      '          ------IRRIGATION------  ----------------FERTILIZER', &
      '-----------------' / &
      '   DATE   AMOUNT  METHOD   COST    NH4   NO3   UREA   COST ', &
      '  METHOD  LOCATION' /  &
      '           (in)            ($/a)  ----(lb N / a)----  ($/a)', &
      '         dpth side' /)
 1020 FORMAT(1X,A10,1X,F5.2,1X,A6,1X,F7.2,F7.1,2F6.1,2X,F6.2,2X,A6,2I4)
 1040 FORMAT( //  &
      '          ------DRIP IRRIGATION------  ----------------',  &
      'DRIP FERTILIZER---------------' /  &
      '  START      STOP     AMOUNT   COST   NH4   NO3   UREA   COST  ', &
      '  METHOD  LOCATION' / &
      '   DATE      DATE              ($/a)  ----(lb N / a)---  ($/a)',	&
      '          side dpth' /)
 1060 FORMAT(1X,A10,1X,A10,2X,F5.2,1X,F7.2,3F7.1,2X,F6.2,2X,A6,2I5)
 1080 FORMAT(I2.2,'/',I2.2,'/',I4.4)
 1800 FORMAT(' ')

      IF(IRRFRT.NE.'             ') THEN
         OPEN(11,FILE=irrfrt,FORM=FILFRM,STATUS='OLD',ERR=300)
         IDRIP = 0
         IDUM = 0
         READ(11,*,END=300,ERR=300) FILNAM,DSCRIP
         IF(IDUM.EQ.0) THEN
            WRITE(lstng,1800)
            WRITE(lstng,1800)
	 	    WRITE(lstng,1000)
            IDUM=1
         ENDIF
  200    CONTINUE
		 id = 1
         READ(11,*,END=300,ERR=300) CDATE,H2OAMT,MTHIRR(ID), &
                                    AMTAMM,AMTNIT,AMTURA, &
                                    MTHFRT,ISDHRZ,ISDDPH
         CALL CALTOJULAN(CDATE,IDUM1,IDUM2,IDUM3,I)

! *** Check to see if fertilizer application is before start
! *** of simulation.  If desired, set fertilizer application
! *** date to start date of simulation.      Kit 9/1/99

!	    if(i.lt.jdstrs) i = jdstrs

		if(i.le.iFDayAW) i = iFDayAW
	    if(i.le.jdstrs) then
		   jdstrs = i
		   iniFertDate = jdstrs
		endif

! ***	Reset fertilizer application date when planting dates are changed

		if(ipdays.ne.0) i = jdstrs

        IF(I.LT.1) I=1
		MTHIRR(I) = MTHIRR(ID)

          IF((MTHIRR(I).EQ.2).OR.(MTHIRR(I).EQ.4)) THEN
            IF(H2OAMT.GT.0.0) THEN
              IF(IDRIP.LT.1) IDRPBDT = I
              MTHDRP = MTHIRR(I)
              DRPH2O = H2OAMT
              IDRPFRT = MTHFRT
              DRPAMM = AMTAMM
              DRPNIT = AMTNIT
              DRPURA = AMTURA
              IDRPHRZ = ISDHRZ
              IDRPDPH = ISDDPH
              IDRIP = 1
            ELSE
              WRITE(lstng,1040)
              CALL JULANTOCAL(IDUM1,IDUM2,IYEAR,IDRPBDT)
              WRITE(DBDATE,1080) IDUM1,IDUM2,IYEAR
              WRITE(lstng,1060) DBDATE,CDATE,DRPH2O,0.0,DRPAMM,DRPNIT, &
                                DRPURA,0.0,FMTHOD(2),IDRPHRZ,IDRPDPH

              IF((I.GT.IDRPBDT).AND.(IDRIP.GT.0)) THEN
                IF(IDRPFRT.NE.4) IDRPFRT = 1
                DO J=IDRPBDT,I-1
                  AMTIRR(J) = AMTIRR(J)+DRPH2O
                  MTHIRR(J) = MTHDRP
                  NFERT(J,1) = 0
                  NFERT(J,2) = NFERT(J,2)+DRPAMM+.5
                  NFERT(J,3) = NFERT(J,3)+DRPNIT+.5
                  NFERT(J,4) = NFERT(J,4)+DRPURA+.5
                  NFERT(J,5) = IDRPFRT
                  NFERT(J,6) = IDRPHRZ
                  NFERT(J,7) = IDRPDPH
                ENDDO
              ENDIF
              IDRIP = 0
            ENDIF
          ELSE
            IDUM1=6
            IF(H2OAMT.GT.0.) THEN
              IDUM1=MTHIRR(I)+1
              AMTIRR(I)=AMTIRR(I)+H2OAMT
            ENDIF
            IDUM2=6
            IF((AMTAMM.GT..01).OR.(AMTNIT.GT..01).OR. &
                                  (AMTURA.GT..01)) THEN
              IDUM2=MTHFRT+1
              NFERT(I,1) = 0
              NFERT(I,2) = NFERT(I,2)+AMTAMM+.5
              NFERT(I,3) = NFERT(I,3)+AMTNIT+.5
              NFERT(I,4) = NFERT(I,4)+AMTURA+.5
              NFERT(I,5) = MTHFRT
              NFERT(I,6) = ISDHRZ
              NFERT(I,7) = ISDDPH
            ENDIF
			CSTIRR = 0.00
			CSTFRT = 0.00
            WRITE(lstng,1020) CDATE,AMTIRR(I),IMTHOD(IDUM1), &
               CSTIRR,AMTAMM,AMTNIT,AMTURA,CSTFRT,FMTHOD(IDUM2), &
               ISDDPH,ISDHRZ

          ENDIF
        GO TO 200
      ENDIF
  300 CONTINUE
      CLOSE(11)
  650 CONTINUE
      RETURN
      END


      SUBROUTINE INSECTCDE

      use common_block

      INTEGER*2 INSMTH(10),INSBDW(10),INUNIT(10)
      INTEGER*4 IPOINT
      REAL CSTINS(10),RTEINS(10)
      CHARACTER INSBRD(10)*14,INSNAM*13,DESCRP*51,INSDTE*10

      IF(INSCDE.NE.'             ') THEN
        OPEN(11,FILE=INSCDE,FORM=FILFRM,STATUS='OLD',ERR=320)
        I=0
  100   CONTINUE
        I=I+1
        READ(11,*,ERR=200,END=300) INSNAM,DESCRP,INSDTE
        DO 160 J=1,10
           READ(11,*,ERR=200,END=200) INSBRD(J),INSMTH(J), &
                     INSBDW(J),CSTINS(J),RTEINS(J),INUNIT(J)
  160   CONTINUE
        READ(11,*,ERR=200,END=200) IPOINT
        GO TO 100
  200   CONTINUE
          PRINTBUF=' error on read in insecticide file'
          CALL WRITERR
  300   CONTINUE
        CLOSE(11)
      ENDIF
  320 CONTINUE        
      RETURN
      END


      SUBROUTINE FUNGICDE

      use common_block

      INTEGER*2 FNGMTH(10),FNGBDW(10),FGUNIT(10)
      INTEGER*4 IPOINT
      REAL CSTFNG(10),RTEFNG(10)
      CHARACTER FNGBRD(10)*14,FNGNAM*13,DESCRP*51,FNGDTE*10

      IF(FNGCDE.NE.'             ') THEN
        OPEN(11,FILE=FNGCDE,FORM=FILFRM,STATUS='OLD',ERR=320)
        I=0
  100   CONTINUE
        I=I+1
        READ(11,*,ERR=200,END=300) FNGNAM,DESCRP,FNGDTE
        DO 160 J=1,10
           READ(11,*,ERR=200,END=200) FNGBRD(J),FNGMTH(J), &
                     FNGBDW(J),CSTFNG(J),RTEFNG(J),FGUNIT(J)
  160   CONTINUE
        READ(11,*,ERR=200,END=200) IPOINT
        GO TO 100
  200   CONTINUE
          PRINTBUF=' error on read in fungicide file'
          CALL WRITERR
  300   CONTINUE
        CLOSE(11)
      ENDIF
  320 CONTINUE
      RETURN
      END


      SUBROUTINE PGRHRBCDE

      use common_block

      INTEGER*2 PGRMTH,PGRBDW,PGUNIT
      REAL RTEPGR
      CHARACTER PGRBRD*14,PGRNAM*13,DESCRP*51,PGRDTE*10

 1000 FORMAT(///  &
      '   DATE      BRAND       METHOD     BAND WTH      COST  ', &
      '   RATE      UNITS' /  &
      '             NAME                     (in)       ($/a)' /)
 1020 FORMAT(1X,A10,2X,A8,4X,A6,5X,I5,6X,F7.2,2X,F7.2,5X,A6)
 1040 FORMAT(13X,     A8,4X,A6,5X,I5,6X,F7.2,2X,F7.2,5X,A6)

      IF(PGRHRB.NE.'             ') THEN
        OPEN(11,FILE=PGRHRB,FORM=FILFRM,STATUS='OLD',ERR=380)

        IPX = 0
        IPRP = 0
        IDEF = 0
        I = 0
          IF((OPSYS.EQ.'dos').OR.(OPSYS.EQ.'DOS')) THEN
            READ(11,ERR=200,END=300) PGRNAM,DESCRP
          ELSE
            READ(11,*,ERR=200,END=300) PGRNAM,DESCRP
          ENDIF
          I=I+1
  100   CONTINUE
		  j = 1
            IF((OPSYS.EQ.'dos').OR.(OPSYS.EQ.'DOS')) THEN
              READ(11,ERR=200,END=300) PGRDTE,PGRBRD,PGRMTH, &
                                       PGRBDW,RTEPGR,PGUNIT
            ELSE
              READ(11,*,ERR=200,END=300) PGRDTE,PGRBRD,PGRMTH,	&
                                         PGRBDW,RTEPGR,PGUNIT
            ENDIF

            CALL CALTOJULAN(PGRDTE,IDUM1,IDUM2,IDUM3,jdpgr)
	  	    if((jdpgr.ge.JDSTRS).and.(jdpgr.le.JDSTPS))then
              IF(I.EQ.1) THEN
                 WRITE(lstng,1000)
			      i = 0
              ENDIF

            CSTPGR = 0.0
            L=2
  120       CONTINUE
              IF((PGRBRD(L:L).NE.'.').AND.(PGRBRD(L:L).NE.' ') &
                  .AND.(PGRBRD(L:L).NE.CHAR(0)).AND.(L.LT.14)) THEN
                L=L+1
                GO TO 120
              ENDIF
            DO 140 K=L,14
              PGRBRD(K:K)=' '
  140       CONTINUE
            IF(((PGRBRD(1:3).EQ.'pix').OR.(PGRBRD(1:3).EQ.'PIX')) &
                                      .AND.(RTEPGR.GT.0.01)) THEN
              IPX=IPX+1
              CALL CALTOJULAN(PGRDTE,IDUM1,IDUM2,IDUM3,PIXDAY(IPX))
              DUMRTE=RTEPGR
              IF(PGUNIT.EQ.1) DUMRTE=RTEPGR*8
              IF(PGUNIT.EQ.2) DUMRTE=RTEPGR/16.
              IF(PGUNIT.EQ.4) DUMRTE=1./RTEPGR
              IF(PGUNIT.EQ.5) DUMRTE=8./RTEPGR
              PIXPPA(IPX)=DUMRTE
              PIXMTH(IPX)=PGRMTH
            ELSEIF(((PGRBRD(1:4).EQ.'prep').OR.	&
                   (PGRBRD(1:4).EQ.'PREP')).AND.(RTEPGR.GT.0.01)) THEN
              IPRP = IPRP + 1
              CALL CALTOJULAN(PGRDTE,IDUM1,IDUM2,IDUM3,PRPDATE(IPRP))
              DUMRTE=RTEPGR
              IF(PGUNIT.EQ.1) DUMRTE=RTEPGR*8
              IF(PGUNIT.EQ.2) DUMRTE=RTEPGR/16.
              IF(PGUNIT.EQ.4) DUMRTE=1./RTEPGR
              IF(PGUNIT.EQ.5) DUMRTE=8./RTEPGR
              PRPPPA(IPRP) = DUMRTE
              PRPMTH(IPRP) = PGRMTH
            ELSEIF(((PGRBRD(1:3).EQ.'DEF').OR. &
                    (PGRBRD(1:5).EQ.'DROPP').OR. &
                    (PGRBRD(1:7).EQ.'HARVADE').OR. &
                    (PGRBRD(1:8).EQ.'GRAMOXON').OR.	&
                    (PGRBRD(1:3).EQ.'def').OR. &
                    (PGRBRD(1:5).EQ.'dropp').OR. &
                    (PGRBRD(1:7).EQ.'harvade').OR. &
                    (PGRBRD(1:8).EQ.'gramoxon')).AND. &
                                            (RTEPGR.GT.0.01)) THEN
              IDEF = IDEF + 1
              CALL CALTOJULAN(PGRDTE,IDUM1,IDUM2,IDUM3,DEFDATE(IDEF))
              DUMRTE=RTEPGR
              IF(PGUNIT.EQ.1) DUMRTE=RTEPGR*8
              IF(PGUNIT.EQ.2) DUMRTE=RTEPGR/16.
              IF(PGUNIT.EQ.4) DUMRTE=1./RTEPGR
              IF(PGUNIT.EQ.5) DUMRTE=8./RTEPGR
              DEFPPA(IDEF)=DUMRTE
              DEFMTH(IDEF)=PGRMTH
            ENDIF

            IF(PGRBRD(1:1).NE.CHAR(0)) THEN
              K=PGRMTH+1
              L=PGUNIT+1
              IF(J.EQ.1) THEN
                WRITE(lstng,1020) PGRDTE,PGRBRD(1:8),HMTHOD(K), &
                                  PGRBDW,CSTPGR,RTEPGR,PGRUNT(L)
              ELSE
                WRITE(lstng,1040) PGRBRD(1:8),HMTHOD(K), &
                                  PGRBDW,CSTPGR,RTEPGR,PGRUNT(L)
              ENDIF
            ENDIF
		endif

        GO TO 100
  200   CONTINUE
        PRINTBUF=' error on read in PGR & Herbicide file'
        CALL WRITERR
        CLOSE(11)
        RETURN
  300   CONTINUE
        CLOSE(11)
        IF(IDEF.GT.0) THEN
          IF((DEFDATE(1).LT.DEFBGN).OR.(DEFBGN.LE.0)) &
                                              DEFBGN = DEFDATE(1)
          DEFDAY = DEFDATE(1)
        ENDIF
        IF(IPRP.GT.0) THEN
          IF((PRPDATE(1).LT.DEFBGN).OR.(DEFBGN.LE.0)) &
                                              DEFBGN = PRPDATE(1)
          PRPDAY = PRPDATE(1)
        ENDIF
      ENDIF
  380 CONTINUE
      RETURN
      END


      SUBROUTINE PLANTMAPS

! *** Original PLANTMAPS subroutine with DNBaker's code ***
! *** MAX # Of: VBranch = 3, FBranch = 30, Fruiting Sites = 5

      use common_block

      INTEGER*2 IPLTHT
      CHARACTER FILNAM*13,DSCRIP*51
      INTEGER*4 POINTR

 1000 FORMAT(I2)

      DO 80 I=1,30
        NPLANT(I)=0
        PLTHT(I)=0.
        MSATYP(I)=0
        ANGBOLL(I)=0
        DO 80 J=1,450
          AVGPLT(I,J)=0.
   80 CONTINUE
      IF(PMAFIL.NE.'             ') THEN
        OPEN(11,FILE=PMAFIL,FORM=FILFRM,STATUS='OLD',ERR=340)
        I12 = 0
   90   CONTINUE

        I=0
        JPDAY=0
  100   CONTINUE
          IF((RUNMODE.EQ.'pmap').OR.(RUNMODE.EQ.'PMAP')) THEN
            CALL MEMPLANT
            IF((PMAFIL.EQ.'iamdone').OR.(PMAFIL.EQ.'IAMDONE')) &
                                                        GO TO 300
            PROFLE=PMAFIL
            IPLTHT=NFRQ
          ELSEIF((OPSYS.EQ.'dos').OR.(OPSYS.EQ.'DOS')) THEN
            READ(11,END=300,ERR=200) FILNAM,DSCRIP,MSDATE,IPLTHT, &
                NDMSMS,NDMSPF,NDV1MS,NDV1PF,NDV2MS, &
                NDV2PF,(PLTMAP(J),J=1,450),POINTR
            IF(I12.GT.0) READ(12,1000,END=120,ERR=120) NGBOLLMS
          ELSE
            READ(11,*,END=300,ERR=200) FILNAM,DSCRIP,MSDATE,IPLTHT,	&
                NDMSMS,NDMSPF,NDV1MS,NDV1PF,NDV2MS,	&
                NDV2PF,(PLTMAP(J),J=1,450),POINTR
            IF(I12.GT.0) READ(12,1000,END=120,ERR=120) NGBOLLMS
          ENDIF
  120     CONTINUE
          CALL CALTOJULAN(MSDATE,IDUM1,IDUM2,KDAY,J2DAY)
          IF(J2DAY.NE.JPDAY) THEN
            IF(I.GT.0) THEN
              NPLANT(I)=KNT
              PLTHT(I)=PLTHT(I)/KNT
              MSANODE(I)=MSANODE(I)/KNT
              IF((MSATYP(I).GT.0).AND.(KNT.GT.1)) THEN
                ANGBOLL(I)=ANGBOLL(I)/KNT
                DO J=1,450
                  AVGPLT(I,J)=AVGPLT(I,J)/KNT
                ENDDO
              ENDIF
            ENDIF
            I=I+1
            KNT=1
            JPDAY=J2DAY
            MSADTE(I)=J2DAY
            PLTHT(I)=IPLTHT
            MSANODE(I)=NDMSMS
            MSATYP(I)=0
            NODPMAP(I,1)=NDMSMS
            NODPMAP(I,2)=NDMSPF
            NODPMAP(I,3)=NDV1MS
            NODPMAP(I,4)=NDV1PF
            NODPMAP(I,5)=NDV2MS
            NODPMAP(I,6)=NDV2PF
            ANGBOLL(I)=NGBOLLMS
            IF(NDMSMS.GT.0) MSATYP(I)=1
            IF(I12.GT.0) MSATYP(I)=2
            DO 140 J=1,450
              IF((PLTMAP(J).EQ.'b').OR.(PLTMAP(J).EQ.'B').OR. &
                 (PLTMAP(J).EQ.'s').OR.(PLTMAP(J).EQ.'S').OR. &
                 (PLTMAP(J).EQ.'g').OR.(PLTMAP(J).EQ.'G')) THEN
                AVGPLT(I,J)=1
              ELSE
                AVGPLT(I,J)=0
              ENDIF
  140       CONTINUE
          ELSE
            KNT=KNT+1
            PLTHT(I)=PLTHT(I)+IPLTHT
            MSANODE(I)=MSANODE(I)+NDMSMS
            IF(NODPMAP(I,1).LT.NDMSMS) NODPMAP(I,1)=NDMSMS
            IF(NODPMAP(I,2).LT.NDMSPF) NODPMAP(I,2)=NDMSPF
            IF(NODPMAP(I,3).LT.NDV1MS) NODPMAP(I,3)=NDV1MS
            IF(NODPMAP(I,4).LT.NDV1PF) NODPMAP(I,4)=NDV1PF
            IF(NODPMAP(I,5).LT.NDV2MS) NODPMAP(I,5)=NDV2MS
            IF(NODPMAP(I,6).LT.NDV2PF) NODPMAP(I,6)=NDV2PF
            ANGBOLL(I)=ANGBOLL(I)+NGBOLLMS
            DO J=1,450
              IF((PLTMAP(J).EQ.'b').OR.(PLTMAP(J).EQ.'B').OR. &
                 (PLTMAP(J).EQ.'s').OR.(PLTMAP(J).EQ.'S').OR. &
                 (PLTMAP(J).EQ.'g').OR.(PLTMAP(J).EQ.'G')) THEN
                 AVGPLT(I,J)=AVGPLT(I,J)+1
              ENDIF
            ENDDO
          ENDIF
        GO TO 100
  200   CONTINUE
          PRINTBUF=' error on read in plant map adj. file'
          CALL WRITERR
          ABEND= .TRUE.
          CLOSE(11)
          RETURN
  300   CONTINUE
        IF(I.GT.0) THEN
          NPLANT(I)=KNT
          PLTHT(I)=PLTHT(I)/KNT
          MSANODE(I)=MSANODE(I)/KNT
          IF((MSATYP(I).GT.0).AND.(KNT.GT.1)) THEN
            ANGBOLL(I)=ANGBOLL(I)/KNT
            DO J=1,450
              AVGPLT(I,J)=AVGPLT(I,J)/KNT
            ENDDO
          ENDIF
        ENDIF
        CLOSE(11)
      ENDIF
  340 CONTINUE
      RETURN
      END


      SUBROUTINE VARIETYS

      use common_block

      CHARACTER VDSCRIP*51

  120 FORMAT(A8)

      OPEN(11,FILE=variety,FORM=FILFRM,STATUS='OLD',ERR=200)
  100 CONTINUE
      READ(11,*,ERR=200,END=300) VDSCRIP
      READ(11,120,ERR=200,END=300) VTYNAM
      READ(11,*,ERR=200,END=300) CALBRT
      VARITY(IVARTY) = VTYNAM
      CLOSE(11)
      RETURN

  200 CONTINUE
        PRINTBUF=' error on read in variety file'
        CALL WRITERR
        CLOSE(11)
        ABEND = .TRUE.
        RETURN
  300 CONTINUE
        PRINTBUF=' end-of-file while reading variety file'
        CALL WRITERR
        CLOSE(11)
        ABEND = .TRUE.
      RETURN
      END


      SUBROUTINE GCMDATA

      use common_block

      CHARACTER GCMNAME*51

      OPEN(11,FILE=gcmfil,FORM=FILFRM,STATUS='OLD',ERR=200)

      READ(11,*,ERR=200,END=300) GCMNAME

 100  continue
      READ(11,*,ERR=200,END=300) igcmonth,(gcm(I),I=1,6)
	  j=igcmonth
	  Radfactor(j)  = gcm(1)
      Tmaxfactor(j) = gcm(2)
	  Tminfactor(j) = gcm(3)
	  Rainfactor(j) = gcm(4)
	  RelHumfactor(j) = gcm(5)
	  Windfactor(j) = gcm(6)
	  if (j.ne.12 ) GO TO 100
      CLOSE(11)
      RETURN
  200 CONTINUE
        PRINTBUF=' error on read in GCM file'
        CALL WRITERR
        CLOSE(11)
        ABEND = .TRUE.
        RETURN
  300 CONTINUE
        PRINTBUF=' end-of-file while reading GCM file'
        CALL WRITERR
        CLOSE(11)
        ABEND = .TRUE.
      RETURN
      END


      SUBROUTINE NUMCHR(FILNAM,NUMBR)
      CHARACTER FILNAM*13
      I=0
   80 CONTINUE
        I=I+1
        IF(FILNAM(I:I).EQ.CHAR(0)) THEN
          DO 90 J=I,13
            FILNAM(J:J)=' '
   90     CONTINUE
        ENDIF
      IF(I.LT.13) GO TO 80
      I=0
      NUMBR=1
  100 CONTINUE
        I=I+1
        IF((FILNAM(I+1:I+1).EQ.'.').OR.(FILNAM(I+1:I+1).EQ.' ')) THEN
          NUMBR=I
        ELSE
          IF(I.LT.9) GO TO 100
          NUMBR=8
        ENDIF
      RETURN
      END


      SUBROUTINE PrintSummary

      use common_block

 1040 FORMAT(/	 &
      5X,'Profile Name:... ',A13,8x,'Start Simulation:', &
                                    I3,'/',I2,'/',I4/	 &
      5X,'Simulation Date: ',A15,6X,'Stop Simulation: ', &
                                    I3,'/',I2,'/',I4/	 &
      5X,'Variety:........ ',A8,13x,'Emergence Date:  ', &
                                    I3,'/',I2,'/',I4/	 &
      5X,'Description:.... ',A51/)
 1120 FORMAT('1',25X,'MAFES-GOSSYM (TM)  SUMMARY', /	&
             30X,'COPYRIGHT (C)  2003' / &
             14X,'Mississippi Agriculture & Forestry Experiment Station' /	&
             18X,'Contact:  Dr. K.R. Reddy and David W. Brand' /	&
             32X,'(662) 325-9462'/)
 1140 FORMAT(/													 &
       '                            PLANT           LIGHT', &
       '         GREEN   OPEN' /								 &
       '    EVENT         DATE DAE HEIGHT NODES LAI  INT.', &
       ' SQARS  BOLLS  BOLLS YIELD' /							 &
       '                            (in)             (%) ', &
       ' (-- 1000''S/ACRE --) (b/a)' / )
 1150 FORMAT(25x,'Soil Dispersion Trigger   ',i5,/)

      OPEN(22,FILE=sumryfle,STATUS='UNKNOWN')
      WRITE(22,1120) 

! ***	get the emergence,start and stop dates

      CALL JULANTOCAL(MO,DAZE,IYEAR,emerge)
      IMoemg = MO
      IDyemg= DAZE
	  iyremg= iyear
      CALL JULANTOCAL(MO,DAZE,IYEAR,jdstrs)
      IMostr = MO
      IDystr = DAZE
      iyrstr = iyear
      CALL JULANTOCAL(MO,DAZE,IYEAR,jdstps)
      IMostp = MO
      IDystp = DAZE
      iyrstp = iyear
      WRITE(22,1040) profle,imostr,idystr,iyrstr, &
	                 RUNDATE,imostp,idystp,iyrstp, &
                     VTYNAM,imoEMG,idyemg,iyremg,PDESCP
      WRITE(22,1150)idanoflg
      WRITE(22,1140) 
      WRITE(lstng,1040) profle,imostr,idystr,iyrstr, &
	                    RUNDATE,imostp,idystp,iyrstp,	 &
                        VTYNAM,imoEMG,idyemg,iyremg,PDESCP
      WRITE(lstng,1140) 
      RETURN
      END


      SUBROUTINE WRITERR

      use common_block

      CHARACTER*1 DUMBO

 1000 FORMAT(A80)
      OPEN(23,FILE=ERRFLE,STATUS='UNKNOWN')
  100 CONTINUE
        READ(23,*,END=200,ERR=100) DUMBO
      GO TO 100
  200 CONTINUE
      BACKSPACE(23)
      WRITE(23,1000) PRINTBUF
      CLOSE(23)
      RETURN
      END


      SUBROUTINE CALTOJULAN(MODYYR,MONTH,IDAY,IYEAR,JDAY)
!  ************************************************************
!  *                                                          *
!  *   DATE SUBROUTINE.  CONVERTS CALENDAR DAY TO JULIAN      *
!  *    DAY AND ALLOWS FOR LEAP YEARS.                        *
!  *                                                          *
!  ************************************************************
 
      CHARACTER*10 MODYYR
      DIMENSION I0(12)
 
      I0(1) = 0
      I0(2) = 31
      I0(3) = 59
      I0(4) = 90
      I0(5) = 120
      I0(6) = 151
      I0(7) = 181
      I0(8) = 212
      I0(9) = 243
      I0(10) = 273
      I0(11) = 304
      I0(12) = 334

      IF(MODYYR.NE.' ') THEN
        MONTH = (ICHAR(MODYYR(1:1))-48)*10 + (ICHAR(MODYYR(2:2))-48)
        IDAY = (ICHAR(MODYYR(4:4))-48)*10 + (ICHAR(MODYYR(5:5))-48)

!        IYEAR = (ICHAR(MODYYR(7:7))-48)*10+(ICHAR(MODYYR(8:8))-48)
! *** Changes to make the model read 4 digit year 

        IYEAR = (ICHAR(MODYYR(7:7))-48)*1000+(ICHAR(MODYYR(8:8))-48)*100 &
                + (ICHAR(MODYYR(9:9))-48)*10 + (ICHAR(MODYYR(10:10))-48)
      ENDIF

      IF((MONTH.LE.0).OR.(IDAY.LE.0)) THEN
        JDAY=0
      ELSE
        if((mod(iyear,4).eq.0).and.(mod(iyear,100).ne.0)) then
          do I=3,12
             I0(I) = I0(I) + 1
          enddo
        endif
        if((mod(iyear,100).eq.0).and.(mod(iyear,400).eq.0)) then
          do I=3,12
             I0(I) = I0(I) + 1
          enddo
        endif

!        IF(IYEAR/4*4.EQ.IYEAR) THEN
!          DO 100 I=3,12
!             I0(I) = I0(I) + 1
!  100     CONTINUE
!        ENDIF

        IF(MONTH.GT.12) MONTH=12
        JDAY = I0(MONTH) + IDAY
      ENDIF
      RETURN
      END

 
      SUBROUTINE JULANTOCAL(MONTH,IDAY,IYEAR,JULIAN)
!  ************************************************************
!  *                                                          *
!  *   DATE SUBROUTINE.  CONVERTS JULIAN TO CALENDAR AND      *
!  *    ALLOWS FOR LEAP YEARS.                                *
!  *                                                          *
!  ************************************************************
 
      INTEGER DACNT(12)

      DACNT(1) = 31
      DACNT(2) = 28
      if((mod(iyear,4).eq.0).and.(mod(iyear,100).ne.0)) DACNT(2) = 29
      if((mod(iyear,100).eq.0).and.(mod(iyear,400).eq.0)) DACNT(2) = 29

!      IF(IYEAR/4*4.EQ.IYEAR) DACNT(2) = 29

      DACNT(3) = 31
      DACNT(4) = 30
      DACNT(5) = 31
      DACNT(6) = 30
      DACNT(7) = 31
      DACNT(8) = 31
      DACNT(9) = 30
      DACNT(10) = 31
      DACNT(11) = 30
      DACNT(12) = 31

      MONTH = 1
      IDAY = JULIAN
      DO 100 I=1,12
        IF(IDAY.LE.DACNT(I)) GO TO 200
        MONTH = MONTH + 1
        IDAY = IDAY - DACNT(I)
  100 CONTINUE
  200 CONTINUE
      RETURN
      END

         
      SUBROUTINE GBLOCK

! ******** SIMULATED BLOCK DATA FOR GOSSYM **** DATE NOVEMBER 9, 1989 ********
! ******** CONVERTED BY WENDELL LADNER AND SUSAN BRIDGES FOR C ***************

      use common_block

      ABEND=.FALSE.
      ABZ=0.
      abz0=0.

! *** added (2) variables to differentiate rain from irrigation  gtheseira 

      actirrg = 0.
      actrain = 0.
  	  add60 = 0.0

      ADDEDN=0.
      AGETOP=0.
      AIRDRI=0.
      AIRDRW=0.
      ALPHA=3.5
      APRES=0.0
      AREA=0.
      AT=.2
      AVAILN=0.
      AVGTSD=0.
      AVGTSP=0.
      AVTEMP=20.

      BDI=0.
      BDRATO=0.
      BDSLOP=0.
      BDW=0.
      BETAI=0.
      BETAW=0.
      BOLL1=0.
      BURCN=0.
      BURMIN=0.
      BURR1=0.
      BURRN=0.

      CD=0.0
      CDBOLL=0.
      CDLEAF=0.
      CDROOT=0.
      CDSQAR=0.
      CDSTEM=0.
      CHARI = 'I'
      CMXIRR = 0.

! *** fertigation variables for FRTLIZ (mg N/cu. cm of irrigation water) 

      conamm = 0.
      connit = 0.
      conura = 0.

      COTXX=0.

!     CO2 = 0

      CPOOL=0.
      CSTRES=1.
      CSTORE=1.0
      CUMEP=0.
      CUMES=0.
      CUMRAN=0.
      CUMSOK=0.

      DCELL=5.
      DAY1PN=0.
      DAY1SN=0.
      DEAD2DAY=0.
      DEFBGN=0
      DEFDAY=0
      DEFKGH=0.
      DAYLNG=13.
      DAYNUM=1
      DAYWTF=0.
      DAZE=0.0
      DIFF0I=0.
      DIFF0W=0.
      DZ=0.0

	  dzprep = 1.0
	  duration = 0.0
	  daysq=0.
	  dayfb=0.
	  daymt=0.
	  daysnf=0.
	  daysnv=0.
	  dd60 = 0.0

      EMERGE=0.
      EP=0.
      ES=0.

      F2=.5
      FBLOOM=0.
      FCINIC=.4
      FCFCTI=0.
      FCINIW=.2
      FERN=0.
      FILFRM='formatted'
      IF((OPSYS.EQ.'dos').OR.(OPSYS.EQ.'DOS')) FILFRM='binary'
      FL = 0.
      FLOSS=0.0
      FLNMIN=1.0E-6
      FNH4=0.
      FNO3=1.
      FRATIO=0.
      FS = 0.
      FSQ=0.
      FSTRES=0.

	  frstbl = 0
	  frstob = 0
	  frstsq = 0
	  fsqfra=0.0

      GEOTR=10.
      GAMMA=0.653
      GBLOS=0.
      GBOLWT=0.
      GBZ2=0.
      GIN=0.
      GINP=0.0
      GSUBR=.375
	  growfac = 0.0

! *** TH2OADD replacements separates water infiltration on the left 
! *** from right of the plant 

      h2oaddl = 0.
      h2oaddr = 0.

      IDAY=0
      IFGIRR=0
      IFGRAIN=0
      INRIM=0
      INT=0
      IPIX=1
      ISCRN=6
      ISQ=0

	  isqday=0
	  iflday=0
	  igcmflg = 0.

	  iobday=0
	  imtday=0
!	  ipdays=0
	  irrflag=0

      KRAIN=0

! *** variables KULCLF and KULCRT are location ids of cultivated 
! *** cells under regular and skiprow configuration	

      kulclf = 7
      kulcrt = 14
      KWIDTH=0

      JDAY=1

      LAI=0.001
      LAMDAC=0.23
      LAMDAS=0.10
      LATUDE=35
      LDAYAW=0
      LDAYIR=0
      LDAYFW=0
      LDAYPW=0
      LEAFCN=.037
      LEAFR1=0
      LEAFRS=0
      leafwt=0.2
!      leafwt=0.02
      LEFABS=0
      LEFCNT=0
      LINE=51
      LMAX=0
      LR=5
      LTYPE=0
      LVSLOS=0
      LYTRES=0.

      MH2O=0
      MMUPN1=0.
      MMUPN2=0.
      MMUPN3=0.
      MO=1

!	mxfbrch = 30
!	mxvbrch = 3
!	mxfsite = 5
!	mxfruts	= 450

	  mxfbrch = 40
	  mxvbrch = 5
	  mxfsite = 6
	  mxfruts	= 1200
	  matday = 0.0

      NAPS=0
      NDLAY=0.
      NEWEP=0.0
      NEWES=0
      NF=0
      NFRQ=1
      NK=20
      NL=40
      NLOSS=0.
      NOITR=5
      NOPEN=0.
      NPOOL=0
      NUMPFN=1
      NR=0.
      NV=0.
      NVBRCH=1
      NYTTYM=0.
      NYTWTF=0

      n06=0

      ORGN=0.

      PDBOLL = 0.
      PDLEAF = 0.0001
      PDROOT = .0001
      PDSQ = 0.
      PDSTEM = .0001
      PI = 3.14159
      PLANTN = 0.
      PFAREA = 0.0
      PFDAL = 0. 
      PFDWL = 0.
      PIN = 0.0
      PLANTW = 0.
      PLEFABS = 0.
      PIXCON = 0.
      PIXLOS = 0.
      PIXPLT = 0.
      PIXDA = 1.
      PIXDN = 1.
      PIXDZ = 1.
	  PIXDPN=1.0
      PLTN = 0.
      PN = 0.
      POLYNA = 0
      POPFAC = 0.
      POPPLT = 41000.
      PQFLR = 0.
      PRPDAY = 0
      PRPKGH = 0.
      PSIAVG = -.175
      PSICMX = -0.5
      PSILD = -0.8
      PSILN = -0.8
      PSIMAX = 0.0
      PSINUM =  0.
      PTSRED = 1.

	  punits=0.0

      RAIN = 0.
      RCH2O = .0002
      REQ1 = 0.
      RESC = .06
      RESN = 0.
      RI = 0.
      RN = 0.
      ROOTCN = .037
      ROOTN = .00450
      ROOTR1 = 0.
      ROOTRS = 0.
      ROOTS = 0.
      ROOTWT = .200
      ROWSP = 101.6
      RSUBO = .0032
      RTP1 = .3
      RTP2 = .1
      RUTOFF = 0.

	  redlfcn = 1.0
	  rnfactor = 1.0

      SBOLL = 0.
      SDWBOL = 0.
      SDWLEF = 0.
      SDWSQR = 0.
      SDWSTM = 0.
      SEEDCN = 0.
      SEEDN = 0.
      SEEDR1 = 0.
      SEND = .FALSE.
      SESI = 0.
      SESII = 0.
      SITES = 0.
      SITEZ = 0.
      SKPFLG =  .FALSE.
      SLEAF = 0.
      SLEAFN = .0074
      SLF = .02
      SPDWBO = 0.
      SPDWLD = 0.
      SPDWLN = 0.
      SPDWRT = 0.
      SPDWSQ = 0.
      SPN = 0.
      SQUAR = 0.
      SQRZ = 0.
      SQWT = 0.
      SROOT = 0.
      SSTEM = 0.
      STEMCN = .037
      STEMN  = .0074
      STEMRS = 0.
      SUMEP = 0.0
      SUMES = 0.0
      SUMSTRS = 0.0
      SUPNO3 = 0.
      SUPNH4 = 0.
      SUMSUB = 0.
      SUBIRR = 0.
      SUPF = 0.


	  strdays = 1
      stemwt=0.2
!      stemwt=0.02
	  str01=0.0
	  str02=0.0
	  str03=0.0
	  str04=0.0
	  str05=0.0
	  str06=0.0
	  sumxirr = 0.0
	  sum_fsq_tavg = 0.0
	  sum_fbl_tavg = 0.0
	  sum_fob_tavg = 0.0
	  ave_fsq_tavg = 0.0
	  ave_fbl_tavg = 0.0
	  ave_fob_tavg = 0.0

      T = 0.
      TAVG = 0.
      TCELL = 1.
      TD = 0.
      TEMP1C = 0.
      TEMP1G = 0.
      TEMP1R = 0.
      TH2O = 0.
      TDAY = 0.
      TDFKGH = 0.
      THETAI = 0.
      THRLN = 0.3E-3
      THTA0I = 0.
      THTA0W = 0.
      THTARC = 0.
      THTARI = 0.
      THTARW = 0.
      THTASC = 0.
      THTASI = 0.
      THTASW = 0.
      TIH2OC = 0.
      TMAX = 0.
      TMIN = 0.
      TNNH4 = 0.
      TNNO3 = 0.
      TNH4UP = 0.
      TNO3UP = 0.
      TNYT = 0.

      UPNH4 = 0.
      UPNO3 = 0.

      V = 0.
      VARNAME = 'MISSING'
      VERSION = '06/20/03'
      VSTRES = 0.

      WCELL = 5.
      WATTBL = 200.
      WATTSM = 0.
      WIND = 88.
      WSTRSD = 1.
      WSTRSN = 1.
      WTSLFD = 0.
      WSTRS = 1.

      XTRAC = 0.
      XTRAN = 0.

      YIELD = 0.
      IYEAR = 74.
      
      Z = 3.0
!      Z = 0.1
      ZPIXD = 0.

! *** variable for leaf area development

 	 areapflf = 0.0
	 areamslf = 0.0
	 areafblf = 0.0

! ********************************************************************
!    INITIALIZATION OF ARRAYS IS LISTED IN ASCENDING ORDER BY THE
!    FIRST SUBSCRIPT ON THE ARRAY.  WITHIN AN INITIALIZATION OF 
!    ARRAYS OF THE SAME SIZE, THEY ARE LISTED ALPHABETICALLY.
! ********************************************************************

! *** LEAF AREA & LEAF WEIGHT INITIALIZED ACCORDING TO COTYLEDON   
! *** DATA FOR 'M-8' COTTON OF CHRISTIANSEN, M. N. (1962) A METHOD 
! *** OF MEASURING AND EXPRESSING EPIGNEOUS SEEDLING GROWTH RATE.  
! *** CROP SCI. 2:487-488.										   


! ***********************(3),(3,30),(3,30,5)**************************
      DO 100 I=1,mxvbrch
         VDELAY(I) = 0.
         NFBR(I)=0
        DO 100 J=1,mxfbrch
	      agemslf(i,j) = 0.0
	      xmslfage(i,j) = 0.0
	      xmslfarea(i,j) = 0.0
            DELAY(I,J)=0.
            MLAREA(I,J)=.04
            MLEAFW(I,J)=0.
            NNOD(I,J)=0
            PDAMLD(I,J) = 0.
            PDAMLN(I,J) = 0.
            PDWMLD(I,J) = 0.
            PDWMLN(I,J) = 0
            PRI(I,J) = '  '
            DO 100 K=1,mxfsite
               AGE(I,J,K)=0.
               AGEABZ(I,J,K)=0.0
               AGEBOL(I,J,K)=0.0
               AVGT(I,J,K)=0.
               BOLTMP(I,J,K)=0.
               BOLWGT(I,J,K)=0.
               BSIZE(I,J,K)=0.0
               DEHISS(I,J,K)=45.
               FCODE(I,J,K)=0.
               FFRUT(I,J,K)=0.
               FRUTP(I,J,K)=0.0
               FSTAVG(I,J,K)=0.5
               LAGE(I,J,K)=0.
               LAREA(I,J,K)=0.04
               LEAFW(I,J,K)=0
               MATURE(I,J,K)=0
               MCODE(I,J,K)=0
               PDADAY(I,J,K)=0.
               PDANYT(I,J,K) = 0.
               PDWBOD(I,J,K) = 0.
               PDWBON(I,J,K) = 0.
               PDWFLD(I,J,K) = 0.
               PDWFLN(I,J,K) = 0.
               PDWSQ(I,J,K) = 0.
               SQRWT(I,J,K) = 0.
               PRT(I,J,K) = '  '
	         agefrlf(i,j,k) = 0.0
	         frlfage(i,j,k) = 0.0
	         frlfarea(i,j,k) = 0.0
 100  CONTINUE

! *******************************(5)**********************************
      DO 110 I=1,5
         DEFDATE(I)=0
         DEFMTH(I)=0
         DEFPPA(I)=0.
         EPAVG(I)=.15
         KULDAY(I)=0
         FMTHOD(I)='      '
         HMTHOD(I)='      '
         IMTHOD(I)='      '
         PRPDATE(I)=0
         PRPMTH(I)=0
         PRPPPA(I)=0.
		 gcm(i) = 0.0
  110 CONTINUE
      FMTHOD(6)='      '
      FMTHOD(7)='      '
      HMTHOD(6)='      '
      HMTHOD(7)='      '
      IMTHOD(6)='      '
      IMTHOD(7)='      '
            
      
      FMTHOD(1)='BDCAST'
      FMTHOD(2)='SDRESS'
      FMTHOD(3)='FOLIAR'

! *** added (2) methods of fertilizer application:              
! ***  1) broadcast/incorporated and 2) fertigation	  gtheseira 

      FMTHOD(4)='bcstnc'
      FMTHOD(5)='frtgtn'
     
      HMTHOD(1)='BANDED'
      HMTHOD(2)='SPKLER'
      HMTHOD(3)='BDCAST'

      IMTHOD(1)='SPKLER'
      IMTHOD(2)='FURROW'
      IMTHOD(3)=' DRIP '

! ***  added (2) methods of irrigation application:			 
! ***  1) alternate drip and 2) alternate furrow   gtheseira 

      IMTHOD(4)='altfur'
      IMTHOD(5)='altdrp'


! *******************************(7)**********************************
      DO 120 I=1,7
         DTAVG(I)=20.
         PGRUNT(I)='      '
  120 CONTINUE
      PGRUNT(1)=' pts/a'
      PGRUNT(2)=' gal/a'
      PGRUNT(3)=' ozs/a'
      PGRUNT(4)=' lbs/a'
      PGRUNT(5)=' a/lb '
      PGRUNT(6)=' a/gal'

      CHAR1(1) = '-X'
      CHAR1(2) = '-*'
      CHAR1(3) = '-$'
      CHAR1(4) = '-A'
      CHAR1(5) = '-A'
      CHAR1(6) = '-A'
      CHAR1(7) = '-B'

      CHAR2(1) = 'X-'
      CHAR2(2) = '*-'
      CHAR2(3) = '$-'
      CHAR2(4) = 'A-'
      CHAR2(5) = 'A-'
      CHAR2(6) = 'A-'
      CHAR2(7) = 'B-'

! **************************(9)(9,40)******************************
       C1(1) = 0.3964D-0
       C1(2) = 0.3631D+1
       C1(3) = 0.3838D-1
       C1(4) = 0.7659D-1
       C1(5) = 0.0000D+0
       C1(6) = -0.2297D+2
       C1(7) = -0.3885D+0
       C1(8) = -0.1587D-0
       C1(9) = -0.01021D-1

      DO 130 I=1,9
         AIRDR(I)=0.
         BD(I)=0.
         BETA(I)=0.
         DIFF0(I)=0.
         FCININ(I)=.0
         FLXMAX(I)=0.
         FLXMIN(I)=0.
         GH2OC(I)=0.
         IPCLAY(I) = 30
         IPSAND(I) = 35
         LDEPTH(I) = 0
         THETA0(I) = 0.
         THETAR(I) = 0.
         THETAS(I) = 0.
         DO 130 J = 1,40
            TSTBD(I,J) = 0.
            TSTIMP(I,J) = 0.
  130 CONTINUE


! *******************************(10)*********************************
      DO 140 I=1,10
         AGEPFN(I)=0.
         AVTPFN(I)=0.
         PIXDAY(I) = 0
         PIXMTH(I) = 0
         PIXPPA(I) = 0.
         PFAL(I) = 0.04
         PFDALD(I) = 0.
         PFDALN(I) = 0.
         PFDWLD(I) = 0.
         PFDWLN(I) = 0.
         PFWL(I) = 0.
  140 CONTINUE

      DO I=1,10
	   agepflf(i) = 0.0
	   pflfage(i) = 0.0
	   pflfarea(i) = 0.0
	   cntpflfcn(i) = 0.0
	   sumpflfcn(i) = 0.0
 	enddo

! *******************************(12)*********************************

      KA(1) = ' '
      KA(2) = '0'
      KA(3) = '1'
      KA(4) = '2'
      KA(5) = '3'
      KA(6) = '4'
      KA(7) = '5'
      KA(8) = '6'
      KA(9) = '7'
      KA(10) = '8'
      KA(11) = '9'
      KA(12) = '*'

! *******************************(13)*********************************
      CHAR3(1) = '-1'
      CHAR3(2) = '-2'
      CHAR3(3) = '-3'
      CHAR3(4) = '-4'
      CHAR3(5) = '-5'
      CHAR3(6) = '-6'
      CHAR3(7) = '-7'
      CHAR3(8) = '-8'
      CHAR3(9) = '-9'
      CHAR3(10) = '-0'
      CHAR3(11) = '-A'
      CHAR3(12) = '-$'
      CHAR3(13) = '-X'
      CHAR4(1) = '1-'
      CHAR4(2) = '2-'
      CHAR4(3) = '3-'
      CHAR4(4) = '4-'
      CHAR4(5) = '5-'
      CHAR4(6) = '6-'
      CHAR4(7) = '7-'
      CHAR4(8) = '8-'
      CHAR4(9) = '9-'
      CHAR4(10) = '0-'
      CHAR4(11) = 'A-'
      CHAR4(12) = '$-'
      CHAR4(13) = 'X-'

! *******************************(14)*********************************
      DO 150 I=1,14
         H2OINT(I)=100.
  150 CONTINUE

      OMA(1)=1.
      RNNH4(1) = 50.
      RNNO3(1) = 10.
      DO 160 I=2,14
         OMA(I)=0.
         RNNH4(I) = 0.
         RNNO3(I) = 0.
  160 CONTINUE

! *******************************(15)*********************************
      DO 170 I = 1,15
         MSADTE(I) = 0
         DO 165 J=1,6
           NODPMAP(I,J)=0
  165    CONTINUE
         TAIR(I) = 0.
 170  CONTINUE

! *******************************(20)*********************************
      DO 180 I=1,20
         TSOLAV(I) = 0.
         PFNODAGE(I) = 0.
         PFNODLTH(I) = 0.
 180  CONTINUE

      BETAK(1)=1.1429E-4      
      BETAK(2)=1.1429E-4      
      BETAK(3)=1.1429E-4
      BETAK(4)=0.8601E-4
      BETAK(5)=0.8601E-4
      BETAK(6)=0.8601E-4
      BETAK(7)=0.6534E-4
      BETAK(8)=0.6534E-4
      BETAK(9)=0.6534E-4

      CONSK(1)=7.7699E-4
      CONSK(2)=7.7699E-4
      CONSK(3)=7.7699E-4
      CONSK(4)=5.7685E-4
      CONSK(5)=5.7685E-4
      CONSK(6)=5.7685E-4
      CONSK(7)=4.1618E-4
      CONSK(8)=4.1618E-4
      CONSK(9)=4.1618E-4

      DO 190 I=10,20
         BETAK(I)=0.4853E-4
         CONSK(I)=3.2446E-4
  190 CONTINUE

! *******************************(30)*********************************
      DO I=1,mxfbrch
	  XMNODAGE(I) = 0.
        XMNODLTH(I) = 0.
	  xnodage(i) = 0.0
	  cntlfcn(i) = 0.0
	  sumlfcn(i) = 0.0
      ENDDO

! ********************(40)(40,20)(40,20,3)(40,21)*********************
      DO 260 I=1,40
         BDL(I)=1.
         FC(I)=.267
         RTEXNT(I) = .FALSE.
         TSMN(I) = 25.
         TSMX(I) = 25.
         TSOILD(I) = 25.
         TSOILN(I) = 25.
         THTS(I) = 0.
         THTR(I) = 0.
         THAD(I) = 0.
	     cntmslfcn(i) = 0.0
	     summslfcn(i) = 0.0
         DO 270 J=1,20
            DIFF(I,J)=258.3
            KHAR(I,J) = ' '
            PSIS(I,J) =-.175
            ROOTSV(I,J) = 0.
            RTIMPD(I,J) = 0.
            RTWTCU(I,J) = 0.
            UPF(I,J) =0.
            VH2OC(I,J) = .267
            VNC(I,J) = 0.
            VNH4C(I,J) = 0. 
            VNO3C(I,J) = 0.
 270     CONTINUE   
         DO 290 J=1,21
            FNL(I,J)=0.
            TUPF(I,J) = .TRUE.
            TTUPF(I,J) = .TRUE.
 290     CONTINUE
 260  CONTINUE

! **************************(41)(41,10)(41,20)(41,21)*******************
      DO 440 I=1,41
         DO 450 J=1,10
            FWU(I,J)=0.0
  450    CONTINUE
         DO 460 J = 1,21
            PUPF(I,J) = 0.
  460    CONTINUE
         DO 465 J = 1,20
             FNU(I,J)=0.
  465    CONTINUE
  440 CONTINUE

! **************************(365)(365,7)(365,5)***********************
      DO 530 I=1,365
         AMTIRR(I) = 0.
         MTHIRR(I) = 0
         RUNOFF(I) = 0.
	   Bloom_tavg(i) = 0.0
         DO 540 J=1,7
            NFERT(I,J)=0
            CLIMAT(I,J)=0.
 540     CONTINUE
 530  CONTINUE

! ********************************(366)*******************************
      DO 610 I = 1, 366
         BOLOSS(I)=0.
         BLUM(I)=0.
         SQLOSS(I) = 0.
         STMWT(I) =  0.
 610  CONTINUE
      RETURN
      END


      SUBROUTINE INITIALIZE
!                    I N I T I A L I Z E
! ************************************************************
! *   INITIAL SETUP CALCULATIONS FOR THE MODEL               *
! *   Mostly soil parameters and CO2 function  Kit 1/20/2000 *
! ************************************************************
 
      use common_block
 
      WCELL = ROWSP/20.0
      ACELLDW = DCELL * WCELL
      VCELL = DCELL * WCELL * TCELL
      LPLOW = 20 / DCELL
      THRLN = THRLN*(ACELLDW/25.)
      POPFAC = 404685.6/POPPLT
 
! *** SAVE INITIAL SOIL VARIABLE 
 
      DIFF0I = DIFF0(1)
      THTA0I = THETA0(1)
      BETAI = BETA(1)
      THTASI = THETAS(1)
      THTARI = THETAR(1)
      AIRDRI = AIRDR(1)
      BDI = BD(1)
      FCFCTI = FCININ(1)
      K = 9
      DO 100 J=1,8
        DIFF0(K) = DIFF0(K-1)
        THETA0(K) = THETA0(K-1)
        BETA(K) = BETA(K-1)
        LDEPTH(K) = LDEPTH(K-1)
        THETAS(K) = THETAS(K-1)
        FCININ(K) = FCININ(K-1)
        THETAR(K) = THETAR(K-1)
        AIRDR(K) = AIRDR(K-1)
        BD(K) = BD(K-1)
        IPSAND(K) = IPSAND(K-1)
        IPCLAY(K) = IPCLAY(K-1)
        K = K - 1
  100 CONTINUE
      LYRSOL=LYRSOL+1
      LDEPTH(1) = 5.05
      KULKNT = 0

! *** DM**2 GROUND AREA/PLANT 

      DELT = 1./NOITR
      KUPPER = DCELL*(NK+1)-KWIDTH
      J=1
      DO 200 L = 1,NL
  120   CONTINUE
        IF((L*DCELL.GT.LDEPTH(J)).AND.(J.LT.9)) THEN
           J = J+1
           GO TO 120
        ENDIF
        IF(J.GT.LYRSOL) J=LYRSOL
        LYRDPH(L) = J
        ARDRCN(L) = ALOG(-15./PSISFC) / ALOG((THETAR(J)-AIRDR(J)) /	&
                                       (FCININ(J)-AIRDR(J)))

! *** FLXMAX AND FLXMIN ARE MAX AND MIN WATER FLUXES. CM**2			
! *** FC(L) AND THTS(L) ARE FIELD CAPASITY AND SATURATED VOLUMETRIC 
! *** WATER CONTENTS, RESPECTIVELY, BY LAYER.   					

        FLXMAX(J)=ABS(DIFF0(J)*((FCININ(J)-THETAR(J))/DCELL)*  &
                     (WCELL*DELT)*EXP(BETA(J)*(FCININ(J)-THETA0(J))))
        FC(L) = FCININ(J)
        BDL(L) = BD(J)
        THTS(L) = THETAS(J)
        THTR(L) = THETAR(J) 
        THAD(L) = AIRDR(J)
        FLXMIN(J) = ABS(DIFF0(J) * ((THETA0(J)-THETAR(J))/DCELL) &
                                                    *WCELL*DELT)
        IF(FLXMIN(J).LE.0.0001) FLXMIN(J) = 0.0001
  200 CONTINUE
 
! *** Initialize root growth variables    kit jun 20, 1996 
  
      DO 210 I=1,40
        DO 210 J=1,20
          DO 210 K=1,3
            RTWT(I,J,K)=0

! *** If skiprow locate plant 1/3 way into slab.  Test is same as
! *** calculation of ROWSP in PROFILE subroutine        gtheseira  

            IF((SKIPWD.GE.1.).AND.(I.EQ.1).AND.(J.EQ.7).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.035/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.2).AND.(J.EQ.7).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.025/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.3).AND.(J.EQ.7).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.015/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.4).AND.(J.EQ.7).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.010/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.5).AND.(J.EQ.7).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.005/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.1).AND.(J.EQ.6).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.010/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.1).AND.(J.EQ.8).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.035/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.2).AND.(J.EQ.8).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.025/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.3).AND.(J.EQ.8).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.015/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.4).AND.(J.EQ.8).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.010/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.5).AND.(J.EQ.8).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.005/POPFAC
            IF((SKIPWD.GE.1.).AND.(I.EQ.1).AND.(J.EQ.9).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.010/POPFAC

! *** Plant is located in middle of slab 

            IF((SKIPWD.LT.1.).AND.(I.EQ.1).AND.(J.EQ.10).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.035/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.2).AND.(J.EQ.10).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.025/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.3).AND.(J.EQ.10).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.015/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.4).AND.(J.EQ.10).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.010/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.5).AND.(J.EQ.10).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.005/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.1).AND.(J.EQ.9).AND.(K.EQ.1))  &
            RTWT(I,J,K)=.010/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.1).AND.(J.EQ.11).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.035/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.2).AND.(J.EQ.11).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.025/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.3).AND.(J.EQ.11).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.015/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.4).AND.(J.EQ.11).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.010/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.5).AND.(J.EQ.11).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.005/POPFAC
            IF((SKIPWD.LT.1.).AND.(I.EQ.1).AND.(J.EQ.12).AND.(K.EQ.1)) &
            RTWT(I,J,K)=.010/POPFAC
  210 CONTINUE

! *** Initialization of: 1) 'extent of rooting' variables to left and 
! *** right of plant and 2) location variables of first cultivated    
! *** cell to left, KULCLF, and to right, KULCRT of plant under       
! *** regular and skiprow configuration.s                  gtheseira  

      IF(SKIPWD.GE.1.0) THEN
        KLL(1)=6
        DO I=2,40
           KLL(I)=7       
        enddo
        KRL(1)=9
        DO  I=2,40
           KRL(I)=8
        enddo
        KULCLF=4
        KULCRT=11
      ELSE 
        KLL(1)=9
        DO I=2,40
           KLL(I)=10
        enddo
        KRL(1)=12
        DO I=2,40
           KRL(I)=11
        enddo
        KULCLF=7
        KULCRT=14
      ENDIF

! *** SOILNO OF 0.06 IS POT. MIN. NIT. FOR GROUP 2 SOILS (FINE LOAMY/CLAY)
! *** and 0.10 FOR GROUP 1 SOILS (COARSE LOAMY/LOAM)

      IF(PSISFC.LT.-0.2) THEN
         SOILNO = 0.06
      ELSE
         SOILNO = 0.10
      ENDIF
      DO 220 L=1,NL
        VNH4C(L,1)=0.
        VNO3C(L,1)=0.
        VNC(L,1)=0.
  220 CONTINUE
      DO 240 I=1,197
        J=(I-1)/15+1
        IF(J.GT.14) J=14
        L=(I-1)/DCELL+1
        IF(L.GT.NL) L=NL

! *** INITIALIZE NO3 AND NH4 ARRAY WITH INITIAL RESIDUAL NO3 OR WITH 
! *** 4.0 LBS/6 INCHES LAYER

        IF(RNNO3(J).LE.0.) RNNO3(J) = 2.0
        IF(RNNH4(J).LE.0.) RNNH4(J) = 0.2
        VNO3C(L,1) = VNO3C(L,1)+RNNO3(J)/15.*.011219
        VNH4C(L,1) = VNH4C(L,1)+RNNH4(J)/15.*.011219

! *** CONVERT VNC FROM % OM TO MG/CM**3.  1000 mg = 1 g 
! *** POTN IS POTENTIAL NINERALIZABLE N IN ORGANIC MATTER. I IS DEPTH, cm

       POTN = (AMAX1(0.0, 0.055 - 0.0241*LOG10(FLOAT(I))))*1.1
       POTN = AMIN1(0.05, POTN)
       VNC(L,1)   = VNC(L,1)+(OMA(J)/100.)*BDL(L)*1000.*POTN*SOILNO
  240 CONTINUE
      DO 260 L=1,NL
        VNO3C(L,1) = VNO3C(L,1)/DCELL
        VNH4C(L,1) = VNH4C(L,1)/DCELL
        VNC(L,1)   = VNC(L,1)/DCELL
        DO 260 K=2,NK
          VNO3C(L,K) = VNO3C(L,1)
          VNH4C(L,K) = VNH4C(L,1)
          VNC(L,K)   = VNC(L,1)
  260 CONTINUE
      J=1
      DO 300 L=1,NL
        IF((L-1)*DCELL.GT.J*15) THEN
          J=J+1
          IF(J.GT.14) J=14
        ENDIF
        N=LYRDPH(L)
        DO 280 K=1,NK
          IF((L*DCELL).LT.WATTBL) THEN
            VH2OC(L,K)=FC(L)*H2OINT(J)/100.
          ELSE
            VH2OC(L,K)=THETAS(N)
          ENDIF
          IF(VH2OC(L,K).LT.THETA0(N)) VH2OC(L,K) = THETA0(N)
          TIH2OC = TIH2OC + VH2OC(L,K)
          DIFF(L,K)=DIFF0(N)*EXP(BETA(N)*(VH2OC(L,K)-THETA0(N)))
          TEMP1G = (VH2OC(L,K)-AIRDR(N))/(FCININ(N)-AIRDR(N))
          PSIS(L,K)= PSISFC * TEMP1G**ARDRCN(L)
  280   CONTINUE
  300 CONTINUE
      TBL = NL * DCELL
      IF(WATTBL.GT.TBL) THEN
        THETAI = FCININ(LYRSOL)*H2OINT(14)/100.
        IF(THETAI.LE.0.0)THETAI = FCININ(LYRSOL)
      ENDIF
	count_of_columns = dfloat(nk)
      TIH2OC = (TIH2OC * DCELL * 10.)/count_of_columns
      TNNO3 = 0.
      TNNH4 = 0.
      DO 320 L=1,NL
         DO 320 K=1,NK
            TNNO3 = TNNO3 + VNO3C(L,K)
            TNNH4 = TNNH4 + VNH4C(L,K)
 320  CONTINUE
      TNNO3 = TNNO3 * vcell/rowsp/0.011219
      TNNH4 = TNNH4 * vcell/rowsp/0.011219
      DAY1SN = TNNO3 + TNNH4
      DAY1PN = ROOTN + STEMN + SLEAFN + SEEDN + BURRN
      WTDAY1 = ROOTWT+STEMWT+GBOLWT+LEAFWT+SQWT+XTRAC+COTXX+RESC
 
! *** GET THE CO2 CORRECTION FACTOR FOR PHOTOSYNTHESIS.
 
	if(co2.eq.0.0) then
         IF(IYEAR.GT.1900) INDXCO2 = IYEAR - 1959
         IF(INDXCO2.LT.1) INDXCO2 = 1
         PNETCOR = 1.0208594 + 0.0021710*INDXCO2 + 0.0000717*INDXCO2**2
	else

! *** New Carbon vs CO2 function. Assume base CO2=320 ppm  Reddy 11/99

	   stdco2 = 320.0
	   co2num = (0.0179*co2*5.7194)/(0.0179*co2 + 5.7194) 
	   co2den = (0.0179*stdco2*5.7194)/(0.0179*stdco2 + 5.7194) 
	   pnetcor = co2num/co2den
	endif

      RETURN
      END


      SUBROUTINE CLYMAT
! ************************************************************
! *                                                          *
! *                  CLIMATE SUBROUTINE                      *
! *                                                          *
! ************************************************************

      use common_block

      DAYNUM = JDSTRS+IDAY-1
      RI   = CLIMAT(DAYNUM,1)
      TMAX = (CLIMAT(DAYNUM,2)-32.) * .5555556
      TMIN = (CLIMAT(DAYNUM,3)-32.) * .5555556
      WIND = CLIMAT(DAYNUM,6)
      RAIN = (CLIMAT(DAYNUM,5) + AMTIRR(DAYNUM)) * 25.4

! ***  new variables to distinguish rain from irrigation  gtheseira 

      actrain = climat(DAYNUM,5) * 25.4
      actirrg = amtirr(DAYNUM) * 25.4

! *** POLLINATION SWITCH DEPENDS ON RAIN.

      IF(RAIN.GT.12.7) THEN
        POLYNA = 0
        IF(AMTIRR(DAYNUM).GT.0.1) THEN
          J = MTHIRR(DAYNUM)+1
          IF(IMTHOD(J).NE.'SPKLER') POLYNA = 1
        ENDIF
      ELSE
        POLYNA = 1
      ENDIF
      IF(AMTIRR(DAYNUM).GT.0.) LDAYIR = DAYNUM
      CALL RRUNOFF

! ***  lower limit of effective rainfall 

      if(actrain.le.1.5)actrain = 0.0
      IF(RAIN.LE.1.5) RAIN = 0.0
      IF(TMAX.LE.0) TMAX = 2.0
      IF(TMIN.LE.0) TMIN = 1.0

! *** CALCULATE SOLAR DECLINATION USING D. W. STEWART'S PROGRAM
! *** WHICH INCORPORATES THE ALGORITHM OF ROBINSON AND RUSSELO.

      DEGRAD=3.1416/180.0
      XLAT=LATUDE*DEGRAD
      DEC=C1(1)
      DO 2 I=2,5
        N=I-1
        J=I+4
        PHI=N*0.01721*DAYNUM
    2 DEC=DEC+C1(I)*DSIN(PHI)+C1(J)*DCOS(PHI)
      DEC=DEC*DEGRAD

! *** CALCULATE DAYLENGTH. 

      DAYLNG=DACOS((-0.014544-DSIN(XLAT)*DSIN(DEC))/(DCOS(XLAT) *  &
             DCOS(DEC)))*7.6394

      WATTSM = RI * 697.45 / (DAYLNG*60.)
      RN = WATTSM * .8 - 26.

! *** RN = NET RADIATION IN WATTS/M**2	

      JOULES = RI*41860.

! *** JOULES = NET RADIATION IN JOULES/M**2

      TMINT  = (CLIMAT(DAYNUM+1,3)-32.) * .5555556
      if(igcmflg.ne.0) then
         CALL JULANTOCAL(MO,DAZE,IYEAR,daynum)
         j = MO
         TMINT  = TMINT	+ Tmaxfactor(j)
	  endif
      IF(TMINT.LE.0.0) TMINT = 1.0
      IF(RAIN.LE.2.54) KRAIN = 0
      IF(RAIN.GT.2.54.AND.RAIN.LE.25.4) KRAIN = 1
      IF(RAIN.GT.25.4) KRAIN = 2

! *** THE FOLLOWING STATEMENTS CALCULATE TDAY AND TNYT  B.ACOCK 

      WATACT = PI * JOULES/3600./2./DAYLNG
      X5 = 0.0945-(WATACT*8.06*10.0**(-5))+(TMAX*6.77*10.0**(-4))
      TMX5WA = TMAX/X5/WATACT
      IF(TMX5WA.GE.1.)TMX5WA = 1.0
      IF(TMX5WA.LE.-1.)TMX5WA = -1.0
      TMAXHR = DAYLNG/PI*(PI-(ASIN(TMX5WA)))

      X6 = ((TMAX-TMIN)/2.)
      X7 = (PI/TMAXHR)
      X8 = 1.5 * PI
      TDAY=X6*(1-(TMAXHR/PI/DAYLNG*COS((X7*DAYLNG)+X8)))+TMIN

! *** CALCULATE AIR TEMPERATURE AT DUSK. 

      TAIRSS=(X6*(1.0+SIN((X7*DAYLNG+X8))))+TMIN

! *** CALCULATE AVERAGE NIGHT-TIME AIR TEMPERATURE. 

      IF(TAIRSS-TMINT.LE.0.5)  TAIRSS = TMINT + 0.5
      TNYT=(TAIRSS-TMINT)/(-ALOG(TMINT/TAIRSS))
      CALL JULANTOCAL(MO,DAZE,IYEAR,DAYNUM)
      TAVG=(TDAY*DAYLNG+TNYT*(24.-DAYLNG))/24.

! *** LTYPE=1 OKRA LEAF.   LTYPE=0 NORMAL LEAF
! *** TEST OF CHANGE IN CANOPY LIGHT INTERCEPTION DUE TO PIX   06/29/90

      ZDUM = Z + ZPIXD * .6
      IF(LTYPE.EQ.1) THEN
         INT=(-2.05595+1.64301*ZDUM+(-.00648851*(ZDUM**2)))/100.
      ELSE
         INT = 1.0756*ZDUM/ROWSP
      ENDIF
      IF(INT.LT.0.) INT = 0.
      IF(INT.GE.0.95) INT = 0.95

! *** INT = FRACTION OF INCIDENT LIGHT INTERCEPTED BY PLANT CANOPY.
! *** BAKER ET. AL. CANOPY ARCHITECTURE IN RELATION TO YIELD.
! *** CHAPTER 3 IN 'CROP PHYSIOLOGY' ED. V. S. GUPTO.

      IF(LAI.GT.LMAX) LMAX = LAI
      IF(LAI.LT.LMAX.AND.LAI.LT.3.1) THEN
         CLAI = 0.0
         IF(Z.LT.ROWSP) THEN
            CLAI = 3.1 * Z/ROWSP
            IF(LAI.LE.CLAI) INT = INT * LAI / CLAI
         ELSE
            INT=INT*LAI/3.1
         ENDIF
      ENDIF
      IF(INT.LT.0.) INT = 0.
      IF(INT.GE.0.95) INT = 0.95
      CALL TMPSOL
      IF(KDAY.EQ.1) AVTPFN(1) = TAVG
      RETURN
      END


      SUBROUTINE RRUNOFF
! ***********************************************************
! *                                                         *
! *   RUNOFF SUBROUTINE.  CALCULATES PORTION OF RAINFALL    *
! *   THAT IS LOST TO RUNOFF, REDUCES RAINFALL TO THAT      *
! *   WHICH IS INFILTRATED IN THE SOIL.  SUBROUTINE USES    *               
! *   THE SOIL CONSERVATION SERVICE METHOD OF ESTIMATING    *
! *   RUNOFF.                                               * 
! *   REFERENCE:  BRADY, NYLE C. 1984. THE NATURE AND       *
! *   PROPERTIES OF SOILS, 9TH ED. MACMILLAN PUBLISHING CO. *
! *   REFERENCE:  SCHWAB, FREVERT, EDMINSTER, AND BARNES.   *
! *   1981. SOIL AND WATER CONSERVATION ENGINEERING, 3RD    *
! *   JOHN WILEY & SONS, INC.                               *
! *                                                         *
! ***********************************************************

      use common_block

      DIMENSION D02(3,3)
      IDUM=IFGRAIN
      IF(DAYNUM.EQ.LDAYIR) IDUM=IFGIRR
      IF(IDUM.EQ.0) THEN
      
! ***  LOOP TO ACCUMULATE 5-DAY ANTECEDENT RAINFALL (MM) WHICH
! *** WILL AFFECT THE SOIL'S ABILITY TO ACCEPT NEW RAINFALL.
! *** ANTECEDENT RAINFALL INCLUDES RAINFALL AND ALL IRRIGATION. 
! *** D11 = 5-DAY TOTAL 

        D11 = 0.0
        I01 = DAYNUM-6
        IF(I01.LT.1) I01=1
        I02=DAYNUM-1
        IF(I02.LT.1) I02=1
        DO 10 K = I01,I02
          D11 = D11+(AMTIRR(K)+CLIMAT(K,5))*25.4
   10   CONTINUE

! *** I03 IS AN INDICATOR OF THE ANTECEDENT MOISTURE CONDITIONS.
! *** 1=CONDITION I, LOW MOISTURE, LOW RUNOFF POTENTIAL. 
! *** 3=CONTITION III, WET CONDITIONS, HIGH RUNOFF POTENTIAL.

        IF (D11.LT.36)THEN 
          I03=1
        ELSEIF (D11.GT.53)THEN 
          I03=3
        ELSE
          I03=2 
        ENDIF

! *** CALCULATE THE FINAL INFILTRATION RATE FOR Ap HORIZON
! *** FINAL INFILTRATION IS ESTIMATED FROM THE PERCENT SAND
! *** AND PERCENT CLAY IN THE Ap LAYER. IF CLAY CONTENT IS  
! *** GREATER THAN 40%, THE SOIL IS ASSUMED TO HAVE A HIGHER 
! *** RUNOFF POTENTIAL (SOIL GROUP C=3). IF CLAY CONTENT IS 
! *** LESS THAN 15% AND SAND IS GREATER THAN 70%, A LOWER 
! *** RUNOFF POTENTIAL IS ASSUMED (SOIL GROUP A=1).	OTHER 
! *** SOILS (LOAMS) ASSUMED MODERATELY LOW RUNOFF POTENTIAL
! *** (SOIL GROUP B=2). NO 'IMPERMEABLE' (GROUP D) SOILS ARE            
! *** ASSUMED.                     REFERENCES: SCHWAB, BRADY.      

        IF(IPSAND(1).GT.70.AND.IPCLAY(1).LT.15)THEN
          I04 = 1
        ELSEIF(IPCLAY(1).GT.35)THEN
          I04 = 3
        ELSE 
          I04=2   
        ENDIF

! *** ASSUME STRAIGHT ROW, GOOD PRACTICE CROPPING PRACTICE.
! *** RUNOFF CURVE NUMBER, UNADJUSTED FOR MOISTURE AND SOIL TYPE. 

        CRVNUM=78.0 

! *** ADJUST CURVE NUMBER FOR SOIL GROUP A,B,C. COEFFICIENTS ADJUST
! *** FOR SOIL GROUPS OTHER THAN 'B', WHICH IS ASSUMED IN THE SCS EQ.
 
        IF(I04.EQ.3)THEN
          D01=1.14
        ELSEIF(I04.EQ.2)THEN
          D01=1.09
        ELSEIF(I04.EQ.1)THEN
          D01=1.0
        ENDIF
 
! *** ADJUST CURVE NUMBER FOR ANTECEDENT RAINFALL CONDITIONS I,II,III.

        D02(1,1)=0.78
        D02(1,2)=0.83
        D02(1,3)=0.87
        D02(2,1)=1.00
        D02(2,2)=1.00
        D02(2,3)=1.00
        D02(3,1)=1.15
        D02(3,2)=1.10
        D02(3,3)=1.07
    
        CRVNUM = CRVNUM*D01*D02(I03,I04)

! *** EFFECTIVE RAINFALL = RAINFALL (OR IRRIGATION) - RUNOFF.   
! *** D03 = MAX POTENTIAL DIFFERENCE BETWEEN RAINFALL AND RUNOFF.

        D03 = 25400./CRVNUM - 254.0
        RUNOFF(DAYNUM) = ((RAIN-0.2*D03)**2)/(RAIN+0.8*D03)
        IF (RAIN.LE.0.2*D03) RUNOFF(DAYNUM)=0.0
        RAIN = RAIN-RUNOFF(DAYNUM)

! *** if there is a rainfall event while irrigation is on
! *** rain runs off							  gtheseira
    
        if(actrain.le.runoff(DAYNUM))then
           actirrg = actirrg + actrain - runoff(DAYNUM)
           actrain = 0.0
        else
           actrain = actrain - runoff(DAYNUM)
        endif

      ENDIF
      RETURN
      END


      SUBROUTINE TMPSOL
! ***********************************************************
!   THIS SUBROUTINE CALCULATES A TEMPERATURE PROFILE IN THE *
!  SOIL. ASSUMES HORIZONTAL HOMOGENEITY OF TEMPERATURE &    *
!  DISREGARDS MOISTURE CONTENT EFFECTS. FIRST, MAXIMUM (H)  *
!  & MINIMUM (L) TEMPERATURES ARE CALCULATED AT 2, 4, 8, &  *
!  16 INCHES DEPTHS BY MULTIPLE REGRESSION EQUATIONS OF     *
!  J.C. MCWHORTER & B.P. BROOKS, JR.  1965. CLIMATOLOGICAL  *
!  AND SOLAR RADIATION RELATIONSHIPS. BULL. 715, MS AGRI.   *
!  EXP. STA., STARKVILLE.  NOTE THAT THE GRID SIZE (D*W)    *
!  IS NOT VARIABLE IN THIS SUBROUTINE, BUT THE LAYER        *
!  THICNESS IS FIXED AT 5 CM.  MAX & MIN SOIL TEMPS FOR     *
!  EACH OF THE LAYERS ARE THEN OBTAINED BY INTERPOLATION &  *
!  EXTRAPOLATION OF THE 2, 4, 8, & 16 INCH TEMPS.           *
!    FINALLY, DAYTIME AND NIGHTIME TEMPS(TSMX & TSMN) ARE   *
!  OBTAINED AS AVERAGE HOURLY VALUES FROM 7 A.M. THRU       *
!  SUNSET, & SUNSET THRU 7 A.M., RESPECTIVELY, USING AN     *
!  ALGORITHM FOR AIR TEMP PUBLISHED BY H. N. STAPLETON,     *
!  D.R. BUXTON, F.L. WATSON, D.J. NOLTING, AND D.N. BAKER.  *
!  UNDATED.  COTTON: A COMPUTER SIMULATION OF COTTON        *
!  GROWTH.  TECH. BULL. 206, AZ AGRI. EXP. STA. TUCSON.     *
! ***********************************************************
  
      use common_block
 
      DO 1 I = 1,6
      J = 8 - I
      JM1 = J - 1
    1 DTAVG(J) = DTAVG(JM1)
      DTAVG(1) = TAVG
      WTAVG = 0.
      DO 2 J = 1,7
    2 WTAVG = WTAVG +  DTAVG(J)
      WTAVG = WTAVG/7.
      WTAVGF = WTAVG*1.8 + 32.

! *** THE NEXT EIGHT EQUATIONS ARE FROM MCWHORTER AND BROOKS.

        T2H  = 1.1962*WTAVGF + 0.27389
        T2L  = 0.960*WTAVGF  + 1.4404
        T4H  = 1.1493*WTAVGF + 1.1452
        T4L  = 0.9126*WTAVGF + 2.9961
        T8H  = 0.9655*WTAVGF + 8.3121
        T8L  = 0.8700*WTAVGF + 7.9217
        T16H = 0.8409*WTAVGF + 13.988
        T16L = 0.8341*WTAVGF + 13.029

! *** GET TEMP OF SOIL (MAX) BY INTERPOLATION OR EXTRAPOLATION.

      T24 = T2H - T4H
      T48 = T4H - T8H
      TSMX(1) = T2H + (.507874) * T24
      TSMX(2) = T4H + (.523622) * T24
      TSMX(3) = T8H + (.769685) * T48
      TSMX(4) = T8H + (.277559) * T48
      T816 = .0492126 * (T8H - T16H)
      DO 6 I=5,20
        TSMX(I) = T8H - (2.18+(I-5)*5.) * T816
 6    CONTINUE

! *** GET TEMP OF SOIL (MIN) BY INTERPOLATION OR EXTRAPOLATION.

      T24 = T2L - T4L
      T48 = T4L - T8L
      TSMN(1) = T2L + (.507874) * T24
      TSMN(2) = T4L + (.523622) * T24
      TSMN(3) = T8L + (.769685) * T48
      TSMN(4) = T8L + (.277559) * T48
      T816 = .0492126 * (T8L - T16L)
      DO 7 I=5,20
        TSMN(I) = T8L - (2.18+(I-5)*5.) * T816
        IF(TSMN(I).LT.TSMX(I)) GO TO 7
        TSMN(I) = (TSMN(I) + TSMX(I))/2.
        TSMX(I) = TSMN(I)
 7    CONTINUE
      DO 8 I=1,20

! *** CONVERT TEMPS TO CENTIGRADE.

        TSMX(I) = (TSMX(I)-32.)*.555556
        TSMN(I) = (TSMN(I)-32.)*.555556
 8    CONTINUE
      ISR = 12 - IFIX(DAYLNG*.5)
      ISS = ISR + IFIX(DAYLNG+0.5)

! *** HOUR OF SUNSET. EQUATIONS DETERMINING RECDAT pp.37 of Stapleton, et.al.

      DO 9 LAYER = 1,20
      TMEAN = (TSMX(LAYER)+TSMN(LAYER)) * .5
      SWINGH = (TSMX(LAYER)-TSMN(LAYER)) * .5
      DO 11 IH=7,15
      RECDAT(IH) = TMEAN - SWINGH*COS(0.3927*(IH-7.))
      IH9 = IH + 9
      RECDAT(IH9) = TMEAN + SWINGH*COS(0.19635*(IH9-15.))
 11   CONTINUE
      DO 12 IH=1,6
 12   RECDAT(IH) = TMEAN - SWINGH*COS(0.19635*(7-IH))
      SHRTD = 0.
      SHRTN = 0.
      DO 13 IH=7,ISS
      SHRTD = SHRTD + RECDAT(IH)

! *** SUM OF HOURLY TEMPS IN DAYTIME.

 13   CONTINUE
      TSOILD(LAYER) = SHRTD/(ISS-6)

! *** AVERAGE TEMP OF SOIL DURING DAYTIME, DEG C.

      ISS1 = ISS + 1
      DO 14 IH=ISS1,24
      SHRTN = SHRTN + RECDAT(IH)

! *** SUM OF HOURLY TEMPS IN NIGHTIME.

 14   CONTINUE
      DO 15 IH=1,6
      SHRTN = SHRTN + RECDAT(IH)
 15   CONTINUE
      TSOILN(LAYER) = SHRTN/(30-ISS)

! *** AVERAGE TEMP OF SOIL DURING NIGHTIME.

 9    CONTINUE
      NLH = NL/2
      DO 16 LAYER = 1, NLH
        IF(LAYER.GT.10.AND.TSOILD(LAYER).LT.22.) TSOILD(LAYER)=22.0
        IF(LAYER.GT.10.AND.TSOILN(LAYER).LT.20.) TSOILN(LAYER)=20.0
        IF(LAYER.GT.15.AND.TSOILD(LAYER).LT.25.) TSOILD(LAYER)=25.0
        IF(LAYER.GT.15.AND.TSOILN(LAYER).LT.25.) TSOILN(LAYER)=25.0
        TSOLAV(LAYER)= (TSOILD(LAYER)*DAYLNG+TSOILN(LAYER)*	 &
                       (24.-DAYLNG))/24.
 16   CONTINUE
      DO 20 LAYER = 1,NL
        DO 20 K=1,20
          IF(LAYER.LE.20) THEN
            SOILT(LAYER,K)=TSOLAV(LAYER)
          ELSE
            SOILT(LAYER,K)=25.0
          ENDIF
 20   CONTINUE

! *** AVERAGE SOIL TEMPERATURE, DEG C.

      RETURN
      END


      SUBROUTINE SOIL
! ************************************************************
! *                                                          *
! *     SOIL SUBROUTINE.  CALLS FRTLIZ, GRAFLO, ET,          *
! *       UPTAKE, CAPFLO, AND NITRIF.                        *
! *                                                          *
! ************************************************************

      use common_block
 
      CALL FRTLIZ

      DO 10 L=1,NL
        DO 10 K=1,NK
          ZUPT(L,K)=0.0
   10 CONTINUE
      NEWES = 0.
      SUMSUB = 0.
      CMXIRR=0.
      IF(RAIN.GT.0.0) THEN
        CUMRAN = CUMRAN + RAIN
        IF((MTHIRR(DAYNUM).EQ.2).OR.(MTHIRR(DAYNUM).EQ.4)) THEN
          IF(ACTRAIN.GT.0.01) THEN
            NOITR = IFIX(AMAX1(1.,ACTRAIN/12.7)+0.5)
            ACTRAIN = ACTRAIN/FLOAT(NOITR)
            DO ITER=1,NOITR
              CALL GRAFLO
              FULPRO=.TRUE.
              CALL CAPFLO
            ENDDO
          ENDIF
          CALL DRIPFLO
          FULPRO=.TRUE.
          NOITR = 1
          CALL CAPFLO
        ELSE
          NOITR = IFIX(AMAX1(1.,RAIN/12.7)+0.5)
          RAIN=RAIN/FLOAT(NOITR)
          ACTRAIN = ACTRAIN/FLOAT(NOITR)
          ACTIRRG = ACTIRRG/FLOAT(NOITR)
          DO ITER=1,NOITR
            CALL GRAFLO
            FULPRO=.TRUE.
            CALL CAPFLO
          ENDDO

! *** Fertigation concentrations need to be reset to zero  
! *** Daily applications are added immedaitely  gtheseira

          CONNIT = 0.
          CONAMM = 0.
          CONURA = 0.
          RAIN = RAIN*FLOAT(NOITR)
        ENDIF
      ENDIF
      CALL ET
      SUPNH4 = 0.
      SUPNO3 = 0.
      SUMEP = 0.
      DO 860 L=1,40
         DO 860 K=1,21
            FNL(L,K) = 0.
            IF(K.LE.20) FNU(L,K) = 0.
            IF(K.LE.11) FWL(L,K) = 0.
            IF(K.LE.10) FWU(L,K) = 0.
 860  CONTINUE
      NOITR = 5
      DO 60 ITER=1,NOITR
        IF(ITER.EQ.1) THEN
          FULPRO=.TRUE.
        ELSE
          FULPRO=.FALSE.
        ENDIF
        CALL EVPSOIL
        IF(DAYNUM.GE.EMERGE) THEN
          CALL UPTAKE
          IF(UPNO3.GT.0.) SUPNO3 = SUPNO3 + UPNO3
          IF(UPNH4.GT.0.) SUPNH4 = SUPNH4 + UPNH4
        ENDIF
        CALL CAPFLO
   60 CONTINUE
      CALL NITRIF
      UPTAKEN=UPTAKEN+((SUPNO3+SUPNH4)/ACELLDW)*vcell/rowsp/0.011219

! *** 0.001 CONVERTS MG OF N TO GRAMS

      SUPNO3 = SUPNO3 * POPFAC * .001
      SUPNH4 = SUPNH4 * POPFAC * .001
      CUMES = CUMES + NEWES
      CUMEP = CUMEP + NEWEP

! *** TOTAL WATER PROFILE

      TH2O = 0.
      TNNO3 = 0.
      TNNH4 = 0.
      DO 149 L=1,NL
        DO 148 K=1,NK
          TNNO3 = TNNO3 + VNO3C(L,K)
          TNNH4 = TNNH4 + VNH4C(L,K)
          TH2O  = TH2O  + VH2OC(L,K)
 148    continue
 
 149  CONTINUE

! *** 0.891*D*W (0.891*ACELLDW) CONVERTS MG N/CM**3 ADDED THE WHOLE SOIL
! *** PROFILE TO LBS N /ACRE. (MG/.01 M**2) * (10000 M**2/2.475 ACRES) *
! *** (0.001 G/MG) * (LBS/454 G) = 0.891 LBS/ACRE

      TNNO3 = TNNO3 * vcell/rowsp/0.011219
      TNNH4 = TNNH4 * vcell/rowsp/0.011219
      SOILN = TNNO3 + TNNH4
	count_of_columns = dfloat(nk)
      TH2O = (TH2O * DCELL *10.)/count_of_columns
      SUBIRR = SUBIRR + (SUMSUB * 10.)/(count_of_columns*WCELL*TCELL)
      RETURN
      END


      SUBROUTINE FRTLIZ
! *****************************************************************
! *    SUBROUTINE ADDS FERTILIZER TO PROFILE. MUST BE CALLED AT   *
! * PLANTING DATE TO INITIALIZE NITROGEN & ORGANIC MATTER         *
! * PROFILE.  MAY BE CALLED FOR SIDE DRESSING.  INPUTS ARE:       *
! * FERN:  FERTILIZER INORGANIC NITROGEN, LBS N/ACRE.             *
! * FNH4:  FRACTION OF INORGANIC N  IN AMMONIA FORM. 0 TO 1       *
! * FNO3:  FRACTION OF INORGANIC N  IN NITRATE FORM. 0 TO 1       *
! * DR:    DISTANCE TO RIGHT OF ROW OF BAND OF FERTILIZER, IN.    *
! *        EQUALS 0 IF BROADCAST.                                 *
! * DD:    DISTANCE BELOW SOIL SURFACE OF BAND OF FERTILIZER,     *
! *        INCHES.  IGNORED IF DR = 0.                            *
! * OMA:   ORGANIC MATTER PLOWED AT BEGINNING OF SEASON, %,       *
! *        MUST BE .GT. 0 TO INITIALIZE N & ORGANIC MATTER ARRAYS.*
! * RNNO3: RESIDUAL N  AS NITRATE IN UPPER 20 CM, LBS/ACRE.       *
! * RNNH4: RESIDUAL N AS AMMONIUM IN UPPER 20 CM, LBS/ACRE.       *
! * 1 LB/ACRE = .0112186 MG/CM**2   OR  .0005609 MG/20CM**3       *
! *****************************************************************

      use common_block
 
! *** OMA .GT. 0. IMPLIES INITIAL FERTILIZATION AT PLANTING DATE &
! *** PLOWDOWN OF ORGANIC MATTER.

      FERN=0.
        I = DAYNUM
        XFERT = NFERT(I,2)+NFERT(I,3)+NFERT(I,4)
        IF(XFERT.GT.0) THEN
           FERN=XFERT+FERN
 
! *** BROADCAST NITROGEN APPLICATION
 
           IF(NFERT(I,5).EQ.0) THEN
              FERTFAC = .011219/DCELL
              NLAYERS = 1
              ADDAMM = NFERT(I,2)*FERTFAC
              ADDNIT = NFERT(I,3)*FERTFAC
              ADDURA = NFERT(I,4)*FERTFAC
              ADDEDN = ADDEDN+((ADDAMM+ADDNIT+ADDURA)*NLAYERS*NK)* &
                                               VCELL/ROWSP/.011219
              DO 100 L=1,NLAYERS
                DO 100 K=1,NK
                  VNO3C(L,K) = VNO3C(L,K)+ADDNIT
                  VNH4C(L,K) = VNH4C(L,K)+ADDAMM
                  VNH4C(L,K) = VNH4C(L,K)+ADDURA
  100         CONTINUE
           ENDIF

! *** FOLIAR APPLICATIONS. FOLIARN IS GRAMS OF FOLIAR N RECEIVED PER PLANT.
! *** 454 CONVERTS LBS TO GRAMS. 70% OF UREA & AMMONIA INTERCEPTED BY CANOPY
! *** IS ADDED TO LEAF N RESERVES. NITRO SUBROUTINE USES FOR PLANT GROWTH.

           IF(NFERT(I,5).EQ.2) THEN
              FERINT = (NFERT(I,2)+NFERT(I,4))*INT
              FERSOI = NFERT(I,2)+NFERT(I,4)-FERINT
              FOLIARN = (FERINT*454.)/POPPLT*.70
              SLEAFN = SLEAFN+FOLIARN
              ADDAMM = FERSOI*.011219/DCELL
              ADDNIT = NFERT(I,3)*.011219/DCELL
              ADDEDN = ADDEDN+((ADDAMM+ADDNIT)*NK)*VCELL/ROWSP/.011219
              DO 120 K=1,NK
                VNO3C(1,K) = VNO3C(1,K)+ADDNIT
                VNH4C(1,K) = VNH4C(1,K)+ADDAMM
  120         CONTINUE
           ENDIF

! ***  SIDEDRESS NITROGEN APPLICATION

           IF(NFERT(I,5).EQ.1) THEN

! *** Relative position of sidedress changes wrt plant position 
! *** under solid or skip row configuration           gtheseira 

              K = (kulclf+3) + (NFERT(I,6) * 2.54 / WCELL + 1.5)
              L = NFERT(I,7) * 2.54 / DCELL + 1.5
              IF(K.LE.0) K=1
              K1=K+1
              IF(K.GE.NK) THEN
                 K=NK
                 K1=NK-1
              ENDIF
              IF(L.LE.0) L=1
              L1=L+1
              IF(L.GE.NL) THEN
                 L=NL
                 L1=NL-1
              ENDIF
              N00=100./ACELLDW+.5
              IF(N00.LT.1) N00=1
              IF(N00.GT.4) N00=4
              FERADD = NFERT(I,3)*.011219*ROWSP/VCELL/N00
              IF(FERADD.GT.0.) THEN
                ADDEDN = ADDEDN + (FERADD*N00)*VCELL/ROWSP/.011219
                VNO3C(L,K)=VNO3C(L,K) + FERADD
                IF(N00.GE.2) VNO3C(L,K1)=VNO3C(L,K1) + FERADD
                IF(N00.GE.3) VNO3C(L1,K)=VNO3C(L1,K) + FERADD
                IF(N00.EQ.4) VNO3C(L1,K1)=VNO3C(L1,K1) + FERADD
              ENDIF
              FERADD = NFERT(I,2)*.011219*ROWSP/VCELL/N00
              IF(FERADD.GT.0.) THEN
                ADDEDN = ADDEDN + (FERADD*N00)*VCELL/ROWSP/.011219
                VNH4C(L,K)=VNH4C(L,K) + FERADD
                IF(N00.GE.2) VNH4C(L,K1)=VNH4C(L,K1) + FERADD
                IF(N00.GE.3) VNH4C(L1,K)=VNH4C(L1,K) + FERADD
                IF(N00.EQ.4) VNH4C(L1,K1)=VNH4C(L1,K1) + FERADD
              ENDIF
              FERADD = NFERT(I,4)*.011219*ROWSP/VCELL/N00
              IF(FERADD.GT.0.) THEN
                ADDEDN = ADDEDN + (FERADD*N00)*VCELL/ROWSP/.011219
                VNH4C(L,K)=VNH4C(L,K) + FERADD
                IF(N00.GE.2) VNH4C(L,K1)=VNH4C(L,K1) + FERADD
                IF(N00.GE.3) VNH4C(L1,K)=VNH4C(L1,K) + FERADD
                IF(N00.EQ.4) VNH4C(L1,K1)=VNH4C(L1,K1) + FERADD
              ENDIF
           ENDIF

! *** BROADCAST NITROGEN APPLICATION WITH INCORPORATION
! *** N fertilizer is immediately available.  Code that locates
! *** pre-emergence N in the first layer and post-emergence in
! *** the plow layer, has been removed. This should be user   
! *** determined--nfert(i,5) = 0.                  gtheseira

           IF(NFERT(I,5).EQ.3) THEN
              FERTFAC = .011219/(DCELL*LPLOW)
              NLAYERS = LPLOW
              ADDAMM = NFERT(I,2)*FERTFAC
              ADDNIT = NFERT(I,3)*FERTFAC
              ADDURA = NFERT(I,4)*FERTFAC
              ADDEDN = ADDEDN+((ADDAMM+ADDNIT+ADDURA)*NLAYERS*NK)* &
                                              VCELL/ROWSP/.011219
              DO L=1,NLAYERS
                DO K=1,NK
                  VNO3C(L,K) = VNO3C(L,K)+ADDNIT
                  VNH4C(L,K) = VNH4C(L,K)+ADDAMM
                  VNH4C(L,K) = VNH4C(L,K)+ADDURA
                ENDDO
              ENDDO
           ENDIF

! ***   NITROGEN FERTILIZER APPLIED IN IRRIGATION 
! *** In FRTLIZ, three variables are calculated representing the 
! *** concentrations of the various N sources im mg/cubic cm of
! *** irrigation water.  The corresponding amounts of N can then 
! *** be added to the appropriate cells in subroutine GRAFLO     

           IF(NFERT(I,5).EQ.4) THEN
              ADDEDN = ADDEDN + NFERT(I,2)+NFERT(I,3)+NFERT(I,4)
              CONAMM = NFERT(I,2) * 0.011219 / (ACTIRRG*.1)
              CONNIT = NFERT(I,3) * 0.011219 / (ACTIRRG*.1)
              CONURA = NFERT(I,4) * 0.011219 / (ACTIRRG*.1)
           ENDIF
        ENDIF
      RETURN
      END

 
      SUBROUTINE GRAFLO
! ************************************************************
! *                                                          *
! *   GRAVITY FLOW OF NO3 AND H2O, AFTER RAIN OR IRRIGATION. *
! *                                                          *
! ************************************************************
! *** RAIN OR IRRIGATION IS IN MM.

      use common_block

      REAL*4 NO3STAY

      DO K=1,NK
        SOAKN(K) = 0.
      END DO

! ************************************************************
! * Calculation of water to be added to the left or right of the
! * plant.  Two variables are introduced (H2OADDL & H2OADDR)
! * which sum to TH2OADD.  Symmetrical irrigation & rain is
! * applied as H2OADDL & H2OADDR; whereas alternate furrow or
! * alternate drip irrigation as applied as H2OADDL.  Because of
! * this, H2OADDL can be only equal to or greater than H2OADDR.
! * Skiprow format is accomplished through convenient use of
! * cultivation index cells to determine initial and terminal
! * column cells in each layer iteration.
! ************************************************************

      IF((MTHIRR(DAYNUM).EQ.2).OR.(MTHIRR(DAYNUM).EQ.4)) THEN

! *** If drip/alternate drip irrigation method  kit

        H2OADDL = ACTRAIN*.1*(KULCLF+3)*WCELL*1
        H2OADDR = ACTRAIN*.1*(NK-(KULCLF+3))*WCELL*1
      ELSEIF(MTHIRR(DAYNUM).EQ.3) THEN

! *** If alternate furrow irrigation method  kit

        H2OADDL = (ACTRAIN+(2*ACTIRRG))*.1*(KULCLF+3)*WCELL*1
        H2OADDR = ACTRAIN*.1*(NK-(KULCLF+3))*WCELL*1
      ELSE

! *** If sprinkler/furrow irrigation method   Kit

        H2OADDL = (ACTRAIN+ACTIRRG)*.1*(KULCLF+3)*WCELL*1
        H2OADDR = (ACTRAIN+ACTIRRG)*.1*(NK-(KULCLF+3))*WCELL*1
      ENDIF

! **** Profile to left of the plant. gtheseira

      L = 1
      DO WHILE ((L.LE.NL).AND.(H2OADDL.GT.0))

        SUMH2O = 0.
        DO K=1,KULCLF+3
          SUMH2O = SUMH2O+VH2OC(L,K)*VCELL
        END DO

        H2ODEF = THTS(L)*VCELL*(KULCLF+3)-SUMH2O
        IF(H2OADDL.GE.H2ODEF) THEN
          H2OMAX = THTS(L)*VCELL
          DO K=1,KULCLF+3
            TOTNO3 = VNO3C(L,K)*VCELL+SOAKN(K)
            TH2OCELL = (H2OADDL-H2ODEF)/(KULCLF+3)+H2OMAX
            PRCENT = H2OMAX/TH2OCELL
            IF(PRCENT.GT.1.) PRCENT = 1.
            NO3STAY = PRCENT*TOTNO3
            SOAKN(K) = TOTNO3-NO3STAY
            VNO3C(L,K) = NO3STAY/VCELL+(THTS(L)-VH2OC(L,K)) &
                                       *CONNIT
            VNH4C(L,K) = VNH4C(L,K)+(THTS(L)-VH2OC(L,K))	&
                                       *CONAMM
            VNH4C(L,K) = VNH4C(L,K)+(THTS(L)-VH2OC(L,K))	&
                                       *CONURA
            VH2OC(L,K) = THTS(L)
          END DO
          H2OADDL = H2OADDL-H2ODEF
        ELSE
          DO K=1,KULCLF+3
            IF(H2OADDL.GT.0) THEN
              H2ODEF = (THTS(L)-VH2OC(L,K))*VCELL
              IF(H2OADDL.LT.H2ODEF) THEN
                VH2OC(L,K) = VH2OC(L,K)+H2OADDL/VCELL
                VNO3C(L,K) = VNO3C(L,K)+(H2OADDL/VCELL*CONNIT)
                VNH4C(L,K) = VNH4C(L,K)+(H2OADDL/VCELL*CONAMM)
                VNH4C(L,K) = VNH4C(L,K)+(H2OADDL/VCELL*CONURA)
                H2OADDL = 0.
              ELSE
                VNO3C(L,K) = VNO3C(L,K)+(THTS(L)-VH2OC(L,K))*CONNIT
                VNH4C(L,K) = VNH4C(L,K)+(THTS(L)-VH2OC(L,K))*CONAMM
                VNH4C(L,K) = VNH4C(L,K)+(THTS(L)-VH2OC(L,K))*CONURA
                VH2OC(L,K) = THTS(L)
                H2OADDL = H2OADDL-H2ODEF
              ENDIF
              VNO3C(L,K) = VNO3C(L,K)+SOAKN(K)/VCELL
              SOAKN(K) = 0.
            ELSE
              IF(L.GT.1) THEN
                VNO3C(L-1,K) = VNO3C(L-1,K)+SOAKN(K)/VCELL
                SOAKN(K) = 0.
              ENDIF
            ENDIF
          END DO
        ENDIF
        IF((L.EQ.NL).AND.(H2OADDL.GT.0.)) THEN
          CUMNSOK=0.
          DO I=1,KULCLF+3
            CUMNSOK=CUMNSOK+SOAKN(I)
          END DO
        ENDIF
        L=L+1
      END DO

! *** Profile to right of the plant.  gtheseira

      L = 1
      DO WHILE ((L.LE.NL).AND.(H2OADDR.GT.0))
        SUMH2O = 0.
        DO K=KULCRT-3,NK
          SUMH2O = SUMH2O+VH2OC(L,K)*VCELL
        END DO
        H2ODEF = THTS(L)*VCELL*(NK-(KULCLF+3))-SUMH2O
        IF(H2OADDR.GE.H2ODEF) THEN
          H2OMAX = THTS(L)*VCELL
          DO K=KULCRT-3,NK
            TOTNO3 = VNO3C(L,K)*VCELL+SOAKN(K)
            TH2OCELL = (H2OADDR-H2ODEF)/(NK-(KULCLF+3))+H2OMAX
            PRCENT = H2OMAX/TH2OCELL
            IF(PRCENT.GT.1.) PRCENT = 1.
            NO3STAY = PRCENT*TOTNO3
            SOAKN(K) = TOTNO3-NO3STAY
            VNO3C(L,K) = NO3STAY/VCELL+(THTS(L)-VH2OC(L,K))	&
                                       *CONNIT
            VNH4C(L,K) = VNH4C(L,K)+(THTS(L)-VH2OC(L,K))	&
                                       *CONAMM
            VNH4C(L,K) = VNH4C(L,K)+(THTS(L)-VH2OC(L,K))	&
                                       *CONURA
            VH2OC(L,K) = THTS(L)
          END DO
          H2OADDR = H2OADDR-H2ODEF
        ELSE
          DO K=NK,KULCRT-3,-1
            IF(H2OADDR.GT.0) THEN
              H2ODEF = (THTS(L)-VH2OC(L,K))*VCELL
              IF(H2OADDR.LT.H2ODEF) THEN
                VH2OC(L,K) = VH2OC(L,K)+H2OADDR/VCELL
                VNO3C(L,K) = VNO3C(L,K)+(H2OADDL/VCELL*CONNIT)
                VNH4C(L,K) = VNH4C(L,K)+(H2OADDL/VCELL*CONAMM)
                VNH4C(L,K) = VNH4C(L,K)+(H2OADDL/VCELL*CONURA)
                H2OADDR = 0.
              ELSE
                VNO3C(L,K) = VNO3C(L,K)+(THTS(L)-VH2OC(L,K))*CONNIT
                VNH4C(L,K) = VNH4C(L,K)+(THTS(L)-VH2OC(L,K))*CONAMM
                VNH4C(L,K) = VNH4C(L,K)+(THTS(L)-VH2OC(L,K))*CONURA
                VH2OC(L,K) = THTS(L)
                H2OADDR = H2OADDR-H2ODEF
              ENDIF
              VNO3C(L,K) = VNO3C(L,K)+SOAKN(K)/VCELL
              SOAKN(K) = 0.
            ELSE
              IF(L.GT.1) THEN
                VNO3C(L-1,K) = VNO3C(L-1,K)+SOAKN(K)/VCELL
                SOAKN(K) = 0.
              ENDIF
            ENDIF
          END DO
        ENDIF
        IF((L.EQ.NL).AND.(H2OADDR.GT.0.)) THEN
          CUMNSOK=0.
          DO I=KULCRT-3,NK
            CUMNSOK=CUMNSOK+SOAKN(I)
          END DO
        ENDIF
        L=L+1
      END DO
      CUMSOK = CUMSOK+(H2OADDL+H2OADDR)/.10/(NK*WCELL*1)
      RETURN
      END


      SUBROUTINE DRIPFLO
!  ************************************************************
!  *                                                          *
!  *    FLOW OF H2O AFTER DRIP IRRIGATION.                    *
!  *                                                          *
!  ************************************************************
!  *** IRRIGATION IS IN MM.

      use common_block

      IF(MTHIRR(DAYNUM).EQ.2) THEN
        KOLUMM = NFERT(DAYNUM,6) * 2.54 / WCELL + 1.5
        IF((KOLUMM.LT.1).OR.(KOLUMM.GT.20)) KOLUMM=1
      ELSE
        KOLUMM = 1
      ENDIF
      LAYERM = NFERT(DAYNUM,7) * 2.54 / DCELL + 1.5
      IF((LAYERM.LT.1).OR.(LAYERM.GT.40)) LAYERM=1
      TH2OADD = ACTIRRG*.1*ROWSP*1
      AMTAMM = NFERT(DAYNUM,2)*.011219*ROWSP
      AMTNIT = NFERT(DAYNUM,3)*.011219*ROWSP
      AMTURA = NFERT(DAYNUM,4)*.011219*ROWSP
      KNT = 0
  100 CONTINUE
        KNT = KNT+1
        IF(KNT.EQ.1) THEN
          L = LAYERM
          K = KOLUMM
          CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
        ELSE
          IF(MTHIRR(DAYNUM).EQ.2) THEN
            IF(KNT.LE.11) THEN
              DO I=1,KNT-1
                KOLUME = KOLUMM+I
                IF(KOLUME.GT.20) KOLUME = KOLUME-20
                KOLUMB = KOLUMM-I
                IF(KOLUMB.LT.1) KOLUMB = 20+KOLUMB
                IF(LAYERE.LE.40) THEN
                  L = LAYERE
                  IF(I.EQ.1) THEN
                    K = KOLUMM
                    CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                  ENDIF
                  IF(I.LT.10) THEN
                    K = KOLUME
                    CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                  ENDIF
                  K = KOLUMB
                  CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                ENDIF
                IF(LAYERB.GT.0) THEN
                  L = LAYERB
                  IF(I.EQ.1) THEN
                    K = KOLUMM
                    CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                  ENDIF
                  IF(I.LT.10) THEN
                    K = KOLUME
                    CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                  ENDIF
                  K = KOLUMB
                  CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                ENDIF
              ENDDO
              L = LAYERM
              K = KOLUMB
              CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
              IF(KNT.LE.10) THEN
                K = KOLUME
                CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
              ENDIF
              IF(KNT.GT.2) THEN
                DO I=1,KNT-2
                  LAYERE = LAYERM+I
                  LAYERB = LAYERM-I
                  IF(LAYERE.LE.40) THEN
                    L = LAYERE
                    K = KOLUMB
                    CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                    IF(KNT.LE.10) THEN
                      K = KOLUME
                      CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                    ENDIF
                  ENDIF
                  IF(LAYERB.GT.0) THEN
                    L = LAYERB
                    K = KOLUMB
                    CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                    IF(KNT.LE.10) THEN
                      K = KOLUME
                      CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                    ENDIF
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ELSE
            DO I=1,KNT-1
              KOLUME = KOLUMM+I
              IF(LAYERE.LE.40) THEN
                L = LAYERE
                IF(I.EQ.1) THEN
                  K = KOLUMM
                  CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                ENDIF
                IF(KOLUME.LT.20) THEN
                  K = KOLUME
                  CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                ENDIF
              ENDIF
              IF(LAYERB.GT.0) THEN
                L = LAYERB
                IF(I.EQ.1) THEN
                  K = KOLUMM
                  CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                ENDIF
                IF(KOLUME.LE.20) THEN
                  K = KOLUME
                  CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                ENDIF
              ENDIF
            ENDDO
            IF(KOLUME.LE.20) THEN
              L = LAYERM
              K = KOLUME
              CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
            ENDIF
            IF(KNT.GT.2) THEN
              DO I=1,KNT-2
                LAYERE = LAYERM+I
                LAYERB = LAYERM-I
                IF(LAYERE.LE.40) THEN
                  L = LAYERE
                  K = KOLUME
                  CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                ENDIF
                IF(LAYERB.GT.0) THEN
                  L = LAYERB
                  K = KOLUME
                  CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDIF
        IF((MTHIRR(DAYNUM).EQ.2).AND.(KNT.GT.11)) THEN
          IF(LAYERE.LE.40) THEN
            L = LAYERE
            K = 0
  120       CONTINUE
              K = K+1
              CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
            IF((TH2OADD.GT.0.).AND.(K.LT.20)) GO TO 120
          ENDIF
          IF(LAYERB.GT.0) THEN
            L = LAYERB
            K = 0
  140       CONTINUE
              K = K+1
              CALL H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
            IF((TH2OADD.GT.0.).AND.(K.LT.20)) GO TO 140
          ENDIF
        ENDIF
        LAYERB = LAYERM-KNT
        LAYERE = LAYERM+KNT
      IF((TH2OADD.GT.0).AND.((LAYERB.GT.0).OR.(LAYERE.LT.41))) &
                                                        GO TO 100
      CUMNSOK = CUMNSOK+AMTAMM+AMTNIT+AMTURA
      CUMSOK = CUMSOK+TH2OADD/(ROWSP*1)*10
      RETURN
      END


      SUBROUTINE H2OMOVE(TH2OADD,AMTAMM,AMTNIT,AMTURA,L,K)
!  ************************************************************
!  *                                                          *
!  *                 H2OMOVE SUBROUTINE                       *
!  *                                                          *
!  ************************************************************
!  *** THIS SUBROUTINE MOVES H2O ONE CELL AT A TIME FOR THE DRIP
!  *** IRRIGATION

      use common_block

      IF(TH2OADD.GT.0) THEN
        H2ODEF = (THTS(L)-VH2OC(L,K))*VCELL
        IF(TH2OADD.GE.H2ODEF) THEN
          VH2OC(L,K) = THTS(L)
          XDUM = H2ODEF/TH2OADD
          TH2OADD = TH2OADD-H2ODEF
          VNO3C(L,K) = VNO3C(L,K)+AMTNIT*XDUM/VCELL
          VNH4C(L,K) = VNH4C(L,K)+AMTAMM*XDUM/VCELL
          VNH4C(L,K) = VNH4C(L,K)+AMTURA*XDUM/VCELL
          AMTNIT = AMTNIT-AMTNIT*XDUM
          AMTAMM = AMTAMM-AMTAMM*XDUM
          AMTURA = AMTURA-AMTURA*XDUM
        ELSE
          VH2OC(L,K) = VH2OC(L,K)+TH2OADD/VCELL
          VNO3C(L,K) = VNO3C(L,K)+AMTNIT/VCELL
          VNH4C(L,K) = VNH4C(L,K)+AMTAMM/VCELL
          VNH4C(L,K) = VNH4C(L,K)+AMTURA/VCELL
          TH2OADD = 0
          AMTNIT = 0
          AMTAMM = 0
          AMTURA = 0
        ENDIF
      ENDIF
      RETURN
      END


      SUBROUTINE ET
!  ************************************************************
!  *                                                          *
!  *           EVAPOTRANSPIRATION SUBROUTINE                  *
!  *                                                          *
!  ************************************************************
!  *** SUBROUTINE TAKEN ALMOST ENTIRELY FROM RITCHIE, A MODEL
!  *** FOR PREDICTING EVAPORATION FROM A ROW CROP WITH INCOMPLETE
!  *** COVER. WATER RESOURCES RESEARCH VOL. 8:1204.

      use common_block
 
      U=6.
      P = RAIN
      RS = RI*.0169491525

! *** RS = SOLAR RADIATION IN MM H2O/DAY.

      TAVM1 = TAVG-1.
      DEL = VP(TAVG) - VP(TAVM1)

! *** DEL=SLOPE OF SATURATION VAPOR PRESSURE CURVE AT MEAN AIR TEMP.

      LAMDA = INT*LAMDAC + (1.-INT)*LAMDAS

! *** LAMDAC & LAMDAS = ALBEDOS OF CROP & SOIL.
! *** INT=INTERCEPTION (FRACTION OF INCIDENT RS)
! *** RNO=NET RADIATION ABOVE CANOPY (MM/DAY)

      RNO=(RS-LAMDA*RS)

! *** TDRY & TWET = DRY AND WET BULB TEMPERATURES.

      TDRY = TAVG
      VPO = VP(TDRY)
      TWET = TMIN

! *** EO=POTENTIAL EVAPORATION RATE ABOVE CANOPY (MM/DAY)
! *** MODIFIED PENMAN EQ.
! *** WIND = WINDSPEED AT 2 METERS (MILES/DAY)
! *** GAMMA=PSYCHROMETER CONSTANT

      VPA = VP(TWET)
      EO=(RNO*DEL/GAMMA+.262*(1.+0.0061*WIND)*(VPO-VPA))/(DEL/GAMMA+1.)

! *** THE FOLLOWING CALCULATES ESO(POTENTIAL EVAP. RATE AT SOIL SURFACE1
! *** RNS=NET RADIATION AT SOIL SURFACE BELOW CANOPY

      RNS=((1.-INT)-(1.-INT)*LAMDAS)*RS
      ESO=DEL*RNS/(DEL+GAMMA)

! *** STAGE I DRYING
! *** SESI=CUMULATIVE STAGE ONE EVAPORATION FROM SOIL SURFACE
! *** U=UPPER LIMIT OF SESI

      IF(SESI.GT.U)GOTO 100

! *** P=RAINFALL

      IF(P.GE.SESI)GOTO 101
      SESI=SESI-P
 99   SESI=SESI+ESO
      IF(SESI.GE.U)GOTO 102
      ES=ESO
      GOTO 110
 102  ES=ESO-.4*(SESI-U)
      SESII=.6*(SESI-U)
      DUMY01 = SESII / ALPHA
      T = DUMY01 * DUMY01
      GO TO 110
 101  SESI=0.
      GO TO 99

! *** STAGE II DRYING

 100  IF(P.GE.SESII)GO TO 103
      T=T+1.
      ES = ALPHA * (SQRT(T)-SQRT(T-1.))
      IF(P.GT.0.)GO TO 104
      IF(ES.GT.ESO)GO TO 105
 106  SESII=SESII+ES-P
      DUMY02 = SESII / ALPHA
      T = DUMY02 * DUMY02
      GO TO 110
 105  ES=ESO
      GO TO 106
 104  ESX=0.8*P
      IF(ESX.LT.ES)GO TO 107
 111  IF(ESX.GT.ESO)GO TO 108
 109  ES=ESX
      GO TO 106
 108  ESX=ESO
      GO TO 109
 107  ESX=ES+P
      GO TO 111
 103  P=P-SESII
      SESI=U-P
      IF(P.GT.U)GO TO 101
      GO TO 99

! *** TRANSPIRATION IS PROPORTIONAL TO LIGHT INTERCEPTION (INT).
! *** THIS REPRESENTS A MODIFICATION TO RITCHIE'S MODEL.

 110  EP = INT * EO
      IF(EP.GT.(EO-ES)) EP=EO-ES

!       E = ES + EP

      AVGPSI = -1. * PSIAVG
      IF(AVGPSI.GT.0.8) AVGPSI = 0.8
      RN = RI*.71536-26.

! *** RFEP = REDUCTION FACTOR FOR EVAPORATION FROM PLANT.  BASED ON
! *** UNPUBLISHED DATA OF BAKER & HESKETH. 1969.

      RFEPN = 749.5831405 + 0.9659065*RN - 54.6600986*TAVG              &
       - 194.6508431*AVGPSI - 0.0010226*RN*RN + 1.0153007*TAVG*TAVG +	&
       29.775978*AVGPSI*AVGPSI + 0.0293687*RN*TAVG						&
       - 4.206856*TAVG*AVGPSI
      RFEPD = 749.5831405 + 0.9659065*RN					&
       - 54.6600986*TAVG - 19.46508431 - 0.0010226*RN*RN +	&
       1.0153007*TAVG*TAVG + .29775978 + 0.0293687*RN*TAVG	&
       - .4206856*TAVG
      RFEP = RFEPN/RFEPD*CALBRT(25)
      IF(RFEP.LT.(0.2*LAI/2.0)) RFEP = 0.2 * LAI/2.0
      IF(RFEP.LT.0.2) RFEP = 0.2
      IF(RFEP.GT.1.0) RFEP = 1.0
      EP = EP * RFEP
      RETURN
      END 

 
      FUNCTION VP(TMP)
 
! *** THIS FUNCTION CALCULATES SATURATION VAPOR PRESSURE
! *** FOR AIR TEMPERATURE.
 
      VP = EXP(1.8282+TMP*(0.07046136-TMP*0.000215743))
      RETURN
      END


      SUBROUTINE EVPSOIL
! ***********************************************************
! *    EVAPORATION OF WATER FROM THE SOIL IS APPROXIMATED   *
! * FROM THE TOTAL VOLUMETRIC WATER CONTENT OF THE EXPOSED  *
! * CELL COLUMNS CALCULATED AS A FRACTION OF INTERCEPTION.  *
! ***********************************************************
! ***  UNITS     ES : MM/DAY   H2OREML/R : CM**3/CM**3   

      use common_block

      NCELLSL = 0
      NCELLSR = 0
      SUMESL = 0.
      SUMESR = 0.
      ESPART = ES/NOITR

! *** Determine parameters for cells to left of plant.  gtheseira

      NKESL=MAX0(1,IFIX((1.0-INT)*(KULCLF+3)))
      H2OREML=(0.1*ESPART*WCELL*TCELL/VCELL)*(KULCLF+3)/NKESL

! *** Determine parameters for cells to right of plant.	 gtheseira

      NKESR=MAX0(1,IFIX((1.0-INT)*(KULCLF+3))+(NK-2*(KULCLF+3)))      
      H2OREMR=(0.1*ESPART*WCELL*TCELL/VCELL)*(NK-(KULCLF+3))/NKESR
      IF (NKESL.GE.1) THEN
        DO K=1,NKESL
          IF (H2OREML.GT.VH2OC(1,K) - AIRDR(1)) THEN
            SUMESL= SUMESL + (VH2OC(1,K) - AIRDR(1))
            VH2OC(1,K) = AIRDR(1)
          ELSE
            SUMESL = SUMESL + H2OREML
            VH2OC(1,K) = VH2OC(1,K) - H2OREML
          ENDIF
          NCELLSL = NCELLSL + 1
        ENDDO
      ENDIF
      IF (NKESR.GE.1) THEN
        DO K=NK+1-NKESR,NK
          IF (H2OREMR.GT.VH2OC(1,K) - AIRDR(1)) THEN
            SUMESR = SUMESR + (VH2OC(1,K) - AIRDR(1))
            VH2OC(1,K) = AIRDR(1)
          ELSE
            SUMESR = SUMESR + H2OREMR
            VH2OC(1,K) = VH2OC(1,K) - H2OREMR
          ENDIF
          NCELLSR = NCELLSR + 1
        ENDDO
      ENDIF
      SUMEST = SUMESL + SUMESR
      NCELLST = NCELLSL + NCELLSR
      NKEST = NKESL + NKESR
      ADJES = (10.*SUMEST*VCELL/(NCELLST*WCELL*TCELL))*NKEST/NK
      NEWES = NEWES + ADJES
      RETURN
      END


      SUBROUTINE UPTAKE
! ************************************************************
! *    UPTAKE OF WATER FROM EACH SOIL CELL IS PROPORTIONAL   *
! * TO THE PRODUCT OF ROOT WEIGHT CAPABLE OF UPTAKE AND THE  *
! * HYDRAULIC CONDUCTIVITY OF THE CELL.  THE SUM OF THE      *
! * UPTAKE FROM THE CELLS EQUALS TRANSPIRATION.  ALL NO3 IN  *
! * THE WATER TAKEN UP BY THE ROOTS IS ALSO TAKEN UP,AND     *
! * ACTIVE UPTAKE OF NO3 IS CALCULATED USING MICHELIS-MENTON *
! * CONSTANTS FOR ROOTS OF 3 DIFFERENT AGE CLASSES.          *
! ************************************************************
! *** EP - TRANSPIRATION BY PLANTS, MM/DAY.
! *** SUPNO3 - SUPPLY OF NITRATE FROM SOIL, MG.
! *** SUPNH4 - SUPPLY OF AMMONIA FROM SOIL, MG.
! *** UPF - UPTAKE FACTOR, GM CM/DAY.
! *** ROOT WEIGHT CAPABLE OF UPTAKE, GM/CELL.

      use common_block
 
! *** SINCE UPTAKE IS CALLED TWO TIMES PER DAY DELT IS 0.5 OR (1.0/2.0)
 
      DELT = 1./NOITR
      DUMY01 = (.10*NK*WCELL*EP)*DELT
      Write(99,1000) DAYNUM, DUMY01, EP, EO, ES
      DO 8 I=1,40
        DO 8 J=1,20
          RTWTCU(I,J) = 0.
          UPF(I,J) = 0.
    8 CONTINUE
      DO 888 I=1,41
        DO 888 J=1,21
          PUPF(I,J)  = 0.
          TUPF(I,J)  = .TRUE.
          TTUPF(I,J) = .TRUE.
  888 CONTINUE
      DO 1 L=1, LR

! *** Determination of initial rooted cell; execution of loop  
! *** from initial to terminal rooted cell.      gtheseira

         KL = KLL(L)
         KR = KRL(L)
         DO 1 K =KL, KR
         RTWTCU(L,K) = (RTWT(L,K,1) + 0.20*(RTWT(L,K,2) + RTWT(L,K,3)))

! *** SUMS THE WEIGHT OF ROOTS 15 DAYS OLD OR LESS IN CELL.
 
      IF(L+2.LE.NL) THEN 
        IF((RTWTCU(L,1).GT.0.).AND.(.NOT.RTEXNT(L+2)))RTEXNT(L+2)=.TRUE.
      ENDIF
 1    CONTINUE

! *** Deletion of entire loop that inserts roots from right half  
! *** of the plant which have grown past NKH.
! *** Definition of NKK not required         gtheseira                  

      SUPF = 0.
      DO 5 L = 1, LR

! *** Determination of initial rooted cell; execution of loop  
! *** from initial to terminal rooted cell.  gtheseira                  

        KL = KLL(L)
        KR = KRL(L)

! *** KR not limited to NKH   

        DO 5 K = KL, KR
          UPF(L,K)=RTWTCU(L,K)*DIFF(L,K)
 
! *** UPTAKE FACTOR FOR EACH CELL, HAS UNITS OF GM CM/DAY.
 
          SUPF = SUPF + UPF(L,K)
 
! *** SUM OF UPTAKE FACTORS IN THE PROFILE.  USED FOR APPORTIONING
! *** UPTAKE AMONG CELLS.
 
 5    CONTINUE
      D01 = 1./SUPF
      DO 556 L=1,LR

! *** Determination of initial rooted cell; execution of loop   
! *** from initial to terminal rooted cell.        gtheseira                   

        KL = KLL(L)
        KR = KRL(L)

! *** KR not limited to NKH       

        DO 555 K = KL, KR
          PUPF(L,K)=UPF(L,K)*D01
  555   CONTINUE
  556 CONTINUE
      DO 656 L=1,LR

! *** Determination of initial rooted cell; execution of loop  
! *** from initial to terminal rooted cell.        gtheseira              

        KL = KLL(L)
        KR = KRL(L)

! *** KR not limited to NKH     

        DO 655 K=KL, KR											  
          IF((PUPF(L,K).LE.0.005).AND.(PUPF(L+1,K).LE.0.005))	   &
              TUPF(L,K) = .FALSE.
          IF((PUPF(L,K).LE.0.005).AND.(PUPF(L,K+1).LE.0.005).AND.  &
             (.NOT.TUPF(L,K))) TTUPF(L,K) = .FALSE.

! *** TTUPF and TUPF calculations for image cells deleted. 

  655   CONTINUE
  656 CONTINUE
      UPNH4 = 0.
      UPNO3 = 0.
      DO 6 L  = 1, LR
      J = LYRDPH(L)

! *** Determination of initial rooted cell; execution of loop  
! *** from initial to terminal rooted cell.         
      KL = KLL(L)
      KR = KRL(L)

! *** KR not limited to NKH       
! *** DUMY01 no longer halved for UPTH2O calculation   

      DO 6 K = KL, KR
      UPTH2O = (UPF(L,K)/SUPF) * DUMY01
      H2OUPT = UPTH2O/ACELLDW
 
! *** UPTAKE OF WATER FROM EACH CELL, CM**3/DAY. EP HAS UNITS OF MM/DAY.
! *** Image column definition not required.                      
 
      H2OUPT = AMAX1(0., AMIN1(VH2OC(L,K)-THETA0(J), H2OUPT))
 
! *** VOLUMETRIC WATER CONTENT OF CELL IS DECREASED BY AMOUNT
! *** OF UPTAKE FROM CELL. VOLUMETRIC WATER CONTENT OF IMAGE CELL IS
! *** ALSO REDUCED.
 
! *** Image cell reduction not required; H2OUPT not doubled.  

      SUMEP = SUMEP + H2OUPT
      UPTH2O = H2OUPT*ACELLDW
      IF(UPTH2O.LE.0.0)GO TO 6
 
! *** CALCULATE AMOUNT OF ACTIVE NO3 AN NH4 UPTAKE IN MG N/CM**3
 
      UPNO3C = UPTH2O*(VNO3C(L,K)/VH2OC(L,K))
      UPNH4C = UPTH2O*((0.10*VNH4C(L,K))/VH2OC(L,K))
 
! *** CALCULATE PASSIVE UPTAKE OF NO3 FROM CELL, MG N/DAY.
 
      AUNO3C = 0.0
      EFOW = (VH2OC(L,K)-THETA0(J))/(FCININ(J)-THETA0(J))
      IF(EFOW.LT.0.) EFOW=0.
      IF(VNO3C(L,K).GT.0.000001) THEN
         MMUPN1 = 9.00/(1.+0.000000248/(VNO3C(L,K)/VH2OC(L,K)))
         MMUPN2 = 9.00/(1.+0.000000248/(VNO3C(L,K)/VH2OC(L,K)))
         MMUPN3 = .900/(1.+0.000000248/(VNO3C(L,K)/VH2OC(L,K)))
         AUNO3C = ((MMUPN1*RTWT(L,K,1)) + (MMUPN2*RTWT(L,K,2)) +  &
                   (MMUPN3*RTWT(L,K,3))) * EFOW * DELT
         IF(AUNO3C.LT.0.0) AUNO3C = 0.0
      ENDIF

! *** Duplicate code for image column cells not required    

! *** NOW UPDATE WATER CONTENT OF CELLS AFTER NITROGEN MOVEMENT 
! AND UPTAKE HAVE BEEN ACCOUNTED FOR.
! *** Image column code not required.                        

      VH2OC(L,K) = VH2OC(L,K) - H2OUPT

      ZUPT(L,K) = ZUPT(L,K) + H2OUPT

! *** REDUCE VOLUMETRIC NO3 CONTENT OF CELLS, BY THE AMOUNT 
! *** OF PASSIVE AND ACTIVE N UPTAKE. MG/CM**3
 

      IF((VNO3C(L,K) - UPNO3C/ACELLDW).LT.0.0) THEN
        UPNO3C = (VNO3C(L,K) - FLNMIN/VCELL) * ACELLDW
        VNO3C(L,K) = FLNMIN/VCELL
        AUNO3C = 0.0
      ELSE
        VNO3C(L,K) = VNO3C(L,K) - UPNO3C/ACELLDW
        IF((VNO3C(L,K) - AUNO3C/ACELLDW) .LT.0.0) THEN
          AUNO3C = (VNO3C(L,K) - FLNMIN/VCELL) * ACELLDW
          VNO3C(L,K) = FLNMIN/VCELL
        ELSE
          VNO3C(L,K) = VNO3C(L,K) - AUNO3C/ACELLDW
        ENDIF
      ENDIF
 
! *** REDUCE VOLUMETRIC NH4 CONTENT OF CELLS, BY THE AMOUNT 
! *** OF PASSIVE N UPTAKE. MG/CM**3

      IF(VNH4C(L,K) - UPNH4C/ACELLDW .LT.0.0) THEN
        UPNH4C = (VNH4C(L,K) - FLNMIN/VCELL) * ACELLDW
        VNH4C(L,K) = FLNMIN/VCELL
      ELSE
        VNH4C(L,K) = VNH4C(L,K) - UPNH4C/ACELLDW
      ENDIF

! *** ACCUMULATE PASSIVE UPTAKE OF NO3 AND NH4 AND THE ACTIVE 
! *** UPTAKE OF NH4 FROM EACH CELL, UNITS AREA MG N/CM**3.
 
      UPNO3 = UPNO3 + UPNO3C + AUNO3C
      UPNH4 = UPNH4 + UPNH4C
 
! *** Duplicate code for image column cells not required  
 
 6    CONTINUE
      IF(ITER.EQ.(NINT(1./DELT))) NEWEP = (SUMEP*DCELL*10.)/NK
1000  format(I6, 2x,F7.5, 2x,F7.5, 2x, F7.5,2x,F7.5)
      RETURN
      END


      SUBROUTINE CAPFLO
! ************************************************************
! *                                                          *
! *  CAPILLARY FLOW OF NO3 AND H2O IN ALL DIRECTIONS.        *
! *                                                          *
! ************************************************************
! *** AFTER ITERATIONS DURING THE DAY HAVE BEEN ADDED.     
! *** FLUX OF H2O TO THE LEFT OUT OF THE CELL, CM**3/CELL/DAY.
! *** FLUX OF H2O UPWARD OUT OF THE CELL, CM**3/CELL/DAY.    
! *** FLUX OF NITROGEN TO THE LEFT OUT OF THE CELL, MG N/CELL/DAY.
! *** FLUX OF NITROGEN UPWARD OUT OF THE CELL, MG N/CELL/DAY. 

       use common_block

! *** SKPFLG CONTROLS THE NUMBER OF ITERATIONS OF CAPFLO THAT 
! *** ARE NECESSARY PER DAY.  IF THERE IS A SUBSTANTIAL AMOUNT  
! *** OF WATER MOVEMENT INTO OR OUT OF ANY CELL, SKPFLG WILL
! *** BE SET TO TRUE.IF IT REMAINS FALSE, THEN CAPFLO WILL NOT
! *** BE CALLED AGAIN IN SOIL.

      SKPFLG = .TRUE.

! *** CALCULATE NECESSARY INDICES

      NLM1 = NL-1
      NKH = NK/2
      NKMM1 = NK-1
      NKHP1 = NKH + 1
      NKHP2 = NKH + 2

! *** CALCULATE FRACTION OF DAY THAT WILL BE PROCESSED DURING 
! *** THIS ITERATION.

      DELT = 1.0 / NOITR

! *** CALCULATE DIFFUSIVITY OF EACH SOIL CELL IN LEFT HALF OF .
! *** PROFILE. ASSUME NO FLOW ACROSS VERTICAL BOUNDARY DUE TO 
! *** SYMMETRY, UNDER ROW AND MIDWAY BETWEEN ROWS. ASSUME NO 
! *** FLUX INTO ROOT ZONE FROM BENEATH. DIFFUSIVITY FUNCTION 
! *** FOUND IN GARDNER AND MAYHUGH. 1966. SSSAP 22:197-201. FDW.
 
      DO 4 L = 1, NL
         J = LYRDPH(L)

! *** Diffusivity has to be calculated for whole slab.

         DO 4 K = 1, NK

! *** CALCULATE AMOUNT OF WATER THAT CAN MOVE. THEN	IF
! *** AMOUNT NEGLIGIBLE, SET DIFFUSIVITY TO DIFF0 FOR LAYER
! *** ELSE CALCULATE IT.

              DUMY02 = VH2OC(L,K)-THETA0(J)
              IF(DUMY02.LT.0.00001)THEN
                   DIFF(L,K)=DIFF0(J)
              ELSE
                   DIFF(L,K) = DIFF0(J)*EXP(BETA(J)*DUMY02)
              ENDIF
    4 CONTINUE

! *** With plant in slab middle, cultivated cells are 7 through 14.
! *** DIFF for these cells redefined as initial DIFF0 value.

      DO 150 K=kulcrt,kulclf
         DUMY02 = VH2OC(1,K)-THTA0I
         IF(DUMY02.LT.0.00001)THEN
              DIFF(1,K)=DIFF0I
         ELSE
              DIFF(1,K) = DIFF0I*EXP(BETAI*DUMY02)
         ENDIF
 150  CONTINUE

! *** CALCULATE POTENTIAL FLUX OF WATER LEFT AND UP IN LEFT HALF OF
! *** SOIL PROFILE.  CALCULATE FLUX OF NITROGEN LEFT AND UP IN BOTH
! *** HALVES OF PROFILE

      DO 5 L  = 1, NL

! *** IF NO ROOTS IN PREVIOUS LAYER AND NOT TIME TO CONSIDER
! *** THE FULL PROFILE, EXIT THIS LOOP

         IF((.NOT.RTEXNT(L)).AND.(.NOT.FULPRO)) GO TO 4406
         J = LYRDPH(L)

! *** Calculate H2O & nitrogen flux in the whole slab.

         DO K = 1, NK

! *** SET WATER TABLE. START WITH THE BOTTOM LAYER AND PROCESS 
! *** THE WATER TABLE FROM THE BOTTOM UPWARD UNTIL A ABOVE THE 
! *** WATER TABLE.

            IT0 = IDAY + 1
            IF(FBLOOM.GT.0.) IT0 = FBLOOM + TD
            IF((L.LT.NL).AND.(WATTBL.LT.DCELL*(L+1))) THEN
               SAVEH2O = VH2OC(L+1,K)
               VH2OC(L+1,K) = THTS(L+1)
               SUMSUB1 = SUMSUB + (satmc - SAVEH2O)*VCELL
               SUMSUB = SUMSUB + (VH2OC(L+1,K) - SAVEH2O)*VCELL
            ENDIF

! *** IF THERE WAS UPTAKE FROM THE CELL ON THE PREVIOUS ITERATION
! *** OR IF THIS IS AN ITERATION OVER ENTIRE PROFILE

            IF(TTUPF(L,K).OR.FULPRO) THEN
               IF(K .NE. 1) THEN

! *** USE HARMONIC MEAN OF DIFFUSIVITY TO	CALCULATE H2O FLOW
! *** K = (2*D1*D2)/(D1+D2)
! *** FWL = FLUX OF WATER TO THE LEFT. CM**3

                  IF(DIFF(L,K)+DIFF(L,K-1).GT.0.00001)THEN
                     FWL(L,K) = (DIFF(L,K)+DIFF(L,K-1))/2.			&
                                 *((VH2OC(L,K)-VH2OC(L,K-1))/WCELL) &
                                 * (DCELL*TCELL*DELT)
                  ELSE
                     FWL(L,K)=0.0
                  ENDIF

! *** H2OMAX IS THE MAXIMUM AMOUNT OF WATER THAT A CELL CAN
! *** RECEIVE FROM ALL DIRECTIONS (LEFT,RIGHT,DOWN,UP) BY 
! *** CAPILLARY ACTION PER ITERATION.  THTS(L) IS THETAS BY 
! *** LAYER

                  WOUT = VH2OC(L,K) - VH2OC(L,K-1)

! *** CALCULATE DIMENSIONLESS FLOW (CC/CC)

                  FLOW = FWL(L,K)/VCELL

! *** THE FOLLOWING CODE WAS INSERTED BY JIM SIEFKER TO PREVENT 
! *** MORE WATER FROM FLOWING INTO A CELL THAN THE CELL CAN HOLD.
! *** IT WILL NOT ALLOW MORE THAN 25% OF THE CAPACITY OF THE CELL
! *** TO COME IN FROM ANY ONE DIRECTION.  THE RATIO TERM PREVENTS 
! *** INSTABILITIES THAT ARISE WHEN TOO MUCH WATER IS MOVED ON
! *** ON ONE ITERATION.
             
                  IF(WOUT.LT.0.0) THEN

! *** CALCULATE THE FLOW OF WATER TO THE RIGHT

                     RATIO = VH2OC(L,K)/VH2OC(L,K-1)
                     H2OMAX = (THTS(L) - VH2OC(L,K)) * (1.0-RATIO) 
                     IF(ABS(FLOW).GT.(0.25 * H2OMAX)) THEN
                        FWL(L,K) = -0.25*H2OMAX*VCELL
                     ENDIF
                  ELSE 

! *** CALCULATE THE FLOW OF WATER TO THE LEFT

                     RATIO = VH2OC(L,K-1)/VH2OC(L,K)
                     H2OMAX = (THTS(L) - VH2OC(L,K-1)) * (1.0-RATIO)
                     IF(ABS(FLOW).GT.(0.25*H2OMAX)) THEN
                        FWL(L,K) = 0.25*H2OMAX*VCELL
                     ENDIF
                  ENDIF

! *** IF MORE THAN MIN WATER FLUX SET SKPFLG TO TRUE

                  IF(ABS(FWL(L,K)).GT.FLXMIN(J)) SKPFLG = .TRUE.

! *** FOLLOWING CODE PREVENTS FLOATING POINT UNDERFLOW
 
                  IF(ABS(FWL(L,K)).LE.FLXMIN(J)) THEN
                     FWL(L,K) = 0.0
                     FNL(L,K) = 0.0
                  ELSE
                     IF(ABS(FWL(L,K)).GT.(ABS(FLXMAX(J))))THEN
                        IF(FWL(L,K).LT.0.)FWL(L,K) = -FLXMAX(J)
                        IF(FWL(L,K).GT.0.)FWL(L,K) =  FLXMAX(J)
                     ENDIF

! *** FLOW OF NO3 TO THE LEFT, OUT OF CELL, IN MG N/CELL/DAY. 
! *** FLOW OF NITROGEN IS PROPORTIONAL TO FLOW OF WATER
                              
                     IF(FWL(L,K) .GT. 0.0) THEN
                        FNL(L,K) = FWL(L,K) * (VNO3C(L,K)/VH2OC(L,K))
                     ELSE
                        FNL(L,K) = FWL(L,K) * (VNO3C(L,K-1)/VH2OC(L,K-1))
                     ENDIF
                     IF(ABS(FNL(L,K)).LT.FLNMIN) FNL(L,K)=0.0

! *** Entire section removed dealing with calculation of FNL on the
! *** NKK side of the slab.

                  ENDIF
               ENDIF

! *** CALCULATE FLUX OF WATER UPWARD INTO EACH CELL FROM NEXT LAYER OF 
! *** CELLS DOWNWARD.	FWU = FLUX OF WATER UPWARD. CM**3             

               IF(L .LT. NL) THEN
                  IF(DIFF(L+1,K)+DIFF(L,K).GT.0.00001) THEN
                     IF(LYRDPH(L+1).EQ.LYRDPH(L)) THEN
                        VH2OCEQ=VH2OC(L,K)
                        FWU(L+1,K)=(DIFF(L+1,K)+DIFF(L,K))/2.		   &
                                   *((VH2OC(L+1,K)-VH2OCEQ)/DCELL)	   &
                                   *(WCELL*TCELL*DELT)
                     ELSE
                        VH2OCEQ=THTR(L+1)+((FC(L+1)-THTR(L+1))*		   &
                              ((PSIS(L,K)/PSISFC)**(1.0/ARDRCN(L+1))))
                        FWU(L+1,K)=(DIFF(L+1,K)+DIFF(L,K))/2.		   &
                                   *((VH2OC(L+1,K)-VH2OCEQ)/DCELL)	   &
                                   *(WCELL*TCELL*DELT)
                     ENDIF
                  ELSE
                     VH2OCEQ=VH2OC(L,K)
                     FWU(L+1,K)=0.0
                  ENDIF     

! *** H2OMAX*.25 IS THE MAXIMUM AMOUNT OF WATER THAT A CELL
! *** CAN RECEIVE FROM EITHER DIRECTION (LEFT,RIGHT,DOWN,UP)
! *** BY CAPILLARY ACTION PERC ITERATION.  THTS(L) IS THETAS 
! *** BY LAYER

                  WOUT = VH2OC(L+1,K) - VH2OCEQ

! *** DIVIDE BY CELL AREA TO CONVERT TO DIMENSIONLESS FLOW (CC/CC)
                    
                  FLOW = FWU(L+1,K)/VCELL
                  IF(WOUT.LT.0.0) THEN

! *** IF WATER MOVEMENT IS DOWNWARD OUT OF THE CELL

                     RATIO = VH2OC(L+1,K)/VH2OCEQ
                     H2OMAX = (THTS(L+1) - VH2OC(L+1,K))*(1.0 - RATIO)
                     IF(ABS(FLOW).GT.(0.25 * H2OMAX)) THEN
                        FWU(L+1,K) = -0.25*H2OMAX*VCELL
                     ENDIF
                  ELSE 

! *** IF WATER MOVEMENT IS UPWARD INTO THE CELL

                     RATIO = VH2OCEQ/VH2OC(L+1,K)
                     H2OMAX = (THTS(L) - VH2OC(L,K))*(1.0 - RATIO) 
                     IF(ABS(FLOW).GT.(0.25 * H2OMAX)) THEN
                        FWU(L+1,K) = 0.25*H2OMAX*VCELL
                     ENDIF
                  ENDIF
                  IF(ABS(FWU(L+1,K)).LE.FLXMIN(J)) THEN
                     FWU(L+1,K) = 0.0
                     FNU(L+1,K) = 0.0
                  ELSE
                     SKPFLG = .TRUE.
                     IF(ABS(FWU(L+1,K)).GT.FLXMAX(J)) THEN
                        IF(FWU(L+1,K).LT.0.) FWU(L+1,K)= -FLXMAX(J)
                        IF(FWU(L+1,K).GT.0.) FWU(L+1,K)=  FLXMAX(J)
                     ENDIF

! *** CALCULATE FLOW OF NITROGEN UPWARD, BOTH PLANES
! *** FLOW OF NO3 UPWARD IN THE WATER, MG N/CELL/DAY.

                     IF(FWU(L+1,K) .GT. 0.0) THEN
                        FNU(L+1,K) = FWU(L+1,K) * (VNO3C(L+1,K)/VH2OC(L+1,K))
                     ELSE
                        FNU(L+1,K) = FWU(L+1,K)	* (VNO3C(L,K)/VH2OC(L,K))
                     ENDIF
                     IF(ABS(FNU(L+1,K)).LT.FLNMIN) FNU(L+1,K) = 0.0
                   ENDIF
                ELSE
                   IF(WATTBL.GT.NL*DCELL) THEN
                      H2OCNT=(IPCLAY(J)-5.)/70.*(FCININ(J)-THETA0(J))+THETA0(J)
                      IF(H2OCNT.GT.FCININ(J)) H2OCNT=FCININ(J)
                      IF(H2OCNT.LT.THETA0(J)) H2OCNT=THETA0(J)
                   ELSE
                      H2OCNT=THTS(NL)
                   ENDIF
                   IF(DIFF(L,K).GT.0.000005) THEN
                      FWU(L+1,K)=DIFF(L,K)*(H2OCNT-VH2OC(L,K))/DCELL*(WCELL*TCELL*DELT)
                   ELSE
                      FWU(L+1,K)=0.0
                   ENDIF

! *** H2OMAX*0.25 IS THE MAXIMUM AMOUNT OF WATER THAT A CELL CAN
! *** RECEIVE FROM EITHER DIRECTION (LEFT,RIGHT,DOWN,UP) BY
! *** CAPILLARY ACTION PER ITERATION. THTS(L) IS THETAS BY LAYER

                   WOUT = H2OCNT - VH2OC(L,K)

! *** DIVIDE BY CELL AREA TO CONVERT TO DIMENSIONLESS FLOW, (CC/CC)

                   FLOW = FWU(L+1,K)/VCELL
                   IF(WOUT.LT.0.0) THEN

! *** IF WATER MOVEMENT IS DOWNWARD OUT OF A CELL

                      RATIO = H2OCNT / VH2OC(L,K)
                      H2OMAX = (THTS(L) - H2OCNT)*(1.0 - RATIO)
                      IF(ABS(FLOW).GT.(0.25 * H2OMAX)) THEN
                         FWU(L+1,K) = -0.25*H2OMAX*VCELL
                      ENDIF
                   ELSEIF(NL*DCELL.GE.WATTBL) THEN

! *** FOR WATER TO MOVE UPWARD THERE MUST BE A H2O TABLE

                      RATIO = VH2OC(L,K) / H2OCNT
                      H2OMAX = (THTS(L) - VH2OC(L,K))*(1.0 - RATIO) 
                      IF(ABS(FLOW).GT.(0.25 * H2OMAX)) THEN
                         FWU(L+1,K) = 0.25*H2OMAX*VCELL
                      ENDIF
                   ELSE
                      FWU(L+1,K) = 0.
                   ENDIF

! *** SUMSUB - Units are cm**3/(nk cells)

                   SUMSUB = SUMSUB+FWU(L+1,K)
                   IF(ABS(FWU(L+1,K)).LE.FLXMIN(J)) THEN
                      FWU(L+1,K) = 0.0
                      FNU(L+1,K) = 0.0
                   ELSE
                      SKPFLG = .TRUE.
                      IF(ABS(FWU(L+1,K)).GT.FLXMAX(J)) THEN
                         IF(FWU(L+1,K).LT.0.) FWU(L+1,K)= -FLXMAX(J)
                         IF(FWU(L+1,K).GT.0.) FWU(L+1,K)=  FLXMAX(J)
                      ENDIF

! *** CALCULATE FLOW OF NITROGEN UPWARD, BOTH PLANES
! *** FLOW OF NO3 UPWARD IN THE WATER, MG N/CELL/DAY.

                      FNU(L+1,K) = FWU(L+1,K) * (VNO3C(L,K)/H2OCNT)
                      IF((ABS(FNU(L+1,K)).LT.FLNMIN).OR.  &
                         (FNU(L+1,K).GT.0.)) FNU(L+1,K) = 0.0
                      CUMNSOK=CUMNSOK+FNU(L+1,K)
                   ENDIF
                ENDIF
            ENDIF
	   enddo
 5    CONTINUE
 4406 CONTINUE

! *** CALCULATE TOTAL FLUX OF WATER NITROGEN FOR EACH CELL
! *** FWICN = FLUX OF H2O INTO THE CELL, NET, CM**2. 
! *** FWICN/ACELLDW= CM**2/CM**2 (DIMENSIONLESS)
! *** FNICN IS FLUX OF NO3 INTO THE CELL, NET, MG N/CELL/DAY.
! *** VNO3C IS VOLUMETRIC CONTENT OF SOIL CELL, MG N/CM**3.

      DO 16 L = 1,NL
         IF((.NOT.RTEXNT(L)).AND.(.NOT.FULPRO)) GO TO 4706
         J = LYRDPH(L)

! *** Calculations run for entire slab.  gtheseira

         DO K = 1,NK

            IF(TTUPF(L,K).OR.FULPRO) THEN
               FWICN = FWL(L,K+1) - FWL(L,K) + FWU(L+1,K) - FWU(L,K)
               VH2OC(L,K) = VH2OC(L,K) + FWICN/VCELL
               IF(VH2OC(L,K).GT.THTS(L)) VH2OC(L,K) =THTS(L)
               IF(VH2OC(L,K).LT.AIRDR(J)) VH2OC(L,K)=AIRDR(J)

! *** CALCULATE N FLUX FOR EACH CELL.  EXTREME COLUMN IN EACH
! *** PLANE MUST BE TREATED SEPARATELY.              

               IF(K .EQ. 1) THEN

! *** CALCULATE FLUX FOR FAR LEFT CELL AND FAR RIGHT CELL

                  FNICN = FNL(L,2) + FNU(L+1,1) - FNU(L,1)
                  VNO3C(L,1) = AMAX1(VNO3C(L,1) + FNICN/VCELL, &
                                     FLNMIN/VCELL)
                  FNICN = FNL(L,NK-1) + FNU(L+1,NK) - FNU(L,NK)
                  VNO3C(L,NK) = AMAX1(VNO3C(L,NK) + FNICN/VCELL, &
                                      FLNMIN/VCELL)
               ELSE
                  FNICN = FNL(L,K+1) - FNL(L,K) + FNU(L+1,K) - FNU(L,K)
                  VNO3C(L,K) = AMAX1(VNO3C(L,K) + FNICN/VCELL,	&
                                     FLNMIN/VCELL)
               ENDIF
            ENDIF
	   enddo
  16  CONTINUE
 4706 CONTINUE

! *** UPDATE THE VALUE OF PSIS FOR EACH CELL.  UPDATE NOT NECESSARY ON
! *** 3RD OR 4TH ITERATION BECAUSE THERE WILL NOT BE AN INTERVENING
! *** CALL TO UPTAKE ON THESE CALLS. IF SKPFLG IS STILL FALSE AT THIS
! *** POINT, THERE WAS NO SIGNIFICANT WATER MOVEMENT SO THIS MAY BE 
! *** THE LAST ITERATION DONE AND SO VALUES MUST BE UPDATED.

      DO 8 L = 1,NL
         J = LYRDPH(L)

! ***DO for full profile width    gtheseira

         DO K = 1, NK
            PSIS(L,K) = -15.
            IF(VH2OC(L,K).GT.THETAR(J)) THEN
               TEMP1C = (VH2OC(L,K)-AIRDR(J))/(FCININ(J)-AIRDR(J))
               PSIS(L,K) = PSISFC * TEMP1C**ARDRCN(L)
               IF(PSIS(L,K).GT.-0.00001)PSIS(L,K)=-0.00001
            ENDIF
            IF(PSIS(L,K).LT.-15.) PSIS(L,K)=-15.
	   enddo
  8   CONTINUE
	  if(idanoflg.eq.1) CALL DANO
      RETURN
      END


      SUBROUTINE DANO

! ***   This is a routine to simulate dilution, dispersion, 
! *** or mixing of Nitrogen concentrations in the GOSSYM
! *** soil profile.  This routine will dilute or spread N
! *** in the cases of very sharp concentration gradients in
! *** wet soils.  This situation arises in drip fertigation
! *** and sidedress N applications.
! ***   This routine reads in a 2-D array of N concentrations and 
! *** a 2-D array of soil water contents, representing the GOSSYM 
! *** profile conditions.  This routine is a temporary way of 
! *** accounting for a N movement problem and is not a true  
! *** representation of dispersion-advection solute transport model.
! ***   This routine will be run after CAPFLO, when the updated
! *** N concentrations and soil water potentials are known.

      use common_block

      dimension conc1(40,20),vh2o(40,20),cmassn(40,20)

      rfact1 = .1
      do i = 1,nl   
        do j = 1,nk   
           conc1(i,j) = vno3c(i,j)
           vh2o(i,j) = vh2oc(i,j)*vcell
           cmassn(i,j) = conc1(i,j)*vh2o(i,j)
           totaln = totaln + cmassn(i,j)
        enddo
      enddo

      do i = nl-1,2,-1
        do j = nk-1,2,-1
          rfact = rfact1 * ((vh2oc(i,j)-thtr(i))/(thts(i)-thtr(i)))
          aven = (cmassn(i-1,j)+cmassn(i,j-1)+cmassn(i+1,j)+  &
                  cmassn(i,j+1))/4.0
          cmassn(i,j) = (1-rfact)*cmassn(i,j) + rfact*aven
          vno3c(i,j) = cmassn(i,j)/vh2o(i,j)
        enddo
      enddo

      return
      end
  

      SUBROUTINE NITRIF
! *************************************************************
! *                                                           *
! *  SIMPLIED VERSION BASED ON KAFKAFI,BAR-YOSEF AND HADAS    *
! *  MODEL. SOIL SCI 125:261-268. (1978).                     *
! *                                                           *
! *************************************************************
 
      use common_block
 
      NLH = NL/2
      DO 140 L=1,NLH
        M=LYRDPH(L)
        TT = TSOLAV(L)
        TFMIN = CONSK(L) + BETAK(L) * EXP(TT/10.)
        TFNIT = (10.0**12.02) * (10.0**(-3573.0/(TT+273)))
        DO 120 K=1,NK

! *** MINERALIZATION OF ORGANIC MATTER AND UREA. 

! *** TFMIN AND WFMIN ARE TEMPERATURE AND WATER DEPENDENT 
! *** FACTORS AFFECTING	THE NITROGEN MINERALIZATION RATE. 
! *** TFMIN WAS DERIVED BY L. BRIONES. WFMIN WAS DERIVED 
! *** FROM THE DATA OF MILLER AND JOHNSON, 1964. 
 
          IF(PSIS(L,K).GE.-0.07) THEN
            WFMIN = 0.58 - 0.30*LOG10(-100*PSIS(L,K))
          ELSE
            IF(PSIS(L,K).GT.-1.0) THEN
              WFMIN = -0.7 + 1.25*LOG10(-100*PSIS(L,K))
            ELSE
              WFMIN = 1.5 - 0.266*LOG10(-100*PSIS(L,K))
            ENDIF
          ENDIF
          IF(WFMIN.GT.1.0)WFMIN = 1.0
          IF(WFMIN.LT.0.0)WFMIN = 0.0
          DMINF = TFMIN * WFMIN 
          IF(DMINF.GT.1.0) DMINF = 1.0
          IF(DMINF.LT.0.0) DMINF = 0.0
          DNMIN = VNC(L,K) * DMINF
          VNC(L,K) = VNC(L,K) - DNMIN
          VNH4C(L,K) = VNH4C(L,K) + DNMIN
          ORGN = ORGN + (DNMIN*vcell/rowsp/0.011219)
 
! *** NITRIFICACION OF NH4. 

! *** TFNIT AND WFNIT ARE TEMPERATURE AND WATER DEPENDENT FACTORS 
! *** AFFECTING THE NITRIFICATION RATE. TFNIT COMES FROM KAFKAFI  
! *** ET.AL. WFNIT WAS DERIVED FROM THE DATA OF SABEY, 1969. 
 
          IF(PSIS(L,K).GT.-0.07) THEN
            WFNIT = 0.051 + 1.125*LOG10(-100*PSIS(L,K))
          ELSE
            WFNIT = 1.341 - 0.389*LOG10(-100*PSIS(L,K))
          ENDIF

          IF(WFNIT.GT.1.0)WFNIT = 1.0
          IF(WFNIT.LT.0.0)WFNIT = 0.0

          DNITF = TFNIT * WFNIT * 0.05 * 15.0
		   
! *** 15.0 IN THE ABOVE EQUATION COMES FROM J.VARCO & FDW FIELD OBS.

          IF(DNITF.GT.1.0) DNITF = 1.0
          IF(DNITF.LT.0.0) DNITF = 0.0
          DNIT = VNH4C(L,K) * DNITF
          VNH4C(L,K) = VNH4C(L,K) - DNIT
          IF(VNH4C(L,K).LE.0.000001)VNH4C(L,K)=0.0
          VNO3C(L,K) = VNO3C(L,K) + DNIT
  120   CONTINUE
  140 CONTINUE
      RETURN
      END


      SUBROUTINE DEFOLIAT
! *************************************************************************
! *                                                                        *
! *      DEFOLIAT SUBROUTINE.  SIMULATES THE EFFECTS OF COMMONLY USED      *
! *        DEFOLIATES AND PREP ON COTTON.                                  *
! **************************************************************************

      use common_block

! IF TODAY IS DATE OF PREP APPLICATION THEN WE MUST CALCULATE % INTERCEPTED

      DO I=1,5
        IF(DAYNUM.EQ.PRPDATE(I)) THEN
          IF(PRPMTH(I).EQ.0) THEN
            PRPKGH=PRPKGH+PRPPPA(I)*.95*1.12085*.75
          ELSE
            PRPKGH=PRPKGH+PRPPPA(I)*INT*1.12085*.75
          ENDIF
          TDFKGH=DEFKGH+PRPKGH
        ENDIF

        IF(DAYNUM.EQ.DEFDATE(I)) THEN
          IF(DEFMTH(I).EQ.0) THEN
            DEFKGH=DEFKGH+DEFPPA(I)*.95*1.12085*.75
          ELSE
            DEFKGH=DEFKGH+DEFPPA(I)*INT*1.12085*.75
          ENDIF
          TDFKGH=DEFKGH+PRPKGH
        ENDIF
      ENDDO
      IF((DEFBGN.GT.0).AND.(DAYNUM.EQ.DEFBGN)) THEN
        LFATDF = 0
        DO 70 I=1,LEFCNT
          IDUM = LEFSRT(I)
          IF(IDUM.GT.0) THEN
            K = IDUM/1000
            L = IDUM/10-K*100
            M = IDUM-IDUM/10*10
            IF(LEAFW(K,L,M).LT..0001) THEN
              LEFSRT(I) = 0
            ELSE
              LFATDF = LFATDF+1
            ENDIF
          ENDIF 
   70   CONTINUE
      ENDIF

      IF((PRPDAY.GT.0).AND.(DAYNUM.GT.PRPDAY))	 &
          AVGTSP=(AVGTSP*(DAYNUM-(PRPDAY+1))+TAVG)/(DAYNUM-PRPDAY)

! *** The constant 4. in the following logic is from DNB,VRR,FDW    
! *** field obs. at MITCHENERS 1989

      IF((DEFBGN.GT.0).AND.(DAYNUM.GT.DEFBGN)) THEN
        AVGTSD=(AVGTSD*(DAYNUM-(DEFBGN+1))+TAVG)/(DAYNUM-DEFBGN)
        DUM = PSILD*10.
        PERDEF = -35.2076+.5255*AVGTSD+7.0586*TDFKGH+				 &
                 1.7510*(DAYNUM-DEFBGN)-2.476*DUM-.03744*DUM**2+	 &
                 .0004198*AVGTSD*TDFKGH*(DAYNUM-DEFBGN)*DUM
        IF((DEFDAY.GT.0).AND.(DAYNUM.GT.DEFDAY)) THEN
          PERDEF = PERDEF*4.
        ELSE
          PERDEF = PERDEF*1.5
        ENDIF
        IF(PERDEF.LT.0.) PERDEF = 0.
        IF(PERDEF.GT.100.) PERDEF = 100.
      ENDIF

      RETURN
      END


      SUBROUTINE PIX
!  *************************************************************************
!  *                                                                       *
!  *              KRREDDY's PIX  SUBROUTINE                  KIT 5/12/99   *
!  *************************************************************************
!  *   THIS SUBROUTINE CALCULATES THE EFFECTS OF MEPIQUAT CHLORIDE         * 
!  * ON DIFFERENT PHYSIOLOGICAL PROCESSES OF COTTON INCLUDING GROWTH AND   *
!  * DEVELOPMENT.  THE DATABASE WAS DEVELOPED USING BASF CORPORATION'S     *
!  * GROWTH REGULATOR PIX AND USDA/ARS SOIL PLANT ATMOSPHERIC RESEARCH     *
!  * UNITS AT MISS.STATE.  THE DATABASE WAS DEVELOPED USING                *
!  * 20 G PIX AI/ACRE AND  AT A RANGE OF TEMPERATURE TREATMENTS.           *
!  *************************************************************************

      use common_block

! IF AN APPLICATION OF PIX HAS BEEN MADE, CALCULATE INTERCEPTED PORTION PIX 
! AND CONVERT INTO SI UNITS (G). TAKE THE AMOUNT PIX APPLIED AND ADD TO THE  
! EXISTING CONCENTRATION OF PIX IN THE PLANT.               

      PIXLOS=DEAD2DAY*PIXCON
      PIXPLT=PIXPLT-PIXLOS
      IF(DAYNUM.EQ.PIXDAY(IPIX)) THEN
        IF(PIXMTH(IPIX).EQ.0) THEN
          PIXPLT = PIXPLT + (19068.0 * PIXPPA(IPIX)*.90) / POPPLT
        ELSE
          PIXPLT = PIXPLT + (19068.0 * PIXPPA(IPIX)*INT) / POPPLT
        ENDIF
        IPIX = IPIX + 1
      ENDIF

! *** BASED ON FALLEN LEAF WEIGHT,SQUAR WEIGHT, BOLL WEIGHT AND DEAD 
! *** ROOT WEIGHT PIXLOS IS CALCULATED AND PIX CONCENTRATION IN THE 
! *** PLANT IS REDUCED EVERDAY IN OTHER SUBROUTINES OF THE PROGRAM.
! *** CALCULATE PIX CONCENTRATION IN COTTON PLANT ON A DRYWEIGHT BASIS 
! *** EVERDAY. ASSUME PIX CONCENTRATION IS SAME IN ALL PLANT PARTS. 
! *** PIX IS SYSTEMIC IN NATURE AND EASILY MOVABLE WITH IN THE PLANT.

	if(plantw.gt.0.0) then
	   pixcon = pixplt / plantw
	else
	   pixcon = 0.0
	endif

! *** CALCULATE STRESS PARAMETERS FOR HEIGHT,NODES AND LEAFAREA AS A 
! *** FUNCTION OF DAILY AVERAGE TEMPERATURE AND PIX CONC. IN PLANT TISSUE.

	conc = pixcon
	if(conc.gt.0.035 ) conc = 0.035
	pixdz = 1.0 - 18.619 * conc
	pixda = 1.0 - 7.774 * conc
	if(conc.gt.0.02 ) conc = 0.02
	pixdpn = 1.0 - 25.167*conc + 590.286*conc*conc

      IF(PIXDZ.GE.1.0)PIXDZ=1.0
      IF(PIXDZ.LE.0.45)PIXDZ=0.45
      IF(PIXDA.GE.1.0)PIXDA=1.0
      IF(PIXDA.LE.0.7)PIXDA=0.7
      IF(PIXDPN.GE.1.0)PIXDPN=1.0
      IF(PIXDPN.LE.0.7)PIXDPN=0.7

      RETURN
      END      


      SUBROUTINE PNET
! ************************************************************
! *                                                          *
! *                    PNET  SUBROUTINE                      *
! *                                     Reddy's new PNET     *
! ************************************************************
 
      use common_block

! PSILIN IS AN INDEX PSI(L) AT 0.3 BARS.
! PSIS AND 27 C. - FULLY TURGID
! PINDEX IS AN INDEX NET P RATE(SAME CONDITIONS).
 
      AVGPSI = PSIAVG

! PSILD IS MINIMUM LEAF WATER POTENTIAL FOR THE DAY. DATA LEADING TO
! PSILD ARE FROM KHARCHE'S DATA. PHD THESIS 1984 MSU.
 
      IF(AVGPSI.LT.-7.0) AVGPSI = -7.0
      DAYTMP = TDAY
      IF(DAYTMP.LT.25.) DAYTMP = 25.
 
! CALCULATE LEAF WATER POTENTIAL DURING THE DAY TIME
 
      IF((VARITY(IVARTY).EQ.'ST1').OR.(VARITY(IVARTY).EQ.'st1').OR.		&
         (VARITY(IVARTY).EQ.'ST2').OR.(VARITY(IVARTY).EQ.'st2').OR.		&
         (VARITY(IVARTY).EQ.'ST3').OR.(VARITY(IVARTY).EQ.'st3').OR.		&
         (VARITY(IVARTY).EQ.'HS26').OR.(VARITY(IVARTY).EQ.'hs26')) THEN
        CALL TEXASPSI
        IF(PSILD.LT.-3.5) PSILD=-3.5
      ELSE
        PSILD = (-1.2992-17.0305*AVGPSI-0.8148*DAYTMP					&
                +0.02173*RI+0.3371*AVGPSI*DAYTMP+0.02019*RI*AVGPSI		&
                -0.00078985*RI*DAYTMP+0.019766*DAYTMP**2)*0.101325		
 
! THE FOLLOWING IF STATEMENTS KEEPS THE PSILD EQUATION WITHIN THE 
! THE RANGE OF KHARCHE'S DATA
 
        PSILD = AMIN1(PSILD,(PSIAVG - 8.0)*0.101325)
      ENDIF
      IF(PSILD.LT.-3.5) PSILD=-3.5
      IF(PSILD.GE.-0.8) PSILD = -0.8
        
!  DATA LEADING TO THIS PTSRED ARE FROM CONTROLLED ENVIRONMENT UNITS
!  (SPAR) THE EFFECT OF WATER STRESS ON CANOPY SENESCENCE AND
!  APPARENT PHOTOSYNTHESIS IN COTTON (MARINI,BAKER,REDDY,AND McKINION
!  CROP SCI 25:798-802 (1985). PSTAND, RSUBL, RSUBO, GSUBR FROM
!  BAKER ET. AL. (1972). SIMULATION OF GROWTH AND YIELD IN COTTON:
!  I. GROSS PHOTOSYNTHESIS, RESPIRATION AND GROWTH. CROP SCI.12:431-435.

! *** Raja's new equation on the effect of water stress; r2 = 0.63   4/20/2001

      stdpsild = -0.8
      PnIndex = (4.7929 - 0.1318*stdpsild - 0.3319*stdpsild**2)
	  PnCurrent =  (4.7929 - 0.1318*psild - 0.3319*psild**2)
	  ptsred = PnCurrent/PnIndex 						
	  ptsred = PnCurrent/PnIndex * calbrt(23)							
      IF(PTSRED.GT.1.0) PTSRED = 1.0
      IF(PTSRED.LT.0.2) PTSRED = 0.2

! *** Nitrogen effect on PNET from unpublished data collected in 
! *** 1992 from pad experiment by K.R. Reddy (PIMA S-6)  3/10/2000

      LEAFCN = leafcn * 100.0
	  if(leafcn.ge.4.5)leafcn = 4.5
      XNIRED = -0.4029+0.5954*LEAFCN-0.0630*LEAFCN**2
      IF(XNIRED.GT.1.0) XNIRED = 1.0
 	  IF(XNIRED.LT.0.2) XNIRED = 0.2

      stress_index = amax1(0.0,ptsred*xnired)
 	  IF(stress_index.lt.0.15) stress_index = 0.15
      IF(stress_index.GT.1.) stress_index = 1.0

      PSTAND =(2.3908 + WATTSM*(1.37379 - WATTSM*0.00054136))*calbrt(24)
      PPLANT = PSTAND * INT *POPFAC *stress_index *0.001 *PNETCOR

!  VALUES BASED ON DATA OF HARPER ET. AL. (1973) CARBON DIOXIDE AND
!  THE PHOTOSYNTHESIS OF FIELD CROPS.  A METERED CARBON DIOXIDE
!  RELEASE IN COTTON UNDER FIELD CONDITIONS.  AGRON. JOUR. 65:7-11.
!  AND ON BAKER (1965)  EFFECTS OF CERTAIN ENVIRONMENTAL FACTORS
!  ON NET ASSIMILATION IN COTTON.  CROP SCI. 5:53-56. FIG 5.

      IF(CO2.EQ.1)PPLANT=PPLANT*1.405

! CO2 IS A FERTILIZATION TRIGGER.WHEN CO2 IS EQUAL TO 1,PPLANT IS
! INCREASED 40.5 % DUE TO  500 PPM CO2 CONCENTRATION.
! This part is not clear. CO2 is a input in ppm everywhere else, 6/1/06, SK

      RSUBL=0.0032125+0.0066875*TDAY
      LYTRES = RSUBL*PPLANT
      BMAIN=(PLANTW-COTXX)*RSUBO
      PTS=PPLANT-LYTRES-BMAIN
      IF(PTS.LE..01)PTS=.01
      PN= (PTS/(1.+GSUBR) * 0.68182) * PIXDPN
      SPN = SPN + PN

!  0.68182 CONVERTS CO2 TO CH2O

      RETURN
      END


      SUBROUTINE TEXASPSI

      use common_block

      RSMAX = RI*41860*4.363E-4/DAYLNG
      VPDMAX = (VP(TMAX)-VP(TMIN))/10.
      DOWN = 70.
      CALL AVGPSIS(DOWN,PSIS70)
      DOWN = 150.
      CALL AVGPSIS(DOWN,PSIS150)
      PSILD = -3.053+5.407*PSIS70-0.370*VPDMAX+0.004*RSMAX+5.878*PSIS150
      RETURN
      END


      SUBROUTINE AVGPSIS(DOWN,PSISDUM)

      use common_block

      XDUM = 0.
      KNT = 0
      MCELL = DOWN/DCELL
      IF(MCELL.GT.37) MCELL=37
      JDUM = MCELL+3
      DO L=1,JDUM
        DO K=1,4
          XDUM = XDUM+PSIS(L,K)
          KNT = KNT+1
        ENDDO
      ENDDO
      PSISDUM = XDUM/KNT/10
      RETURN
      END


      SUBROUTINE GROWTH
! ************************************************************
! *                                                          *
! *                  GROWTH SUBROUTINE                       *
! *                                            kit 5/11/99   *
! ************************************************************

      use common_block
 
      DAYTYM = DAYLNG / 24.
      NYTTYM = 1. - DAYTYM
      IF(TDAY.LT.13.5)TDAY=13.5
      IF(TNYT.LT.13.5)TNYT=13.5
      SPDWSQ = 0.0
      SPDWBO = 0.0
      SPDWLD = 0.
      SPDWLN = 0.
      PFAREA = 0.
      AREA = 0.

! *** Water stress factor for stem growth
 
      call stem_growth_water_stress

! *** Temperature and water stress factors affecting leaf and fruit growth
 
      call leaf_growth_water_stress
      call fruit_growth_water_stress

! *** Calculate potential growth of prefruiting leaves  

	  call PfLeafArea
  
! *** Calculate potential square (PDWSQ) and boll (PDWBOL) growth 

      IF(FCODE(1,1,1).gt.0) then 
         FRATIO = GBOLWT/PLANTW
         DO K=1,NVBRCH
            NBRCH = NFBR(K)
            DO L=1,NBRCH
               NNID = NNOD(K,L)
               DO M=1,NNID
		        call pot_sqr_growth(k,l,m) 
		        call pot_bol_growth(k,l,m) 
	         enddo
	      enddo
	   enddo
 
! *** Calculate potential main stem leaf growth. Assume that main  
! *** stem leaf is initiated at the same time as fruiting branch leaf. 
 
         call MSLeafArea
         call FBLeafArea
      endif

! *** CALCULATE POTENTIAL STEM GROWTH

	  call pot_stm_growth 

! *** Call RUTGRO to calculate potential root growth (SPDWRT).

      CALL RUTGRO(0)

!      SPDWRT = SPDWRT * POPFAC

      SPDWRT = SPDWRT * POPFAC
      SPDWRD = SPDWRT*DAYTYM
      SPDWRN = SPDWRT*NYTTYM

! *** Calculate root/shoot correction factors for potential root, 
! *** stem amd leaf growth . Source: AVI BEN-PORATH PHD THESIS.

 	  call rootshootratio

! *** Calculate carbohydrate demand of each organ and the whole plant

      CDSTEM = PDSTMD + PDSTMN
      CDROOT = SPDWRD + SPDWRN
      CDLEAF = SPDWLD + SPDWLN
      CDSQAR = SPDWSQ
      CDBOLL = SPDWBO
      CD = CDSTEM + CDLEAF + CDROOT + CDSQAR + CDBOLL
      CPOOL = PN + RESC
 
! *** Calculate carbohydrate stress (CSTRES)
 
      call Carbon_stress
 
! *** Reduce the potential weight of each organ 
 
      PDSTEM = CDSTEM * CSTRES
      PDLEAF = CDLEAF * CSTRES
      PDROOT = CDROOT * CSTRES
      PDBOLL = CDBOLL * CSTRES
      PDSQ   = CDSQAR * CSTRES

      IF(PDSTEM.LT.0.) PDSTEM=0.
      IF(PDLEAF.LT.0.) PDLEAF=0.
      IF(PDROOT.LT.0.) PDROOT=0.
      IF(PDBOLL.LT.0.) PDBOLL=0.
      IF(PDSQ.LT.0.) PDSQ=0.
 
! *** Call NITRO to get today's N stress

      CALL NITRO
      IF(LEAFWT.gt.0.0) then
 
! *** Calculate excess carbohydrate created by N stress
 
         CH2OX = (CDSTEM+CDLEAF)*CSTRES*(1.-NV) +			&
                 (CDSQAR+CDBOLL) *CSTRES*(1.-NF)  +			&
                 (CDROOT*CSTRES *(1.-NR))
         XTRAC = XTRAC + CH2OX
 
! *** Initialize LEAFWT to 0.2 (Initial leaf weight). 
! *** Assume cotyledons fall off at first square.

         IF(ISQ.GT.0) THEN
            LEAFWT = 0.
         ELSE
            LEAFWT = 0.2
!            LEAFWT = 0.02
         ENDIF
  
! *** Calculate nutrition stress

         VSTRES = NV * CSTRES
         FSTRES = NF * CSTRES
         RSTRES = NR * CSTRES

         SDWSQR = 0.0
         SDWBOL = 0.0
         SDWLEF = 0.0
  
! *** Calculate actual growth of organs
 
         IF(FCODE(1,1,1).gt.0) then
            SQWT = 0.0
            GBOLWT = 0.0
            DO K = 1,NVBRCH
               NBRCH = NFBR(K)
               DO L = 1,NBRCH
                  NNID = NNOD(K,L)
                  DO M = 1,NNID
 
! *** Calculate actual square growth
 
                     call Actual_square_growth(k,l,m)
 
! *** Calculate actual boll growth
 
				   call Actual_boll_growth(k,l,m)
 
! *** Calculate actual main stem and fruiting branches leaf growth

				   call Actual_leaf_growth(k,l,m)

! *** Calculate actual main stem and fruiting branches leaf areas

				   call Actual_leaf_area(k,l,m)
 100              enddo
 101           enddo
 102        enddo
	   endif
 
! *** CALCULATE prefruiting LEAF GROWTH AND NEW LEAF AREA
 
 26      DO J=1,NUMPFN
            call Actual_prefruitingleaf_growth(j)
            call Actual_prefruitingleaf_area(j)
 28      enddo

         AREA = AREA + PFAREA
         LAI = AREA/POPFAC
 
! *** CALCULATE stem and root GROWTH 
 
         SDWSTM = CDSTEM * VSTRES
         RCH2O  = CDROOT * RSTRES 

! *** SLEAF, SBOLL, SQUAR, SSTEM ARE THE CUMMULATIVE AMOUNT 
! *** OF CH0 ALLOCATED TO THE DIFFERENT SINKS.
 
         SLEAF = SLEAF + SDWLEF
         SBOLL = SBOLL + SDWBOL
         SQUAR = SQUAR + SDWSQR
         SSTEM = SSTEM + SDWSTM
         SROOT = SROOT + RCH2O

         STEMWT = STEMWT + SDWSTM
         STMWT(KDAY) = STEMWT

         CALL RUTGRO(1)

! *** NOTE THAT STORED CH2O IS USED IN CALCULATION OF TOTAL 
! *** PLANT WEIGHT.  CALCULATE VERTICAL GROWTH

 40      IF (PLANTW.LE.1.0) VSTRES = 1.

         CALL PLANTHEIGHT

 510     FORMAT(1X,I5,10F7.2)
	  endif

      RETURN
      END


      subroutine stem_growth_water_stress
! ************************************************************
! *   Based on Marani et.al., 1985 Effects of water stress   *
! * on canopy senescence and carbon exchange rates in cotton *
! * Crop Sci. 25:798-802.                                    *
! *                                            kit 12/14/99  *
! ************************************************************

      use common_block
 
      PSILIN = CALBRT(16)
      WSTST1 = 13.16 + 9.007*PSILIN + 1.452*PSILIN*PSILIN
      WSTST2 = 13.16 + 9.007*PSILD + 1.452*PSILD*PSILD
      IF(WSTST2.LE.0.1) WSTST2=0.1
      WSTRST = WSTST2/WSTST1
      IF(WSTRST.GT.1.05) WSTRST = 1.05
      IF(WSTRST.LT.0.00) WSTRST = 0.00
      RETURN
      END


	  subroutine leaf_growth_water_stress
! ************************************************************
! *   Based on Marani et.al., 1985 Effects of water stress   *
! * on canopy senescence and carbon exchange rates in cotton *
! * Crop Sci. 25:798-802. Leaf growth is affected by water   *                                     *
! * stress at <= -1.0 bar. 	rsquare = 0.73						       *
! *                                            kit 12/14/99  *
! ************************************************************

      use common_block
 
      PSILIN = CALBRT(33)
      WSTRL1 = 44.89+33.98*PSILIN+6.38*PSILIN*PSILIN
      WSTRL2 = 44.89+33.98*PSILD+6.38*PSILD*PSILD
      IF(WSTRL2.LE.0.1) WSTRL2=0.1
      WSTRLF = WSTRL2/WSTRL1
      IF(WSTRLF.GT.1.05) WSTRLF = 1.05
      IF(WSTRLF.LT.0.00) WSTRLF = 0.00
      day_lfstress = DAYTYM * WSTRLF * PIXDA * calbrt(32)
      eve_lfstress = NYTTYM * WSTRLF * PIXDA * calbrt(32)
      RETURN
      END

 
      subroutine fruit_growth_water_stress
! *************************************************************
! *   Based on Marani and Baker, 1981 Science Report to US    *
! * Israel Binational Foundation. Development of a Predictive *
! * Dynamic Simulation Model of growth and yield in Acala.    *
! *                                            kit 12/14/99   *
! *************************************************************

      use common_block

      WSTRSD = -2.5/(PSIAVG - 1.6)+(0.0005*PSIAVG*TDAY) - 0.001*RN
      WSTRSN = -2.5/(PSIAVG - 1.6)+(0.0005*PSIAVG*TNYT)
      IF(WSTRSD.LT.0.0) WSTRSD = 0.0
      IF(WSTRSD.GT.1.0) WSTRSD = 1.0
      IF(WSTRSN.LT.0.0) WSTRSN = 0.0
      IF(WSTRSN.GT.1.0) WSTRSN = 1.0
      WSTRS = (WSTRSD*DAYTYM) + (WSTRSN*NYTTYM)
      WFDAY1 = DAYTYM * CALBRT(52) * AMIN1(1.0,0.5+WSTRSD)
      WFNYT1 = NYTTYM * CALBRT(52) * AMIN1(1.0,0.5+WSTRSN)
      WFDAY2 = DAYTYM * CALBRT(57) * AMIN1(1.0,0.3+WSTRSD)
      WFNYT2 = NYTTYM * CALBRT(57) * AMIN1(1.0,0.3+WSTRSN)
	  eve_water_index = wfnyt1
	  day_water_index = wfday1

! *** Calculate water and temperature stress factors for fruit growth

      day_expfac = AMAX1(0.,(0.0160791*TDAY-0.2120865) * WFDAY2)
      eve_expfac = AMAX1(0.,(0.0160791*TNYT-0.2120865) * WFNYT2)
      day_hitemp = AMAX1(0.,(2.7328500-0.082857*TDAY)  * WFDAY1)
      eve_hitemp = AMAX1(0.,(2.7328500-0.082857*TNYT)  * WFNYT1)
      day_lotemp = AMAX1(0.,(0.0312500*TDAY-0.508125) * WFDAY1)
      eve_lotemp = AMAX1(0.,(0.0312500*TNYT-0.508125) * WFNYT1)
      RETURN
      END


	subroutine pot_sqr_growth(k,l,m)
! ************************************************************
! * CALCULATE POTENTIAL SQUARE GROWTH(PDWSQ); rsquare = 0.98 *
! ************************************************************

      use common_block
 
      IF(FCODE(K,L,M).EQ.1) THEN

! *** Raja Reddy's new square growth equation      5/17/2000 kit

	   if(avgt(k,l,m).lt.17.0) avgt(k,l,m) = 17.0 
	   if(avgt(k,l,m).gt.35.0) avgt(k,l,m) = 35.0 
	   PDWSQ(k,l,m) = -0.01017 + avgt(k,l,m) 				   &
	                   *(0.0012532-0.000019484*avgt(k,l,m)) 
         SPDWSQ = SPDWSQ + PDWSQ(K,L,M) * calbrt(17)
      ENDIF
      RETURN
      END


	subroutine pot_bol_growth(k,l,m)
! ************************************************************
! *     CALCULATE POTENTIAL BOLL GROWTH(PDWBOL) Original     *
! ************************************************************

      use common_block
 
! *** CALCULATE EXPONENTIAL BOLL GROWTH(PDWBOL)
 
      IF((FCODE(K,L,M).EQ.2).OR.(FCODE(K,L,M).EQ.7)) THEN
         PHASE1=DEHISS(K,L,M)*0.15
         IF(AGEBOL(K,L,M).LE.PHASE1) THEN
            PDWBOD(K,L,M) = BOLWGT(K,L,M)*day_expfac*FFRUT(K,L,M)
            PDWBON(K,L,M) = BOLWGT(K,L,M)*eve_expfac*FFRUT(K,L,M)
         ELSE
 
! *** D(LN(W))/DT = B = (1/W)(DW/DT) = F(TAVG) LINEAR BOLL GROWTH PHASE.
! *** ACTIVE BOLL GROWTH PHASE.
 
            C00=1.0
            IF((VARITY(IVARTY).EQ.'PIMA').OR.(VARITY(IVARTY).EQ.'pima')) &
	      	    C00=.75
            IF(TDAY.GT.28.5) THEN
               PDWBOD(K,L,M) = day_hitemp*FFRUT(K,L,M)*C00
               PDWBON(K,L,M) = eve_hitemp*FFRUT(K,L,M)*C00
            ELSE
               PDWBOD(K,L,M) = day_lotemp*FFRUT(K,L,M)*C00
               PDWBON(K,L,M) = eve_lotemp*FFRUT(K,L,M)*C00
            ENDIF
         ENDIF
         IF(FFRUT(K,L,M).GT.0.001)THEN
            BSIZE(K,L,M) = BOLWGT(K,L,M)/FFRUT(K,L,M)
         ELSE
            BSIZE(K,L,M) = 0.0
         ENDIF
         IF(BSIZE(K,L,M).GT.CALBRT(12)) THEN
            PDWBOD(K,L,M) = PDWBOD(K,L,M) * 0.25
            PDWBON(K,L,M) = PDWBON(K,L,M) * 0.25
         ENDIF
         SPDWBO = SPDWBO + PDWBOD(K,L,M) + PDWBON(K,L,M)
      ENDIF
      RETURN
      END


	subroutine pot_bol_growthN(k,l,m)
! ************************************************************
! *          CALCULATE POTENTIAL BOLL GROWTH(PDWBOL)         *
! ************************************************************

      use common_block
 
      IF((FCODE(K,L,M).EQ.2).OR.(FCODE(K,L,M).EQ.7)) THEN

 ! *** Raja Reddy's new boll growth equation      5/17/2000 kit

 	   if(tday.lt.17.0) tday = 17.0 
	   if(tday.gt.35.0) tday = 35.0 
 	   if(tnyt.lt.17.0) tnyt = 17.0 
	   if(tnyt.gt.35.0) tnyt = 35.0 
	   pdwbod(k,l,m) = (-0.23303 + tday*(0.02797-0.0005852*tday))	 &
	                    * FFRUT(K,L,M) * daytym * CALBRT(18)
!	.                     * FFRUT(K,L,M) * day_water_index * CALBRT(18)

	   pdwbon(k,l,m) = (-0.23303 + tnyt*(0.02797-0.0005852*tnyt)) 	 &
                        * FFRUT(K,L,M) * nyttym * CALBRT(18)
!	.                     * FFRUT(K,L,M) * eve_water_index * CALBRT(18)
         IF(FFRUT(K,L,M).GT.0.001)THEN
            BSIZE(K,L,M) = BOLWGT(K,L,M)/FFRUT(K,L,M)
         ELSE
            BSIZE(K,L,M) = 0.0
         ENDIF
         IF(BSIZE(K,L,M).GT.CALBRT(12)) THEN
            PDWBOD(K,L,M) = PDWBOD(K,L,M) * 0.25
            PDWBON(K,L,M) = PDWBON(K,L,M) * 0.25
         ENDIF
         SPDWBO = SPDWBO + PDWBOD(K,L,M) + PDWBON(K,L,M)
      ENDIF
      RETURN
      END


	subroutine pot_stm_growth
! ************************************************************
! *         CALCULATE POTENTIAL STEM GROWTH(PDSTMD/N)        *
! ************************************************************

      use common_block
 
      IF(KDAY.LE.CALBRT(13)) THEN
!         PDSTMD = (0.1 + 0.02 * KDAY) * DAYTYM * CALBRT(14)
!         PDSTMN = (0.1 + 0.02 * KDAY) * NYTTYM * CALBRT(14)
         FACDAY = DAYTYM *AMIN1(1.0,0.2+WSTRST) * CALBRT(14)
         FACNYT = NYTTYM *AMIN1(1.0,0.2+WSTRST) * CALBRT(14)
         PDSTMD = (0.1 + 0.02 * KDAY) * FACDAY
         PDSTMN = (0.1 + 0.02 * KDAY) * FACNYT 
      ELSE
!         PDSTMD=(0.2 + 0.06 * (STEMWT-STMWT(KDAY-32)))*DAYTYM*CALBRT(15)
!         PDSTMN=(0.2 + 0.06 * (STEMWT-STMWT(KDAY-32)))*NYTTYM*CALBRT(15)
         FACDAY = DAYTYM *AMIN1(1.0,0.2+WSTRST) * CALBRT(15)
         FACNYT = NYTTYM *AMIN1(1.0,0.2+WSTRST) * CALBRT(15)
         PDSTMD=(0.2 + 0.06 * (STEMWT-STMWT(KDAY-32))) * FACDAY
         PDSTMN=(0.2 + 0.06 * (STEMWT-STMWT(KDAY-32))) * FACNYT
      ENDIF
      RETURN
      END


	subroutine rootshootratio
! ************************************************************
! * Calculate correction factors to potential stem, leaf and *
! * root growth based on the root/shoot ratio                *
! ************************************************************

      use common_block

      REAL NRUT
 
! *** RATIO IS ROOTWT/(STEMWT+LEAFWT) RATIO UNDER NO WATER STRESS
! *** CONDITIONS. DATA FOR THE EQUATION BELOW COMES FROM AVI 
! *** BEN-PORATH PHD THESIS.

      RATIO  = 0.2072+0.606506*EXP(-0.13*(STEMWT+LEAFWT))

! *** RATIO (DURING THE DAY AND NIGHT) IS INCREASED AS A FUNCTION OF
! *** WATER STRESS. IT COULD BE INCREASED AS MUCH AS 25% IF WSTRS IS 0.

      RATIOD = RATIO * (1.00 + 0.25*(1.0-WSTRSD)) * CALBRT(22)
      RATION = RATIO * (1.00 + 0.25*(1.0-WSTRSN)) * CALBRT(22)

! *** RUTSUT IS THE SIMULATED ROOT/SHOOT RATIO.

      RUTSUT = ROOTWT/(STEMWT+LEAFWT)
      DTOP = 1.
      NTOP = 1.
      DRUT = 1.
      NRUT = 1.

! *** IF ROOT/SHOOT (RUTSUT) RATIO IS BELOW THE TARGET (RATIO) 
! *** INCREASED ROOT GROWTH POTENTIAL AND DECREASE STEM AND 
! *** LEAF GROWTH POTENTIAL.

      IF(RUTSUT.LT.RATIOD) THEN
         DRUT = RATIOD/RUTSUT
         DTOP = RUTSUT/RATIOD
         SPDWRD = SPDWRD * DRUT
         PDSTMD = PDSTMD * DTOP
         SPDWLD = SPDWLD * DTOP
      ENDIF
      IF(RUTSUT.LT.RATION) THEN
         NRUT = RATION/RUTSUT
         NTOP = RUTSUT/RATION
         SPDWRN = SPDWRN * NRUT
         PDSTMN = PDSTMN * NTOP
         SPDWLN = SPDWLN * NTOP
      ENDIF

! *** IF ROOT/SHOOT RATIO IS ABOVE THE TARGET RATIO DECREASE ROOT 
! *** GROWTH POTENTIAL AND INCREASE STEM AND LEAVES GROWTH POTENTIAL.
 
      IF(RUTSUT.GT.RATIOD) THEN
         DRUT = RATIOD/RUTSUT
         DTOP = RUTSUT/RATIOD
         SPDWRD = SPDWRD * DRUT
         PDSTMD = PDSTMD * DTOP
         SPDWLD = SPDWLD * DTOP
      ENDIF
      IF(RUTSUT.GT.RATION) THEN
         NRUT = RATION/RUTSUT
         NTOP = RUTSUT/RATION
         SPDWRN = SPDWRN * NRUT
         PDSTMN = PDSTMN * NTOP
         SPDWLN = SPDWLN * NTOP
      ENDIF
      RETURN
      END


      SUBROUTINE Carbon_stress
! ************************************************************
! ***         CALCULATE CARBOHYDRATE stress (CSTRES)		   *
! ************************************************************

      use common_block
 
      CSTRES = CPOOL / CD
      IF(CSTRES.LE.1.) THEN
         RESC = 0.
      ELSE
         CSTRES = 1.
         RESC = CPOOL - CD
         DUMY11 = 0.3*LEAFWT
         IF(RESC.GT.DUMY11) THEN
            XTRAC  =  XTRAC + RESC - DUMY11
            RESC   = DUMY11
            CSTORE = 1.0
         ELSE
            CSTORE=RESC/DUMY11
         ENDIF
      ENDIF
      RETURN
      END
 

      SUBROUTINE Actual_square_growth(k,l,m)
! ************************************************************
! ***             CALCULATE ACTUAL SQUARE GROWTH	    	   *
! ************************************************************

      use common_block
 
      IF(FCODE(K,L,M).EQ.1) THEN
         DWSQ = PDWSQ(K,L,M)*FSTRES
         SQRWT(K,L,M) = SQRWT(K,L,M) + DWSQ
         SDWSQR = SDWSQR + DWSQ
         SQWT = SQWT + SQRWT(K,L,M)
      ENDIF
      RETURN
      END
 

      SUBROUTINE Actual_boll_growth(k,l,m)
! ************************************************************
! ***             CALCULATE ACTUAL BOLL GROWTH	    	   *
! ************************************************************

      use common_block
 
      IF(FCODE(K,L,M).EQ.2.OR.FCODE(K,L,M).EQ.7) THEN
         DWBOLL = (PDWBOD(K,L,M)+PDWBON(K,L,M))*FSTRES
         BOLWGT(K,L,M) = BOLWGT(K,L,M) + DWBOLL
         SDWBOL = SDWBOL + DWBOLL
         GBOLWT = GBOLWT + BOLWGT(K,L,M)
      ENDIF
      RETURN
      END


      SUBROUTINE Actual_leaf_growth(k,l,m)
! ************************************************************
! ***             CALCULATE ACTUAL LEAF GROWTH	    	   *
! ************************************************************

      use common_block
 
! *** Calculate actual main stem leaf growth

      IF(M.EQ.1) THEN
         DWL = (PDWMLD(K,L)*DTOP+PDWMLN(K,L)*NTOP)*VSTRES
         SDWLEF = SDWLEF + DWL
         MLEAFW(K,L) = MLEAFW(K,L) + DWL
         LEAFWT = LEAFWT + MLEAFW(K,L)
      ENDIF

! *** Calculate actual fruiting branches leaf growth
 
      DWL  = (PDWFLD(K,L,M)*DTOP+PDWFLN(K,L,M)*NTOP)*VSTRES
      SDWLEF = SDWLEF + DWL
      LEAFW(K,L,M) = LEAFW(K,L,M) + DWL
      LEAFWT = LEAFWT + LEAFW(K,L,M)
      RETURN
      END


      SUBROUTINE Actual_leaf_area(k,l,m)
! ************************************************************
! ***             CALCULATE ACTUAL LEAF AREA  	    	   *
! ************************************************************

      use common_block

! *** Calculate actual main stem leaf area

      IF(M.EQ.1) THEN
         DAL = (PDAMLD(K,L)*DTOP+PDAMLN(K,L)*NTOP)*VSTRES
         MLAREA(K,L) = MLAREA(K,L) + DAL
         AREA = AREA + MLAREA(K,L)
      ENDIF

! *** Calculate actual fruiting branches leaf area

      DAL = (PDADAY(K,L,M)*DTOP+PDANYT(K,L,M)*NTOP)*VSTRES
 	  LAREA(K,L,M) = LAREA(K,L,M) + DAL
      AREA = AREA + LAREA(K,L,M)
      RETURN
      END


      SUBROUTINE Actual_prefruitingleaf_growth(j)
! ************************************************************
! ***          CALCULATE PREFRUITING LEAF GROWTH	    	   *
! ************************************************************

      use common_block
 
      PFDWL = (PFDWLD(J)*DTOP + PFDWLN(J)*NTOP) * VSTRES
      PFWL(J) = PFWL(J) + PFDWL
      SDWLEF = SDWLEF + PFDWL
      LEAFWT = LEAFWT + PFWL(J)
      RETURN
      END


      SUBROUTINE Actual_prefruitingleaf_area(j)
! ************************************************************
! ***          CALCULATE PREFRUITING LEAF AREA	    	   *
! ************************************************************

      use common_block
 
      PFDAL = (PFDALD(J)*DTOP + PFDALN(J)*NTOP) * VSTRES
	  PFAL(J) = PFAL(J) + PFDAL
      PFAREA = PFAREA + PFAL(J)
      RETURN
      END
 

      SUBROUTINE PLANTHEIGHT
! ************************************************************
! *   Based on Reddy, et. al., 1997.  Modeling temperature   *
! * effects on cotton internode and leaf growth.  Crop Sci.  *
! *                                            kit 5/11/99   *
! ************************************************************

      use common_block

      DZ = 0.
      if(tavg.lt.17.0) tavg = 17.0
      if(tavg.gt.36.0) tavg = 36.0

!*** rsquare = 0.96 for internode elongation duration equation

      DURATION = -0.04312+0.007383*TAVG-0.0001046*TAVG**2
      IF(DURATION.LT.0.) DURATION = 0.

! *** Calculate reduction to internode length due to defoliants

      IF((PRPDAY.GT.0).AND.(DAYNUM.GT.PRPDAY)) THEN
        DUM = DAYNUM-PRPDAY
        DZCON = .07151+.008617*TAVG-.01844*DUM+.001176*TAVG*DUM
        DZPREP = -.4685+.04319*TAVG-.005356*DUM+.3196*			&
                 PRPKGH-.1063*PRPKGH**2-.0006625*TAVG**2
        DZPREP = DZPREP/DZCON
        IF(DZPREP.GT.1.5) DZPREP=1.5
        IF(DZPREP.LT..001) DZPREP=.001
      ELSE
        DZPREP = 1.
      ENDIF

! *** Reduction in plant height due to water stress from SPAR data 
! *** collected by K. R. Reddy in 1991; 4 cm is the unstressed height  
! *** rsquare = 0.77 for the equation


	if(psild.le.-2.5) psild = -2.5
      H2OIndex = (4.2904 - PSILD*(0.6491+0.9737*PSILD))/4.
      IF(H2OIndex.LT.0.15) H2OIndex = 0.15
      IF(H2OIndex.GT.1.0) H2OIndex = 1.0

      N00 = NUMPFN + NFBR(1)
      DO I=1,N00
         XMNODAGE(I) = XMNODAGE(I)+DURATION * CALBRT(44) 
         IF(XMNODAGE(I).LT.1.0) THEN
	      xnodage(i) = xnodage(i) + 1.0 
	      cntlfcn(i) = cntlfcn(i) + 1.0
	      sumlfcn(i) = sumlfcn(i) + leafcn

! *** Calculate initial internode length

            IF(XMNODLTH(I).LT.0.0001) THEN

! *** Calculate reduction to initial internode length due to temperature
! *** Assume that maximum initial internode length at 27oC = 0.75 cm.
         
	    	rednodlth = (-0.6853 + tavg*(0.1077-0.002031*tavg))/0.75

! *** rsquare = 0.11 for above equation

	        if (rednodlth.lt.0.001) rednodlth = 0.001
	        if (rednodlth.gt.1.0) rednodlth = 1.0
              IF(i.lt.14) THEN
                 XMNODLTH(I)=(0.05738+0.05605*I)*rednodlth*CALBRT(47)

! *** rsquare = 0.93 for above equation

                 IF(XMNODLTH(I).LT.0.1) XMNODLTH(I) = 0.1
              ELSE
                 XMNODLTH(I)=(1.3589-0.0407*I)*rednodlth * CALBRT(48)

! *** rsquare = 0.91 for above equation

                 IF(XMNODLTH(I).LT.0.3) XMNODLTH(I) = 0.3
              ENDIF

! *** Calculate reduction to initial internode length due to nitrogen 

			if(i.gt.1) then
	           avelfcn = sumlfcn(i-1) / cntlfcn(i-1)
	           faclfn = avelfcn / 0.035
			   if(faclfn.gt.1.0) faclfn = 1.0
		   	   redlfcn = -2.9953 + faclfn*(8.0622-4.0669*faclfn)
 		       if(redlfcn.lt.0.2) redlfcn = 0.2
		       if(redlfcn.gt.1.0) redlfcn = 1.0
	        else 
		       redlfcn = 1.0
	        endif
	        stress_index = H2OIndex * redlfcn * pixdz * dzprep
              if(stress_index.lt.0.15) stress_index = 0.15
              if(stress_index.GT.1.) stress_index = 1.0
              DZ = DZ+XMNODLTH(I) * stress_index  * CALBRT(45) 
            ELSE
	        faclfn = leafcn/0.035
			if(faclfn.gt.1.0) faclfn = 1.0
			redlfcn = -2.9953 + faclfn*(8.0622-4.0669*faclfn)
 		    if(redlfcn.lt.0.2) redlfcn = 0.2
		    if(redlfcn.gt.1.0) redlfcn = 1.0
	        stress_index = H2OIndex * redlfcn * pixdz * dzprep
              if(stress_index.lt.0.15) stress_index = 0.15
              if(stress_index.GT.1.) stress_index = 1.0

! *** Calculate the daily increment of internode length elongation
! *** Calculate the relative rate of elongation (yintrcpt)
! *** and rate of reduction (slope).      rsquare = 0.97

              YINTRCPT = -0.001427 + 0.0166*TAVG
              SLOPE = 0.02479 - 0.001994*TAVG

! *** Calculate the current rate of internode length elongation 

              Stem_Elongation_Rate = YINTRCPT + SLOPE*XNODAGE(I)
			if(Stem_Elongation_Rate.le.0.0)Stem_Elongation_Rate = 0.0

! *** Calculate the current internode length 

              Potential_Height = Stem_Elongation_Rate * XMNODLTH(I)

              GROWFAC = Potential_Height * stress_index * CALBRT(49)
              IF(GROWFAC.gt.0.6) GROWFAC = 0.6
              DZ = DZ + GROWFAC	
              XMNODLTH(I) = XMNODLTH(I) + GROWFAC 
            ENDIF
         ENDIF
      ENDDO

      Z = Z + DZ

      RETURN
      END


      SUBROUTINE PfLeafArea

      use common_block

      do i=1,numpfn

! ***	Calculate the leaf age

         PfLfAge(i)=PfLfAge(i)+Duration_leaf_expansion(tavg)*calbrt(40) 
         if(PfLfAge(i).lt.1.0) then
	      agepflf(i) = agepflf(i) + 1.0
	      cntpflfcn(i) = cntpflfcn(i) + 1.0
	      sumpflfcn(i) = sumpflfcn(i) + leafcn

! *** Initiated once the node is produced (ADDPFNOD routine) except 
! *** for the first node.  Calculate the initial leaf area at unfolding

            if(PfLfArea(i).le.0.001) then
	
! *** Calculate running average of leaf nitrogen (max leafcn = 0.035)

	         if(i.eq.1) then
			    avelfcn = 0.035
	         else
			    avelfcn = sumpflfcn(i-1) / cntpflfcn(i-1)
	         endif

! *** Convert leaf area, cm2 to dm2 and accumulate the daily increment
! *** in leaf area to determine the total current leaf area 
! *** PfAL is the initial leaf area at unfolding = PfLfArea

               call Potential_PreFruiting_Leaf_Area(i)
               PfAL(i) = PfLfArea(i)*CNReduction_to_Leaf_Area(i,avelfcn)  &
	                     * calbrt(41)
	         raday = 1.0
	         ranyt = 1.0
            else

! *** Calculate the current rate of leaf expansion during the day and night

               raday = (Rate_Leaf_expansion(tday) + Rate_of_reduction(tday)*agepflf(i))
               ranyt = (Rate_Leaf_expansion(tnyt) + Rate_of_reduction(tnyt)*agepflf(i))
               if(raday.LT.0.003) raday = 0.003
               if(ranyt.LT.0.003) ranyt = 0.003
            endif

 ! *** Calculate the current leaf area 

            raday = raday * PfAL(i)
            ranyt = ranyt * PfAL(i)

! *** Calculate the corresponding daily potential change in leaf weight

            PFDALD(i) = raday * day_lfstress * calbrt(34)	 &
	            	  * CNReduction_to_Leaf_Area(i,leafcn)
            PFDALN(i) = ranyt * eve_lfstress * calbrt(34)	 &
		              * CNReduction_to_Leaf_Area(i,leafcn)
            PFDWLD(i) = PFDALD(i) * Sp_leaf_weight(tday) * CALBRT(36) 
            PFDWLN(i) = PFDALN(i) * Sp_leaf_weight(tnyt) * CALBRT(36) 
            SPDWLD    = SPDWLD + PFDWLD(i)
            SPDWLN    = SPDWLN + PFDWLN(i)
	   else
            PFDALD(i) = 0.0
            PFDALN(i) = 0.0
            PFDWLD(i) = 0.0
            PFDWLN(i) = 0.0
         endif
      enddo
      return
      end


      SUBROUTINE MSLeafArea

      use common_block

      do k=1,nvbrch
         nbrch = nfbr(k)
         do l=1,nbrch

! ***	Calculate the leaf age

            XMSLfAge(k,l) = XMSLfAge(k,l) + Duration_leaf_expansion(tavg)*calbrt(40) 
            if(XMSLfAge(k,l).lt.1.0) then
	         agemslf(k,l) = agemslf(k,l) + 1.0
		     nodecount = numpfn + l
	         cntmslfcn(nodecount) = cntmslfcn(nodecount) + 1.0
	         summslfcn(nodecount) = summslfcn(nodecount)+leafcn

! *** Initiated once new node is produced (ADDMSNOD routine)

               if(XMSLfArea(k,l).le.0.001) then
                  call Potential_MainStem_Leaf_Area(k,l)
	            if(cntmslfcn(nodecount-1).le.0.0) cntmslfcn(nodecount-1) =1.0
	            avelfcn = summslfcn(nodecount-1)/cntmslfcn(nodecount-1)
                  MLArea(k,l) = XMSLfArea(k,l) * calbrt(42)				  &
				         * CNReduction_to_Leaf_Area(nodecount,avelfcn)
	            raday = 1.0
	            ranyt = 1.0
               else

! *** Calculate the current rate of leaf expansion during the day and night

                  raday = (Rate_Leaf_expansion(tday) 			  &
			               + Rate_of_reduction(tday)*agemslf(k,l))
                  ranyt = (Rate_Leaf_expansion(tnyt) 			  &
			               + Rate_of_reduction(tnyt)*agemslf(k,l))
                  if(raday.LT.0.003) raday = 0.003
                  if(ranyt.LT.0.003) ranyt = 0.003
               endif

! *** Calculate the current leaf area 

               raday = raday * MLArea(k,l)
               ranyt = ranyt * MLArea(k,l)

! *** Calculate the corresponding daily potential change in leaf weight

               PDAMLD(K,L) = raday * day_lfstress * calbrt(35)		  &
	            	     * CNReduction_to_Leaf_Area(nodecount,leafcn)
               PDAMLN(K,L) = ranyt * eve_lfstress * calbrt(35)		  &
	            	     * CNReduction_to_Leaf_Area(nodecount,leafcn) 
               PDWMLD(K,L) = PDAMLD(K,L) * Sp_leaf_weight(tday)		  &
					                * CALBRT(36) 
               PDWMLN(K,L) = PDAMLN(K,L) * Sp_leaf_weight(tnyt)		  &
					                * CALBRT(36) 
               SPDWLD = SPDWLD + PDWMLD(K,L)
               SPDWLN = SPDWLN + PDWMLN(K,L)
	      else
               PDAMLD(K,L) = 0.0
               PDAMLN(K,L) = 0.0
               PDWMLD(K,L) = 0.0
               PDWMLN(K,L) = 0.0
	      endif
         enddo
      enddo
      return																     
      end


      SUBROUTINE FBLeafArea

      use common_block

      do k=1,nvbrch
         nbrch = nfbr(k)
         do l=1,nbrch
            nnid = nnod(k,l)
            do m=1,nnid

! ***	Calculate the leaf age

               FRLfAge(k,l,m) = FRLfAge(k,l,m)						   &
			                 +Duration_leaf_expansion(tavg)*calbrt(40)
               if(FRLfAge(k,l,m).lt.1.0) then
	            agefrlf(k,l,m) = agefrlf(k,l,m) + 1.0
			    nodecount = numpfn + l
	            cntmslfcn(nodecount) = cntmslfcn(nodecount) + 1.0
	            summslfcn(nodecount) = summslfcn(nodecount) + leafcn

! *** Initiated once the new node is produced (ADDMSNOD routine)

                  if(FRLfArea(k,l,m).le.0.001) then
				 call Potential_Fruiting_Branch_Leaf_Area(k,l,m)
	             avelfcn=summslfcn(nodecount-1)/cntmslfcn(nodecount-1)
                   LArea(k,l,m) = FRLfArea(k,l,m) * calbrt(43)		    &
				        * CNReduction_to_Leaf_Area(nodecount,avelfcn)
	             raday = 1.0
	             ranyt = 1.0
                  else

! *** Calculate the current rate of leaf expansion during the day and night

                     raday = (Rate_Leaf_expansion(tday) 		  &
	 		                + Rate_of_reduction(tday)			  &
                            * agefrlf(k,l,m))
                     ranyt = (Rate_Leaf_expansion(tnyt) 		  &
			                + Rate_of_reduction(tnyt)			  &
                            * agefrlf(k,l,m))
                     if(raday.LT.0.003) raday = 0.003
                     if(ranyt.LT.0.003) ranyt = 0.003
                  endif
 	   
! *** Calculate the current leaf area 

                  raday = raday * LArea(k,l,m)
                  ranyt = ranyt * LArea(k,l,m)
			    
! *** Calculate the corresponding daily potential change in leaf weight

                  PDADAY(K,L,M) = raday * day_lfstress * calbrt(35)	   &
				        * CNReduction_to_Leaf_Area(nodecount,leafcn)
                  PDANYT(K,L,M) = ranyt * eve_lfstress * calbrt(35)	   &
				        * CNReduction_to_Leaf_Area(nodecount,leafcn)
                  PDWFLD(K,L,M) = PDADAY(K,L,M) * Sp_leaf_weight(tday) &
							    * CALBRT(36) 
                  PDWFLN(K,L,M) = PDANYT(K,L,M) * Sp_leaf_weight(tnyt) &
							    * CALBRT(36) 
                  SPDWLD = SPDWLD + PDWFLD(K,L,M)
                  SPDWLN = SPDWLN + PDWFLN(K,L,M)
	         else
                  PDADAY(K,L,M) = 0.0
                  PDANYT(K,L,M) = 0.0
                  PDWFLD(K,L,M) = 0.0
                  PDWFLN(K,L,M) = 0.0
               endif
            enddo
         enddo
      enddo
      return
      end


      FUNCTION TempFact_leaf(tmp)
! ************************************************************
! *** THIS FUNCTION Calculate the temperature reduction 	   *
! *** factor (assume max area = 12.65 sq cm at 27oC)		   *
! ************************************************************

      if(tmp.lt.17.0) tmp = 17.0
      if(tmp.gt.36.0) tmp = 36.0
	tempfact_leaf = (-18.599 + tmp*(2.186-0.0381*tmp))/12.65
	if(tempfact_leaf.lt.0.01) tempfact_leaf = 0.01
	if(tempfact_leaf.gt.1.0) tempfact_leaf = 1.0
      RETURN
      END
	

      FUNCTION Duration_leaf_expansion(tmp)
! ************************************************************
! *** THIS FUNCTION Calculates the daily increment of leaf   *
! *** expansion duration			rsquare = 0.95  		   *
! ************************************************************
 
      if(tmp.lt.17.0) tmp = 17.0
      if(tmp.gt.36.0) tmp = 36.0
      Duration_leaf_expansion = -0.09365+tmp*(0.01070-0.0001697*tmp)
      if(Duration_leaf_expansion.lt.0.0) Duration_leaf_expansion = 0.0
      if(Duration_leaf_expansion.gt.1.0) Duration_leaf_expansion = 1.0
      RETURN
      END


      FUNCTION Rate_Leaf_expansion(tmp)
! ************************************************************
! *** THIS FUNCTION Calculate the daily increment of leaf    *
! *** area expansion (RLER)	    	rsquare = 0.95 		   *
! ************************************************************

      if(tmp.lt.17.0) tmp = 17.0
      if(tmp.gt.36.0) tmp = 36.0
      Rate_Leaf_expansion = -0.03390 + 0.02041*tmp
      RETURN
      END


      FUNCTION Rate_of_reduction(tmp)
! ************************************************************
! *** THIS FUNCTION Calculate the daily rate of reduction    *
! *** of leaf area expansion (ror)      	rsquare = 0.98	   *
! ************************************************************

      if(tmp.lt.17.0) tmp = 17.0
      if(tmp.gt.36.0) tmp = 36.0
      Rate_of_reduction = 0.01341 - 0.001879*tmp
      RETURN
      END


      FUNCTION Sp_leaf_weight(tmp)
! ************************************************************
! *               CALCULATE Specific leaf weight             *
! ************************************************************
  
      if(tmp.lt.17.0) tmp = 17.0
      if(tmp.gt.36.0) tmp = 36.0
      Sp_leaf_weight = 0.509 + 313.122 * exp(-0.377 * tmp)

! *** Original GOSSYM function

!      if(tmp.le.13.5) tmp = 13.5
!      Sp_leaf_weight = 1.0/(-0.62142855 + tmp*(0.1093651 
!     .                 - tmp*0.00137566))
      RETURN
      END


      SUBROUTINE  Potential_PreFruiting_Leaf_Area(i)
! ************************************************************
! * CALCULATE potential prefruiting leaf area; rsquare =0.91 *
! ************************************************************
 
      use common_block
  
      if(tavg.lt.17.0) tavg = 17.0
      if(tavg.gt.36.0) tavg = 36.0
      PfLfArea(i) = (6.061+1.8069*i) * TempFact_leaf(tavg)/100.0
      if(PfLfArea(i).LT.0.045) PfLfArea(i) = 0.045
      RETURN
      END


      SUBROUTINE  Potential_MainStem_Leaf_Area(k,l)
! ************************************************************
! *  CALCULATE potential main stem leaf area; rsquare = 0.95 *
! ************************************************************
 
      use common_block
  
      if(tavg.lt.17.0) tavg = 17.0
      if(tavg.gt.36.0) tavg = 36.0
      XMSLfArea(k,l) = (18.3812-0.523*l) * TempFact_leaf(tavg)/100.0
      if(XMSLfArea(k,l).LT.0.03) XMSLfArea(k,l) = 0.03
      RETURN
      END


      SUBROUTINE  Potential_Fruiting_Branch_Leaf_Area(k,l,m)
! ************************************************************
! * CALCULATE potential fruiting branch leaf area; r2 = 0.98 *
! ************************************************************
 
      use common_block
  
      if(tavg.lt.17.0) tavg = 17.0
      if(tavg.gt.36.0) tavg = 36.0
      FRLfArea(k,l,m) = (13.457-1.179*m) * TempFact_leaf(tavg)/100.0
      if(FRLfArea(k,l,m).LT.0.045)FRLfArea(k,l,m) = 0.045
      RETURN
      END


      FUNCTION CNReduction_to_Leaf_Area(nodenumber,avelfcn)
! ************************************************************
! *** THIS FUNCTION Calculate the N-reduction to leaf area   *
! ************************************************************

 	faclfn = avelfcn / 0.035
	if(faclfn.gt.1.0) faclfn = 1.0
	redlfcn = -3.5288 + faclfn * (8.0488-3.5200*faclfn)
 	if(redlfcn.lt.0.2) redlfcn = 0.2
	if(redlfcn.gt.1.0) redlfcn = 1.0
	if(nodenumber.eq.1) then
	   CNReduction_to_Leaf_Area = 1.0 
	else
	   CNReduction_to_Leaf_Area = redlfcn
	endif
      RETURN
      END


      SUBROUTINE RUTGRO(KALL)
! ********************************************************************
! *   THIS SUBROUTINE CALCULATES THE GROWTH (IN TERMS OF DRY MATTER) *
! * OF ROOTS IN EACH CELL FOR THE DAY.  FIRST, THE POTENTIAL GROWTH  *
! * PDWRT  FOR THE EXISTING SOIL WATER POTENTIAL (PSIS) AND          *
! * TEMPERATURE (TSOILD & TSOILN) IS CALCULATED FOR EACH SOIL CELL,  *
! * BASED ON THE WEIGHT OF ROOTS CAPABLE OF GROWTH IN EACH CELL      *
! * RTWTCG, THEN THE ACTUAL GROWTH IS DETERMINED, BASED ON THE       *
! * CARBOHYDRATE SUPPLY FOR ROOT GROWTH AND THE POTENTIAL GROWTH FOR *
! * THE CELL. THE ACTUAL GROWTH OCCURING FOR A GIVEN CELL MAY OCCUR  *
! * WITHIN THE CELL OR IN THE CELLS TO THE RIGHT OR LEFT & BELOW.    *
! *    GROWTH IN THE 4 AVAILABLE CELLS IS BASED ON RELATIVE WATER    *
! * POTENTIALS OF THE FOUR, WITH A HEAVIER WEIGHTING GIVEN TO        *
! * DOWNWARD GROWTH.                                                 *
! *    THIS SUBROUTINE DRAWS HEAVILY ON THE IDEAS AND THEORIES OF    *
! * DR. M.G. HUCK, USDA-ARS, AUBURN, AL.  THIS IS ESPECIALLY TRUE    *
! * WITH REGARDS TO SLOUGHING.  C.F. 'A MODEL FOR SIMULATING ROOT    *
! * GROWTH AND WATER UPTAKE ', M.G. HUCK, F.W.T. PENNING DE VRIES,   *
! * AND M.G. KEIZER.  IN PRESS.                                      *
! ********************************************************************
 
      use common_block

      IF(KALL.NE.1) THEN

! *** calculate potential root growth   kit boone

         call potrtgro

      ELSE

! *** THE SECOND TIME RUTGRO IS CALLED FROM GROWTH IT RETURNS HERE 
! *** AND DISTRUBUTES THE ACTUAL (ADJUSTED) CARBOHYDRATE AVAILABLE 
! *** FOR ROOT GROWTH. 

          call actrtgro

! *** calculate average soil water potential (PSIAVG)  kit boone

          call soilpsi

      ENDIF
      RETURN
      END


      SUBROUTINE POTRTGRO
! *******************************************************************
! * THE POTENTIAL ROOT GROWTH (PDWRT) FOR THE EXISTING SOIL WATER   *
! * POTENTIAL (PSIS) AND SOIL TEMPERATURE (TSOILD AND TSOILN) IS    *
! * CALCULATED FOR EACH SOIL CELL, BASED ON THE WEIGHT OF ROOTS     *
! * CAPABLE OF GROWTH IN EACH  CELL (RTWTCG).                       *
! *******************************************************************

      use common_block

! *** variable initialization.
  
      DO I=1,NL
         DO J=1,NK
             DWRT(I,J)=0.
             PDWRT(I,J)=0.
             RTWTCG(I,J)=0.
         enddo
      enddo
      SPDWRT = 0.
      PDWCUTL = 0.0
      PDWCUTR = 0.0
      DO L = 1,LR

! *** KL and KR represent initial rooted cell to the left and right
! *** of the plant in a given layer.  Loops for a layer run from         
! *** initial to terminal rooted cell. 

         KL = KLL(L)
         KR = KRL(L)
   
         DO K = KL,KR
            RTWTCG(L,K) = RTWT(L,K,1) + RTWT(L,K,2)
         enddo
      enddo
         
! *** calculate root impedance factors  

      CALL RIMPED

      DO L = 1, LR
         LP1 = L + 1 - (L/NL)

! *** calculate root growth exponent (ROOTXP) as a function 
! *** of soil temperature 

         TSDL = TSOILD(L)
         TSNL = TSOILN(L)
         IF(TSDL.GT.30.)TSDL = 30.
         IF(TSNL.GT.30.)TSNL = 30.
         IF(TSDL.LT.13.5)TSDL = 13.5
         IF(TSNL.LT.13.5)TSNL = 13.5
         RUTDAY =(-0.2120865+0.016079*TSDL)*DAYTYM
         RUTNYT =(-0.2120865+0.016079*TSNL)*NYTTYM
         ROOTXP= RUTDAY + RUTNYT 

         KL = KLL(L)
         KR = KRL(L)
         DO K = KL, KR
            IF(K.LE.(KULCLF+3)) THEN
               KMOD=(KULCLF+3)-(K-KLL(L))
            ELSE
               KMOD=K
            ENDIF
            PDWRT(L,KMOD)=RTWTCG(L,KMOD)*ROOTXP

! *** POTENTIAL ROOT GROWTH IS REDUCED IN CELLS WITH LOW NITROGEN CONTENT
! *** calculate reduction in root growth due to nitrogen stress.

! *** calculate growth reduction on the left side of the plant  ***

            IF(KMOD.LE.(KULCLF+3)) THEN
               IF(VNO3C(L,KMOD).LT.CALBRT(20).OR.			&
     			             PSIS(L,KMOD).LT.-1.0) THEN
                  FULGRO = PDWRT(L,KMOD)
                  PDWRT(L,KMOD) = PDWRT(L,KMOD) * CALBRT(21)
                  PDWCUTL = PDWCUTL + (FULGRO - PDWRT(L,KMOD))
               ELSE
                  PDWRT(L,KMOD) = PDWRT(L,KMOD) + PDWCUTL*0.10
                  PDWCUTL = PDWCUTL * 0.90
               ENDIF
            ELSE

! *** calculate growth reduction on the right side of the plant  ***

               IF(VNO3C(L,KMOD).LT.CALBRT(20).OR.			 &
                    		  PSIS(L,KMOD).LT.-1.0) THEN
                  FULGRO = PDWRT(L,KMOD)
                  PDWRT(L,KMOD) = PDWRT(L,KMOD) * CALBRT(21)
                  PDWCUTR = PDWCUTR + (FULGRO - PDWRT(L,KMOD))
               ELSE
                  PDWRT(L,KMOD) = PDWRT(L,KMOD) + PDWCUTR*0.10
                  PDWCUTR = PDWCUTR * 0.90
               ENDIF
            ENDIF

            KP1 = KMOD + 1 - (KMOD/NK)
            KM1 = KMOD - 1 + (1/KMOD)

! *** calculate root growth reduction due to soil impedance

            TEST = RTIMPD(L,KMOD)
            IF(TEST.GE.RTIMPD(L,KM1)) TEST = RTIMPD(L,KM1)
            IF(TEST.GE.RTIMPD(L,KP1)) TEST = RTIMPD(L,KP1)
            IF(TEST.GE.RTIMPD(LP1,KMOD)) TEST = RTIMPD(LP1,KMOD)
            RTPCT= (104.6 - 3.53*TEST/1.0216)*.01
            IF(RTPCT.GT.1.0) RTPCT = 1.0
            IF(RTPCT.LT.0.5) RTPCT = 0.5
            PDWRT(L,KMOD) = PDWRT(L,KMOD)*RTPCT
            SPDWRT = SPDWRT + PDWRT(L,KMOD)
         enddo
      enddo
      return
      end


      SUBROUTINE ACTRTGRO
! *****************************************************************
! *   THIS SUBROUTINE CALCULATES THE ACTUAL ROOT GROWTH, BASED    *
! * ON THE CARBOHYDRATE SUPPLY FOR ROOT GROWTH AND THE POTENTIAL  *
! * GROWTH FOR THE CELL.  THE ACTUAL GROWTH OCCURING FOR A GIVEN  *
! * CELL MAY OCCUR WITHIN THE CELL OR IN THE CELLS TO THE RIGHT   *
! * OR LEFT AND BELOW.  GROWTH IN THE 4 AVAILABLE CELLS IS BASED  *
! * ON RELATIVE WATER POTENTIALS OF THE 4 CELLS, WITH A HEAVIER   *
! * WEIGHTING GIVEN TO DOWNWARD GROWTH.          kit bone         *
! *****************************************************************
 
      use common_block

      EFAC1 = 0.0
      EFACL = 0.0
      EFACR = 0.0
      EFACD = 0.0
      SRWP = 0.0

! *** calculate the root growth correction factor.

      IF(SPDWRT.GT.0.0) THEN 
         RGCF = RCH2O / SPDWRT
      ELSE
         RGCF = 0.0
      ENDIF 
      DO L = 1, LR
         KL = KLL(L)
         KR = KRL(L)
         DO K = KL, KR
            IF(K.LE.(KULCLF+3)) THEN
               KMOD=(KULCLF+3)-(K-KLL(L))
            ELSE
               KMOD=K
            ENDIF

! *** calculate the root weight distribution in each age class.


            IF(KDAY.GT.5) THEN
               RTWT(L,Kmod,2)=RTWT(L,Kmod,2)+RTP1*RTWT(L,Kmod,1)
               RTWT(L,Kmod,1) = RTWT(L,Kmod,1) * (1.-RTP1)
               IF(KDAY.GT.15) THEN
                  RTWT(L,Kmod,3)=RTWT(L,Kmod,3)+RTP2*RTWT(L,Kmod,2)
                  RTWT(L,Kmod,2) = RTWT(L,Kmod,2) * (1.-RTP2)
               ENDIF
            ENDIF
            DWRT(L,KMOD) = RGCF * PDWRT(L,KMOD)
            DWRT(L,K) = RGCF * PDWRT(L,K)
         enddo
      enddo
      LRT = LR
      NLR = LR
      DO 240 L=1,NLR

! *** define the direction indices and coefficients.
 
         LDC = GEOTR 
         LD1 = L + 1 - L/NL
         KL = KLL(L)
         KR = KRL(L)
         DO 240 K=KL,KR
            IF(K.LE.(KULCLF+3)) THEN
               KMOD=(KULCLF+3)-(K-KLL(L))
            ELSE
               KMOD=K
            ENDIF
            IF(ROOTSV(L,Kmod).GE.THRLN) THEN
  
! *** define the lateral boundaries of root growth.  Assume 
! *** that at the edges the root growth loops back in.

               KR1 = Kmod + 1 - Kmod/NK
               KL1 = Kmod - 1 + 1/Kmod

! *** IRC = 1 = ILC allows looping back of root growth.

               IRC = 1
               ILC = 1

! *** calculate the combined effect of soil impedance 
! *** and moisture content on root growth.

               call rtswrf(l,k,efac1,efacl,efacr,efacd,srwp)

               RTWT(L,KMOD,1)=RTWT(L,KMOD,1)+DWRT(L,KMOD)*(EFAC1/SRWP)
               RTWT(L,KL1,1)=RTWT(L,KL1,1)+DWRT(L,KMOD)*(EFACL/SRWP)
               RTWT(L,KR1,1)=RTWT(L,KR1,1)+DWRT(L,KMOD)*(EFACR/SRWP)
               RTWT(LD1,KMOD,1)=RTWT(LD1,KMOD,1)					&
                               +DWRT(L,KMOD)*(EFACD/SRWP)

! *** Update initial rooted cell if appropriate to do so.  ***

               IF(KMOD.EQ.KLL(L).AND.KL.GT.1) KLL(L)=KLL(L)-1
               IF(KMOD.EQ.KR.AND.KR.LT.NK) KRL(L)=KRL(L)+1
               IF(L.NE.LR.OR.LR.GE.NL) GO TO 240

! *** LRT is updated if K is initial column; not K = 1.  ***

               IF(KMOD.EQ.KL) LRT = LR + 1

! *** If initial column of next layer is to right of or below current cell,  
! *** then set initial column of next layer to value of current cell. 

               IF(KMOD.LE.KLL(L+1)) KLL(L+1) = KMOD

! *** If initial column of next layer is to left of or below current cell,  
! *** then set terminal column of next layer to value of current cell. 

               IF(KMOD.GE.KRL(L+1)) KRL(L+1) = KMOD
               GO TO 240
            ENDIF
            RTWT(L,KMOD,1) = RTWT(L,KMOD,1) + DWRT(L,KMOD)
  240 CONTINUE

! *** Cultivated cells now located at sides of slab.  ***

      DO 280 J=1,5         
         IF(KULDAY(J).NE.DAYNUM) GO TO 280
         DO K=1,KULCLF
            DO I=1,3
               RTWT(1,K,I) = 0.0
            enddo
         enddo
         DO K=KULCRT,NK
            DO I=1,3
               RTWT(1,K,I) = 0.0
            enddo
         enddo
 280  continue

      LR = LRT
      ROOTS = 0.
      WTSLFD = 0.

! *** calculate the amount of sloughed roots.

      DO L = 1, LR
         KL = KLL(L)
         KR = KRL(L)
         DO K = KL, KR
            WTBSLF = RTWT(L,K,2)
            RTWT(L,K,2) = WTBSLF*(1. - SLF)
            WTSLFD = WTSLFD + (WTBSLF-RTWT(L,K,2))
            ROOTSV(L,K) = RTWT(L,K,1)+RTWT(L,K,2)+RTWT(L,K,3)
            ROOTS = ROOTS + ROOTSV(L,K)
         enddo
      enddo

      ROOTWT = ROOTS * POPFAC
      RUTOFF = RUTOFF + (WTSLFD * POPFAC)

! *** ADJUST ROOTN AND CALCULATE NITROGEN LOSS FROM DEAD ROOTS

      ROOTN  = ROOTN  - (WTSLFD * POPFAC) * ROOTCN
      NLOSS  = NLOSS  + (WTSLFD * POPFAC) * ROOTCN
      PIXLOS = (WTSLFD * POPFAC) * PIXCON
      PIXPLT = PIXPLT - PIXLOS
      return
      end


      subroutine rtswrf(l,k,efac1,efacl,efacr,efacd,srwp)

! *****************************************************************
! *** calculate the combined effect of soil impedance and moisture 
! *** content on root growth  
! *****************************************************************

      use common_block

      STR1 = (104.6 - 3.53*RTIMPD(L,Kmod)/1.0216)*.01
      IF(STR1.GT.1.)  STR1 = 1.
      IF(STR1.LT.0.)  STR1 = 0.
      STRL = (104.6 - 3.53*RTIMPD(L,KL1)/1.0216)*.01
      IF(STRL.GT.1.)  STRL = 1.
      IF(STRL.LT.0.)  STRL = 0.
      STRR = (104.6 - 3.53*RTIMPD(L,KR1)/1.0216)*.01
      IF(STRR.GT.1.)  STRR = 1.
      IF(STRR.LT.0.)  STRR = 0.
      STRD = (104.6 - 3.53*RTIMPD(LD1,Kmod)/1.0216)*.01
      IF(STRD.GT.1.)  STRD = 1.
      IF(STRD.LT.0.)  STRD = 0.

      SWFAC = (1./PSIS(L,Kmod)**3) + (ILC/PSIS(L,KL1)**3) +		 &
              (IRC/PSIS(L,KR1)**3) + (LDC/PSIS(LD1,Kmod)**3)
      WEFAC1 = (1.0/PSIS(L,Kmod)  **3)/SWFAC
      WEFAC1 = (1.0/PSIS(L,K)  **3)/SWFAC
      WEFACL = (ILC/PSIS(L,KL1)**3)/SWFAC
      WEFACR = (IRC/PSIS(L,KR1)**3)/SWFAC
      WEFACD = (LDC/PSIS(LD1,Kmod)**3)/SWFAC
      SRIMPD = STR1 + STRL + STRR + STRD
      IF(SRIMPD.GT.0.0) THEN
         RIMP1 = STR1/SRIMPD
         RIMPL = STRL/SRIMPD
         RIMPR = STRR/SRIMPD
         RIMPD = STRD/SRIMPD
      ELSE
         RIMP1 = 1.
         RIMPL = 1.
         RIMPR = 1.
         RIMPD = 1.
      ENDIF
 
      EFAC1 = WEFAC1 * RIMP1
      EFACL = WEFACL * RIMPL
      EFACR = WEFACR * RIMPR
      EFACD = WEFACD * RIMPD
      SRWP = EFAC1 + EFACL + EFACR + EFACD
      return
      end

! *** inserted on june 19, 1996
! *** end of rewritten rutgro code


       subroutine soilpsi
! **************************************************************
! * Calculate average soil water potential.    kit june 1996   *
! **************************************************************

      use common_block

      DIMENSION NCELLS(9), TVH2OC(9), ARDX(9), AVGH2O(9)
      DIMENSION PAWRCU(40,20), SORT(800,3), INDX(800)

      PSITOT = 0.
      PSINUM = 0.
      SUMH2O = 0.
      NUMCWR = 0

! *** PAWRCU IS PLANT AVAILABLE WATER TIMES ROOT WEIGHT CAPABLE OF UPTAKE
! *** NUMCWR IS NUMBER OF CELLS WITH ROOTS

      J = 0
      AVAILN = 0.0
      DO 360 L = 1, LR
         KL = KLL(L)
         KR = KRL(L)
         DO 360 K = KL, KR
            AVAILN = AVAILN + VNO3C(L,K) + VNH4C(L,K)
            PAWRCU(L,K) = AMAX1(0.,VH2OC(L,K)-THTR(L))*UPF(L,K)
            NUMCWR = NUMCWR + 1
            IF (J .LE. 600) THEN
               J = J + 1
               SORT(J,1) = VH2OC(L,K)-THTR(L)
               SORT(J,2) = AMAX1(0.0,ZUPT(L,K))
               SORT(J,3) = PAWRCU(L,K)
               INDX(J) = ((L-1)*NK) + K
            ENDIF
  360 CONTINUE

      N = MIN(J,600)
 
! *** ORDER SORT IN DESCENDING ORDER BY PAWRCU VALUES 

      CALL SORTIT(N,SORT,INDX)

! *** CALCULATE RUNNING AVERAGE OF EP FOR THE LAST FIVE DAYS. THIS 
! *** IS DONE TO SMOOTH THE SHARP FLUCTUATIONS IN DAILY EP. MOVE ALL 
! *** VALUES DOWN ONE SLOT, PUT NEW VALUE IN FIRST SLOT, AVERAGE 

      DO 1 J =5,2,-1
         EPAVG(J) = EPAVG(J-1)
    1 CONTINUE
      EPAVG(1) = EP
      WEP = 0.
      DO 2 J =1,5
         WEP = WEP + EPAVG(J)
    2 CONTINUE
      WEP = WEP/5.

! *** CALCULATE THE AVERAGE VOLUMETRIC WATER CONTENT OF THE ROOTED 
! *** CELLS IN EACH SOIL HORIZON. ICUR IS THE INDX OF THE CURRENT 
! *** HORIZON. DETERMINE HOW MANY CELLS WITH ROOTS WILL TAKE TO  
! *** SATISFY THE TRANSPIRATION DEMAND (WEP). TWAT IS PAW OF A 
! *** CELL (m) SUMH2O IS THE ACCUMULATED TWAT UP TO THAT ITERATION. 

      DO 365 I = 1, 9
         NCELLS(I) = 0
         TVH2OC(I) = 0.0
 365  CONTINUE

      DO 350 J = 1, NUMCWR

! *** CALCULATE INDICES OF THIS CELL

         IROW = (INDX(J) / NK) + 1
         ICOL = MOD(INDX(J), NK)      
         if(icol.eq.0)icol = 20

! *** FIND HORIZON INDX OF THIS CELL

         IF (PSIS(IROW,ICOL) .GT. -15) THEN
            ICUR = LYRDPH(IROW)
            ARDX(ICUR) = ARDRCN(IROW)
            TWAT = (SORT(J,2)*ACELLDW)/(WCELL*.1*NK)
            SUMH2O = SUMH2O + TWAT
            PSINUM = PSINUM + 1.0
            NCELLS(ICUR) = NCELLS(ICUR) + 1
            TVH2OC(ICUR) = VH2OC(IROW,ICOL) + TVH2OC(ICUR)
            IF(SUMH2O.GT.WEP)GO TO 370
         ENDIF
  350 CONTINUE
  370 CONTINUE
    
! *** CALCULATE AVERAGE WATER CONTENT OF EACH HORIZON WITH ROOTS

      DO 385 I = 1,9
         IF (NCELLS(I) .GT. 0) AVGH2O(I) = TVH2OC(I)/NCELLS(I)
 385  CONTINUE
     
! *** CALCULATE PSI FOR AVERAGE WATER CONTENT, CALCULATE WEIGHTED 
! *** AVERAGE OF PSI VALUES

      PSIAVG = 0.0
      DO 395 N = 1,9
         IF (NCELLS(N) .GT. 0) THEN
            TEMP1G = (AVGH2O(N)-AIRDR(N))/(FCININ(N)-AIRDR(N))
            PSI= PSISFC * TEMP1G**ARDX(N)
            IF (PSI .LT. -15.0) PSI = -15.0
            PSIAVG = PSIAVG + (PSI*(NCELLS(N)/PSINUM))
         ENDIF
 395  CONTINUE
      PSIPCNT = PSINUM/NUMCWR
   
      RETURN
      END


      SUBROUTINE SORTIT(N,SORT,INDX)

! *** SORTS THE 2-D ARRAY SORT IN DESCENDING ORDER BY THE 3RD COL
! *** ALGORITHM FOR THE HEAP SORT USED IS IN NUMERICAL RECIPES,
! *** BY PRESS, FLANNERY, TEUKOLSKY, VETTERLING,  PP. 231

      REAL SORT(800,3)
      INTEGER INDX(800)

         L = N/2+1
         IR = N
  10     CONTINUE
           IF (L .GT. 1) THEN
              L=L-1
              SORT3 = SORT(L,3)
              SORT2 = SORT(L,2)
              SORT1 = SORT(L,1)
              INDX1 = INDX(L)
           ELSE
              SORT3 = SORT(IR,3)
              SORT2 = SORT(IR,2)
              SORT1 = SORT(IR,1)
              INDX1 = INDX(IR)
              SORT(IR,3) = SORT(1,3)
              SORT(IR,2) = SORT(1,2)
              SORT(IR,1) = SORT(1,1)
              INDX(IR) = INDX(1)
              IR = IR-1
              IF (IR .EQ. 1) THEN
                  SORT(IR,3) = SORT3
                  SORT(IR,2) = SORT2
                  SORT(IR,1) = SORT1
                  INDX(IR) = INDX1
                  RETURN
              ENDIF
            ENDIF
           I = L
           J = L + L
  20  IF (J .LE. IR) THEN
         IF (J .LT. IR) THEN
             IF (SORT(J,3) .GT. SORT(J+1,3)) J = J + 1
         ENDIF
         IF (SORT3 .GT. SORT(J,3)) THEN
             SORT(I,3) = SORT(J,3)
             SORT(I,2) = SORT(J,2)
             SORT(I,1) = SORT(J,1)
             INDX(I) = INDX(J)
             I = J
             J = J + J
         ELSE
             J = IR + 1
         ENDIF
         GO TO 20
      ENDIF
      SORT(I,3) = SORT3
      SORT(I,2) = SORT2
      SORT(I,1) = SORT1
      INDX(I) = INDX1
      
      GO TO 10

      END

! *** inserted 3/31/99 kit

      SUBROUTINE RIMPED
! ********************************************************************
! *    THIS SUBROUTINE CALCULATES ROOT IMPEDENCE BASED UPON THE BULK *
! * DENSITY AND WATER CONTENT. THIS IS BASED UPON DATA FROM ARTICLES *
! * BY R.B. CAMPBELL, D.C.REICOSKY AND C.W.DOTY J.OF SOIL AND WATER  *
! * CONS. 29:220-224,1974 AND H.M.TAYLOR AND H.R.GARDNER. SOIL SCI.  *
! * SOIL SCI.96:153-156,1963.                                        *
! *    A LINEAR TABLE LOOK-UP PROCEDURE IS USED. ASSUME ALL          *
! * CURVES ARE READ AT THE SAME BD.                                  *
! ********************************************************************
 
      use common_block

      DO 99 L = 1,NL
         J = LYRDPH(L)
         JJ = 1
   26    IF(BD(J).GT.TSTBD(1,JJ)) THEN
            JJ = JJ+1
            IF(JJ.LT.INRIM) GO TO 26
         ENDIF
         IF(JJ.GT.INRIM) JJ=INRIM
         DO 98 K = 1,NK
            TEST1=VH2OC(L,K)/BD(J)
            IK = 1
   32       IF(TEST1.GT.GH2OC(IK)) THEN
               IK = IK+1
               IF(IK.LT.NCURVE) GO TO 32
            ENDIF
 
! *** SOIL CELL H2O LESS THAN TEST H2O

            IF(JJ.EQ.1) THEN
               IF((IK.EQ.1).OR.(TEST1.EQ.GH2OC(IK))) THEN
                  RTIMPD(L,K)=TSTIMP(IK,JJ)
                  GO TO 98
               ENDIF
            ELSE
               IF((IK.EQ.1).OR.(TEST1.EQ.GH2OC(IK))) THEN
                  RTIMPD(L,K)=TSTIMP(IK,JJ-1)-(TSTIMP(IK,JJ-1)-		   &
                              TSTIMP(IK,JJ))*((TSTBD(IK,JJ-1)-BD(J))/  &
                             (TSTBD(IK,JJ-1)-TSTBD(IK,JJ)))
                  GO TO 98
               ENDIF
            ENDIF
 
! *** CALCULATE SOIL STRENGTH FOR VALUES OF BD LESS THAN TABLE VALUES
 
            IF(JJ.EQ.1) THEN										   
               RTIMPD(L,K)=TSTIMP(IK-1,JJ)-(TSTIMP(IK-1,JJ) 		   &
            			   -TSTIMP(IK,JJ))*((TEST1-GH2OC(IK-1))		   &
                           /(GH2OC(IK)-GH2OC(IK-1)))
            ELSE
 
! *** FOR VALUES OF BD AND H2O BETWEEN TABLE VALUES
 
           TEMP1R=TSTIMP(IK,JJ-1)-(TSTIMP(IK,JJ-1)-TSTIMP(IK,JJ))*		&
              ((TSTBD(IK,JJ-1)-BD(J))/(TSTBD(IK,JJ-1)-TSTBD(IK,JJ)))
           TEMP2=TSTIMP(IK-1,JJ-1)-(TSTIMP(IK-1,JJ-1)-TSTIMP(IK-1,JJ))*	&
           ((TSTBD(IK-1,JJ-1)-BD(J))/(TSTBD(IK-1,JJ-1)-TSTBD(IK-1,JJ)))
           RTIMPD(L,K)=TEMP2+(TEMP1R-TEMP2)*((TEST1-GH2OC(IK-1))/		&
                       (GH2OC(IK)-GH2OC(IK-1)))
            ENDIF
   98    CONTINUE
   99 CONTINUE
      JJ = 1
 126  IF(BDI.GT.TSTBD(1,JJ)) THEN
         JJ = JJ + 1
         IF(JJ.LT.INRIM) GO TO 126
      ENDIF
      IF(JJ.GT.INRIM) JJ=INRIM
         DO 198 K=1,3
            TEST1 = VH2OC(1,K)/BDI
            IK = 1
 132        IF(TEST1.GT.GH2OC(IK)) THEN
               IK = IK + 1
               IF(IK.LT.NCURVE) GO TO 132
            ENDIF
            IF(IK.GT.NCURVE) IK=NCURVE
            IF(JJ.LE.1) THEN
            IF((IK.EQ.1).OR.(TEST1.EQ.GH2OC(IK))) THEN
               RTIMPD(1,K) = TSTIMP(IK,JJ)
               GO TO 198
            ENDIF
         ELSE
            IF((IK.EQ.1).OR.(TEST1.EQ.GH2OC(IK))) THEN
               RTIMPD(1,K) = TSTIMP(IK,JJ-1) - (TSTIMP(IK,JJ-1)-	  &
                        TSTIMP(IK,JJ)) * ((TSTBD(IK,JJ-1)-BDI)/	      &
                       (TSTBD(IK,JJ-1)-TSTBD(IK,JJ)))
            GO TO 198
         ENDIF
      ENDIF
      IF(JJ.LE.1) THEN
         RTIMPD(1,K) = TSTIMP(IK-1,JJ) - (TSTIMP(IK-1,JJ)-			  &
          TSTIMP(IK,JJ))*((TEST1-GH2OC(IK-1))/(GH2OC(IK)-GH2OC(IK-1)))
      ELSE
         TEMP1R = TSTIMP(IK,JJ-1) - (TSTIMP(IK,JJ-1)-TSTIMP(IK,JJ)) *	 &
                ((TSTBD(IK,JJ-1)-BDI)/(TSTBD(IK,JJ-1)-TSTBD(IK,JJ)))
         TEMP2 = TSTIMP(IK-1,JJ-1)-(TSTIMP(IK-1,JJ-1)-TSTIMP(IK-1,JJ))*	 &
             ((TSTBD(IK-1,JJ-1)-BDI)/(TSTBD(IK-1,JJ-1)-TSTBD(IK-1,JJ)))
         RTIMPD(1,K) = TEMP2 + (TEMP1R-TEMP2) * ((TEST1-GH2OC(IK-1))/	 &
                      (GH2OC(IK)-GH2OC(IK-1)))
      ENDIF
 198  CONTINUE
      RETURN
      END


      SUBROUTINE NITRO
!  ************************************************************
!  *                                                          *
!  *                   NITRO SUBROUTINE                       *
!  *                                                          *
!  ************************************************************
!  *  IN THIS SUBROUTINE, THE MAXIMUM AND MINIMUM N           *
!  *  CONCENTRATIONS FOR THE VARIOUS ORGANS ARE AS REPORTED   *
!  *  BY JONES ET. AL. (1974) DEVELOPMENT OF A NITROGEN       *
!  *  BALANCE FOR COTTON GROWTH MODELS: A FIRST APPROXIMATION *
!  *  CROP SCI. 14:541-546.                                   *
!  ************************************************************
! 
      use common_block
 
      NF = 1.
      SEEDR1 = 0.
      BURADD = 0.
      SEDADD = 0.
      PLANTN = SLEAFN + ROOTN + STEMN + SEEDN + BURRN
      NV = 1.
      NR = 1.
 
! F2 = RESERVE NITROGEN AVAILABILITY COEFFICIENT.
! THE MIN N CONC. IN LEAVES IS .015, IN STEMS IS .009, IN STEMS & ROOTS
! IS .009, IN BURRS IS .006 ACCORD. TO JONES & HESKETH 73

      LEAFRS = (SLEAFN- .015 * LEAFWT) * F2
      STEMRS = (STEMN - .009 * STEMWT) * F2
      ROOTRS = (ROOTN - .009 * ROOTWT) * F2
 
! RESN IS RESERVE NITROGEN. UNITS ARE G/PLANT
! NPOOL IS TOTAL NITROGEN AVAILABLE FOR GROWTH.
 
      RESN =LEAFRS + STEMRS + ROOTRS
      IF(RESN.LT.0.) RESN = 0.
      NPOOL = (SUPNO3*0.2258) + (SUPNH4*0.6364) + RESN

! *** conversion factor from NO3 to N is 0.2258; NH4 to N is 0.6364
! THE FOLLOWING REPRESENTS POTENTIAL GROWTH REQUIREMENT
 
      LEAFR1 = 0.035 * PDLEAF
      STEMR1 = 0.020 * PDSTEM
      ROOTR1 = 0.020 * PDROOT
      NPOOL = NPOOL - ROOTR1

      IF(NPOOL.LE.0.0) THEN
        NR = (NPOOL+ROOTR1)/ROOTR1
        NPOOL = 0.0
      ENDIF

      REQV   = LEAFR1 + STEMR1 
 
!  IF THERE IS A GREEN BOLL ON PLANT, CALCULATE SEED & BURR REQS
 
      IF(GBOLWT.GT.0.) THEN
          SEEDR1 = PDBOLL * .045 * .416
          BURMIN = PDBOLL * .006 * .278
          BURR1  = PDBOLL * .014 * .278
          BOLL1  = BURMIN + SEEDR1
          REQ1   = BOLL1  + BURR1
          BURADD = BURR1  + BURMIN
          SEDADD = SEEDR1
      ENDIF
 
      IF((REQV+REQ1).GT.NPOOL) THEN
         IF(GBOLWT.GT.0.) THEN
            IF((REQV+REQ1-BURR1).LE.NPOOL) THEN
               BURADD = NPOOL-REQV-REQ1+BURR1
            ELSE
               IF(BOLL1.LE.NPOOL) THEN
                  NV = (NPOOL-BOLL1) / REQV
                  IF(NV.GT.1.) NV = 1.
                  BURADD = BURMIN
               ELSE
                  NV = 0.
                  BURADD = NPOOL*BURMIN / (BURMIN+SEEDR1)
                  SEDADD = NPOOL*SEEDR1 / (BURMIN+SEEDR1)
                  NF = NPOOL/BOLL1
                  IF(NF.GT.1.) NF = 1.
               ENDIF
            ENDIF
         ELSE
            NV = NPOOL / REQV
            IF(NV.GT.1.) NV = 1.
         ENDIF
      ENDIF

! VEGETATIVE GROWTH SECTION

      SLEAFN = SLEAFN + LEAFR1 * NV
      STEMN  = STEMN  + STEMR1 * NV
      ROOTN  = ROOTN  + ROOTR1 * NR

! BOLL GROWTH SECTION

      BURRN = BURRN + BURADD
      SEEDN = SEEDN + SEDADD
      PLTN  = SLEAFN + STEMN + ROOTN + BURRN + SEEDN
      XTRAN = (SUPNO3 + SUPNH4) - (PLTN - PLANTN)
 
!  THE PART OF XTRAN SUPPLIED BY RESN IS (PLTN-PLANTN)-(SUPNO3 + SUPNH4)
!  IF THIS  IS +, THEN SOME CAME FROM RESN. XTRAN MUST BE NEGATIVE
!  IF (PLTN-PLANTN) -(SUPNO3+SUPNH4) IS - , THEN ALL CAME FROM SUPNO3
!  SUPNH4 & XTRAN MUST BE +.
 
      VEGWT  = LEAFWT + STEMWT + ROOTWT
      IF(XTRAN.GT.0.) THEN
          SLEAFN = SLEAFN + XTRAN * (LEAFWT/VEGWT)
          STEMN  = STEMN  + XTRAN * (STEMWT/VEGWT)
          ROOTN  = ROOTN  + XTRAN * (ROOTWT/VEGWT)
      ELSE
          IF(RESN.LE.0.) THEN
              XTRAN = 0.0
          ELSE
              SLEAFN = SLEAFN + XTRAN * (LEAFRS/RESN)
              STEMN  = STEMN  + XTRAN * (STEMRS/RESN)
              ROOTN  = ROOTN  + XTRAN * (ROOTRS/RESN)
          ENDIF
      ENDIF
      XTRAN = 0.
      STEMCN = STEMN / STEMWT
 
      IF(LEAFWT.GT.0.0) THEN
          LEAFCN = SLEAFN / LEAFWT
          ROOTCN = ROOTN / ROOTWT
          IF(GBOLWT.GT.0.) THEN
              XXWT = COTXX + GBOLWT
              SEEDCN = SEEDN / (XXWT*.416)
              BURCN  = BURRN / (XXWT*.278)
          ENDIF
          TNO3UP = TNO3UP + SUPNO3
          TNH4UP = TNH4UP + SUPNH4
      ENDIF
      RETURN
      END


      SUBROUTINE MATBAL
!  **************************************************************
!  *                                                            *
!  *            MATERIALS BALANCE SUBROUTINE                    *
!  *                                                            *
!  **************************************************************
!  *  IN THIS SUBROUTINE, AN INVENTORY OF THE CARBOHYDRATE      *
!  *  AND NITROGEN IS MADE FOR DIAGNOSTIC PURPOSES.             *
!  *  CARBOHYDRATES PRODUCED IN PNET ARE ALLOCATED TO GROWING   *
!  *  POINTS IN GROWTH. IF THE DIFFERENCE OF AVAILABLE AND      *
!  *  ALLOCATED IS ZERO, THE CARBOHYDRATE BALANCE IS CORRECT.   *
!  *  IF THE DIFFERENCE IS NEGATIVE, MORE CARBOHYDRATE HAS      *
!  *  BEEN ALLOCATED THAN WHAT WAS PRODUCED. IF POSITIVE, NOT   *
!  *  ALL CARBOHYDRATE PRODUCED HAS BEEN ALLOCATED.             *
!  *  SIMILAR LOGIC IS USED FOR SOIL AND PLANT NITROGEN.        *
!  *  DEBUGING SHOULD BE DONE IF BALANCES ARE EITHER POSITIVE   *
!  *  OR NEGATIVE.                                              *
!  *                                                            *
!  **************************************************************

      use common_block

! CALCULATE CARBOHYDRATE BALANCE IN EACH ORGAN AND IN THE WHOLE
! PLANT. UNITS ARE GMS/PLANT.

! RUTBAL = ROOTS C BALANCE
! SALBAL = STEM AND LEAF C BALANCE
! FRUBAL = FRUIT C BALANCE
! CHOBAL = IS WHOLE PLANT C BALANCE

      RUTBAL = (SROOT + 0.2) - (ROOTWT + RUTOFF)
      SALBAL = (SSTEM + SLEAF + 0.4) - (STEMWT + LEAFWT + LEFABS)
      FRUBAL = (SQUAR + SBOLL) - (SQWT + PQFLR + GBOLWT + GBLOS + COTXX)
      AVAIL  = SPN + WTDAY1
      USED   = ROOTWT + STEMWT + GBOLWT + LEAFWT + SQWT   + COTXX +	  &
               XTRAC  + GBLOS  + LEFABS + PQFLR  + RUTOFF + RESC
      CHOBAL = AVAIL - USED
 
! CALCULATE SOIL NITROGEN BALANCE. UNITS ARE LBS/ACRE
! 
! DAY1SN  = RESIDUAL AND/OR INITIAL SOIL N
! ADDEEDN = ADDED N AS FERTILIZER
! ORGN    = AVAILABLE N IN ORGANIC MATTER
! SOILN   = SOIL N TODAY
! UPTAKEN = N REMOVED BY PLANT FROM THE SOIL
 
      AVAIL = DAY1SN + ADDEDN + ORGN
      USED  = SOILN + UPTAKEN + (CUMNSOK/ROWSP/.011219)
      SNBAL = AVAIL - USED
 
! CALCULATE PLANT NITROGEN BALANCE. UNITS ARE MG N/GRM DRY MATTER.
 
      AVAIL = DAY1PN + TNO3UP + TNH4UP + FOLIARN
      USED  = ROOTN + STEMN + SLEAFN + SEEDN + BURRN + NLOSS
      PNBAL = AVAIL - USED
 
! CALCULATE WATER BALANCE. UNIT ARE MM
 
      AVAIL  = TIH2OC + CUMRAN + SUBIRR
      USED   = TH2O + CUMEP + CUMES + CUMSOK
      H2OBAL = AVAIL - USED
 
! CALCULATE PLANT WEIGHT AND DEAD WEIGHT FOR PLANT

      PLANTW = ROOTWT + STEMWT + GBOLWT + LEAFWT + SQWT + XTRAC +	 &
               COTXX + RESC
      DEAD2DAY = LEFABS + PQFLR + GBLOS + RUTOFF - DEADWT
      IF(DEAD2DAY.LT.0.) DEAD2DAY = 0.
      DEADWT = LEFABS + PQFLR + GBLOS + RUTOFF

      RETURN
      END


      SUBROUTINE PLTMAPS
!  *************************************************************
!  *                                                           *
!  *                     PLTMAP SUBROUTINE                     *
!  *                                                           *
!  *************************************************************
!  *                                                           *
!  * FCODE                          MCODE                      *
!  *                                                           *
!  * 1 SQUARE                       1 SQUARE                   *
!  * 2 GREEN BOLL                   2 GREEN BOLL               *
!  * 3 MATURE BOLL                  3 MATURE BOLL              *
!  * 4 FRUIT ABCISED AS A SQUARE    4 FRUIT ABZ AS A SQUARE    *
!  * 5 FRUIT ABCISED AS A BOLL      5 FRUIT ABZ AS A BOLL      *
!  * 6 NOT USED                     6 ABZ AS PINHEAD OR POLYNA *
!  * 7 YOUNG GREEN BOLL             7 YOUNG GREEN BOLL         *
!  *************************************************************
 
      use common_block
 
      POLOSS = 0.
      AGEADD = 0.
      DSITES = SITEZ-SITES
      SITES  = SITEZ
      SITEZ  = 0
      PINHED = 0.
 
! *** Calculate morphogenetic delays due to N- and C-stress
 
	call N_and_C_Delay

! *** Update age of pre-fruiting node

      DO J=1,NUMPFN
        AGEPFN(J) = AGEPFN(J) + 1
        AVTPFN(J) = (AVTPFN(J)*(AGEPFN(J)-1)+TAVG) / AGEPFN(J)
      ENDDO
      AVTEMP = ((KDAY-1)*AVTEMP + TAVG)/ KDAY

! *** If first square has not occurred, calculate the time from emergence 
! *** to first square, TSQ

	if(isq.le.0)then
	   call Time_to_First_Square
	   sum_fsq_tavg = sum_fsq_tavg + tavg

! *** If first square has not occurred, decide whether to add a PFruit node

         IF(KDAY.LT.IFIX(TSQ)) THEN

	      IF(AGEPFN(NUMPFN).LE.66.) THEN

! *** Increment the age of each node and calculate its running ave.temp.

               IF(NUMPFN.LT.9) THEN
 	            call PreFruiting_Node_Time_Interval
                  IF(AGEPFN(NUMPFN).GE.(pfti+CDLAYV+NDLAY))CALL ADDPFNOD
	         endif
	      endif
	      return
	   endif   

! *** HAVE FIRST SQUARE, INITIALIZE THAT NODE AREA, AVGT, AND LEAF.

         CALL FSTSQTODAY
         idum0 = 0
      ENDIF
 
      PIN = 0.0457*DSITES
      PINHED = PIN

! *** DECIDE WHETHER TO ADD A VEGETATIVE BRANCH.

      IF((NVBRCH.LT.mxvbrch).AND.(INT.LT.0.9)) THEN
 	     call Vegetative_Branch_Time_Interval
	     call Time_to_First_Square
	     dum = tsq
         IF(NVBRCH.EQ.2) DUM=0.0
         IF(AGE(NVBRCH,1,1).GE.(vbti+CDLAYV+NDLAY+DUM))	   &
	        call Add_A_Vegetative_Branch
      ENDIF

      DO 30 K = 1,NVBRCH
 
! *** DECIDE WHETHER TO ADD A FRUITING BRANCH TO THIS VEGETATIVE BRANCH.
 
         NBRCH = NFBR(K)
         IF(NBRCH.LT.mxfbrch) THEN
            VDELAY(K) = VDELAY(K) + ((CDLAYV  + NDLAY )/ PIXDN)
 	        call Main_Stem_Node_Time_Interval(k,nbrch)
            IF(AGE(K,NBRCH,1).GE.(xmsti+VDELAY(K))) CALL ADDMSNOD(K)
         ENDIF

! *** DECIDE WHETHER TO ADD A NODE TO THIS FRUITING BRANCH OF THIS
! *** VEGETATIVE BRANCH. 

         NBRCH = NFBR(K)
         DO L=1,NBRCH
            NNID = NNOD(K,L)
            IF(NNID.LT.mxfsite) THEN
               DELAY(K,L) = DELAY(K,L) + ((CDLAYF+NDLAY) / PIXDN)
 	           call Fruiting_Branch_Node_Time_Interval(k,l,nnid)
               IF(AGE(K,L,NNID).GE.(fbnti+DELAY(K,L)))				&
			  
! *** Add A fruiting branch.

                  call ADDFBNOD(k,l,nnid)
            ENDIF
 
! *** AGE ALL EXISTING NODES AND UPDATE AVERAGE TEMPERATURE OF EACH.
! *** AVGT(K,L,M) IS BOLL RUNNING AVERAGE TEMPERATURE	SINCE BOLLSET
 
         DO 50 M=1,NNID
            AGE(K,L,M) = AGE(K,L,M) + 1.0

!		   Adj_temp = TAVG + Canopy_Temp_Differential(psild)
!            AVGT(K,L,M) = (AVGT(K,L,M)*(AGE(K,L,M)-1) + Adj_temp)/
!     .                     AGE(K,L,M)

		    AVGT(K,L,M) = (AVGT(K,L,M)*(AGE(K,L,M)-1)+TAVG)/AGE(K,L,M)
            AVTNOD = AVGT(K,L,M)
            AGENOD = AGE(K,L,M)

! *** Set temperature used in BLOOM calculation to include canopy temperature
! *** differential; used in the calculation of heatindex.

            Bloom_tavg(daynum) = avtnod

! *** AGE ALL LEAVES

            AGEFAC = AMAX1((1.0-WSTRS),(1.0-NV))*CALBRT(38) 
            LAGE(K,L,M)=LAGE(K,L,M) + 1.0 + AGEFAC
 
! *** CALCULATE BOLL TEMP
 
 	        call Boll_Temperature(k,l,m)

! *** FCODE =          1  2  3  4  5  6  7

            IF(FCODE(K,L,M).GT.0) THEN

! *** Squares are present

               IF(FCODE(K,L,M).EQ.1) THEN
	            call Time_to_Bloom(avtnod)
				if((k.eq.1).and.(l.eq.1).and.(m.eq.1)) then
				   if(iflday.eq.0.and.ifbl.lt.kday) 		   &
	                  sum_fbl_tavg = sum_fbl_tavg + tavg
				endif
 
! *** IF SQUARE IS OLD ENOUGH, MAKE IT A Young GREEN BOLL (FCODE=7).

                  IF(AGENOD.GE.BLOOM) THEN

! *** Estimate boll weight and update square weight 

	               call Update_Boll_and_Square_Weight(k,l,m)

                     IF(GBOLWT.GT.0.0.AND.FBLOOM.LE.1.) 		&
							        	call First_Bloom_Today

! *** Estimate bloom loss due to a rain event during pollination

			       call Loss_due_Pollination(k,l,m)
                  ENDIF

! *** Young green bolls are present

               ELSEIF(FCODE(K,L,M).EQ.7) THEN
 
! *** AGEABZ IS THE AGE AFTER WHICH GREEN BOLLS CANNOT BE ABSCISED

	            call Boll_Safe_Age(k,l,m)

! *** Green bolls are present

               ELSEIF(FCODE(K,L,M).EQ.2) THEN
	            call Time_to_Open_Boll(k,l,m)

! *** IF GREEN BOLL IS OLD ENOUGH, MAKE IT AN OPEN BOLL AND SET FCODE TO 3.

                  IF(AGEBOL(K,L,M).GE.DEHISS(K,L,M)) THEN
	               call Seed_and_Burr_Ncontent(k,l,m)
                     COTXX  = COTXX + BOLWGT(K,L,M)
                     GBOLWT = GBOLWT - BOLWGT(K,L,M)

                     GINP = (50.54 - .6755*BOLTMP(K,L,M)) / 100.0

! *** NOPEN IS NUMBER OF OPEN BOLLS
 
                     NOPEN=NOPEN+FFRUT(K,L,M)
                     IF(NOPEN.GT.0.0) GIN=(GINP+GIN)/NOPEN

! *** YIELD=500 LB. BALES/ACRE OF LINT. K=(453.6 G/LB)*(500 LB./BALE)

                     if(bolwgt(k,l,m).le.0.0001) then
				      yield = yield 
				   else
                         YIELD = YIELD + (GINP * (BOLWGT(K,L,M)*0.75)	   &
                                  * POPPLT/226800.)
				   endif

! *** Estimate fiber length and stength

	               call Fiber_Length_and_Strength(avtnod)
                  ENDIF
               ENDIF
               SITEZ = SITEZ + 1
            ENDIF
 50       CONTINUE
 40     enddo	
 30   CONTINUE

      RETURN
      END


	subroutine N_and_C_Delay
! ****************************************************************
! *** Calculate morphogenetic delays due to N and C stresses   ***
! *** Well watered experiments (AAA), Bruce and Romkens (1965) ***
! ****************************************************************
 
	use common_block

      NDLAY = CALBRT(19) * (1.0 - NV)
      IF (NDLAY.GT.1.0) NDLAY = 1.0
      IF (NDLAY.LT.0.0) NDLAY = 0.0
      CDLAYF =  CALBRT(1) +FSTRES *(CALBRT(2) +CALBRT(3) *FSTRES)
      IF (CDLAYF.GT.1.0) CDLAYF = 1.0
      IF (CDLAYF.LT.0.0) CDLAYF = 0.0
      CDLAYV = 1.0 +FSTRES * CALBRT(4) + FSTRES**2 * CALBRT(5)
      IF (CDLAYV.GT.1.0) CDLAYV = 1.0
      IF (CDLAYV.LT.0.0) CDLAYV = 0.0
      RETURN
      END


	subroutine Time_to_First_Square
! **************************************************************
! *** Calculate the time from emergence to first square, TSQ ***
! *** THE NEW TSQ EQN IS FROM REDDY'S 1988 SPAR EXPERIMENT.  ***
! **************************************************************
 
	use common_block

      IF(AVTEMP.LT.17.0) AVTEMP = 17.0
      IF((VARITY(IVARTY).EQ.'PIMA').OR.(VARITY(IVARTY).EQ.'pima')) THEN
         IF(AVTEMP.GT.32.0) AVTEMP = 32.0
         TSQ = (247.722-16.484*AVTEMP+0.306*AVTEMP**2)*CALBRT(29)
      ELSE
         IF(AVTEMP.GT.35.0) AVTEMP = 35.0
         TSQ = (190.338-11.372*AVTEMP+0.194*AVTEMP**2)*CALBRT(29)
      ENDIF
      RETURN
      END


 	subroutine PreFruiting_Node_Time_Interval
! ************************************************************
! *** Calculate the time interval between prefruiting node ***
! ************************************************************
 
	use common_block

      AT = AVTPFN(NUMPFN)
      IF(AT.lT.17.0) AT = 17.0
      IF((VARITY(IVARTY).EQ.'PIMA').OR.(VARITY(IVARTY).EQ.'pima')) THEN
         IF(AT.GT.33.0) AT = 33.0
         pfti = (57.993 - 3.5871*AT+0.058894*AT**2)*CALBRT(26)
      ELSE
         IF(AT.GT.32.0) AT = 32.
         pfti = (41.205-2.67219*AT+0.0459705*AT**2)*CALBRT(26)
      ENDIF
      RETURN
      END


 	subroutine Vegetative_Branch_Time_Interval
! ************************************************************
! *** Calculate the time interval between prefruiting node ***
! ************************************************************
 
	use common_block

      AT = (AVGT(NVBRCH,1,1)*AGE(NVBRCH,1,1)+TAVG)/(AGE(NVBRCH,1,1)+1)
      IF(AT.LT.17.0) AT = 17.0
      IF((VARITY(IVARTY).EQ.'PIMA').OR.(VARITY(IVARTY).EQ.'pima')) THEN
         IF(AT.GT.37.0) AT = 37.0
         vbti = 17.365 - 0.87003*AT + 0.012265*AT**2
      ELSE
         IF(AT.GT.32.0) AT = 32.0
         vbti = 41.205 - 2.67219*AT + 0.0459705*AT**2
      ENDIF
      RETURN
      END


      SUBROUTINE Add_A_Vegetative_Branch
! *****************************************************************
! *** Add other vegetative branches mimicking that of main stem ***
! *****************************************************************

	use common_block

      NVBRCH = NVBRCH + 1
      FFRUT(NVBRCH,1,1) = 1.
      FCODE(NVBRCH,1,1) = 1
      IF(PINHED.GE.1.) THEN
         PINHED = PINHED - 1.
         MCODE(NVBRCH,1,1) = 6
      ENDIF

! *** Calculate the age, weight and fruiting branch leaf area at unfolding 

      FRLfAge(NVBRCH,1,1) = FRLfAge(NVBRCH,1,1)						  &
			              + Duration_leaf_expansion(tavg)*calbrt(40)
      if(FRLfAge(NVBRCH,1,1).lt.1.0) 								  &
	                   agefrlf(NVBRCH,1,1) = agefrlf(NVBRCH,1,1) + 1.0
 
! *** Assume that at node initiation, the leaf area is negligible
!	Larea(NVBRCH,1,1) = 0.04 
!	LA_WT_CF=(Sp_leaf_weight(tday)+Sp_leaf_weight(tnyt))*CALBRT(36)/2.0   
!	LEAFW(NVBRCH,1,1) = 0.04 * Lf_WT_CF
! *** The above 3 statements needs further refinements

	LEAFW(NVBRCH,1,1) = 0.001
      STEMWT = STEMWT -  LEAFW(NVBRCH,1,1)
      SLEAFN = SLEAFN + (LEAFW(NVBRCH,1,1) * 0.035)
      STEMN  = STEMN  - (LEAFW(NVBRCH,1,1) * 0.035)

! *** Calculate the age, weight and leaf area on main stem node at unfolding

      XMSLfAge(NVBRCH,1) = XMSLfAge(NVBRCH,1) 						   &
		                 + Duration_leaf_expansion(tavg)*calbrt(40) 	   
      if(XMSLfAge(NVBRCH,1).lt.1.0) 								   &
	                  agemslf(NVBRCH,1) = agemslf(NVBRCH,1) + 1.0
 
! *** Assume that at node initiation, the leaf area is negligible
!	MLarea(NVBRCH,1) = 0.04 
!	LA_WT_CF=(Sp_leaf_weight(tday)+Sp_leaf_weight(tnyt))*CALBRT(36)/2.0   
!	MLEAFW(NVBRCH,1) = 0.04 * Lf_WT_CF
! *** The above 3 statements needs further refinements

 	MLEAFW(NVBRCH,1) = 0.001
      STEMWT = STEMWT -  MLEAFW(NVBRCH,1)
      SLEAFN = SLEAFN + (MLEAFW(NVBRCH,1) * 0.035)
      STEMN  = STEMN  - (MLEAFW(NVBRCH,1) * 0.035)
      NFBR(NVBRCH)     = 1
      NNOD(NVBRCH,1)   = 1
      LEFCNT=LEFCNT+1
      LEFSRT(LEFCNT)=NVBRCH*1000+10+1
 
      RETURN
      END


 	subroutine Main_Stem_Node_Time_Interval(k,nbrch)
! ************************************************************
! *** Calculate the time interval between prefruiting node ***
! ************************************************************
 
	use common_block

! *** THE DATA FOR THIS ESTIMATE OF THE EFFECT OF TEMPERATURE ON TIME
! *** BETWEEN NODES ON A FRUITING BRANCH ARE PUBLISHED IN HESKETH, J.D.,
! *** D.N. BAKER, AND W.G. DUNCAN (1972) THE SIMULATION OF GROWTH AND
! *** YIELD IN COTTON II. ENVIRONMENTAL CONTROL OF MORPHOGENESIS. CROP
! *** SCI. 12:436-439 - EXCEPT FOR THE FACTOR .51. THIS FACTOR IS 
! *** REQUIRED TO MATCH THE PHYTOTRON GROWTH RATES TO THE FIELD DATA OF 
! *** BRUCE, R.R. AND ROMKENS. (1965) FRUITING AND GROWTH IN RELATION 
! *** TO SOIL MOISTURE TENSION.AGRON. J. 57:135-140.
 
      AT = (AVGT(K,NBRCH,1)*AGE(K,NBRCH,1)+TAVG)/(AGE(K,NBRCH,1)+1)
      IF(AT.LT.17.0) AT=17.0
      IF(K.EQ.1) THEN
         IF(AT.GT.35.0) AT=35.0
         xmsti = (35.919 - 2.252*AT + 0.0375*AT**2)*CALBRT(28)
      ELSE
         IF(AT.GT.32.) AT=32.
         xmsti = (61.416 - 3.833*AT + 0.0652*AT**2)*CALBRT(27)
      ENDIF
 
! *** THE ABOVE EQUATIONS PROVIDES THE SAME MORPHOGENETIC RATE ON  
! *** SECONDARY VEGETATIVE BRANCHES AS ON FRUITING BRANCHES.

      RETURN
      END
 

 	subroutine Fruiting_Branch_Node_Time_Interval(k,l,nnid)
! ****************************************************************
! *** Calculate the time interval between fruiting branch node ***
! ****************************************************************
 
	use common_block

      AT = (AVGT(K,L,NNID)*AGE(K,L,NNID)+TAVG)/(AGE(K,L,NNID)+1)
      IF(AT.LT.17.) AT=17.
      IF(AT.GT.32.) AT=32.
      fbnti = (61.416 - 3.833*AT + 0.0652*AT**2)*CALBRT(27)
      RETURN
      END


 	subroutine Boll_Temperature(k,l,m)
! **************************************************************
! *** Calculate boll temperature                             ***
! **************************************************************
 
	use common_block

      IF(BOLWGT(K,L,M).GT.0.0) THEN
         DUM=TAVG 
!         DUM=TAVG + Canopy_Temp_Differential(psild)
         IF((LAI.LE.CALBRT(37)).AND.(KDAY.GT.100)) 			   &
                  DUM = DUM * (1.129-0.043*LAI)
         BOLTMP(K,L,M) = (BOLTMP(K,L,M) * AGEBOL(K,L,M)+DUM)/  &
                         (AGEBOL(K,L,M)+1)
         AGEBOL(K,L,M) = AGEBOL(K,L,M) + 1.0 
      ENDIF
      RETURN
      END


      FUNCTION Canopy_Temp_Differential(psild)
! ************************************************************
! *** Calculates the canopy temperature differential as a    *
! *** function of midday leaf water potential  KR Reddy 1997 *
! ***		                  rsquare = 0.51			       *
! ************************************************************

	Canopy_Temp_Differential = -4.7739 - 4.2762*psild
	if(Canopy_Temp_Differential.lt.-5.0) Canopy_Temp_Differential=-5.0
	if(Canopy_Temp_Differential.gt.8.0) Canopy_Temp_Differential = 8.0
      RETURN
      END
	

	subroutine Time_to_Bloom(avtnod)
! **************************************************************
! *** Calculate the time from first square to bloom, BLOOM   ***
! **************************************************************
 
	use common_block

      AT = AVTNOD
      IF(AT.LT.17.0) AT = 17.0
      IF(AT.GT.32.0) AT = 32.0
      IF((VARITY(IVARTY).EQ.'PIMA').OR.(VARITY(IVARTY).EQ.'pima')) THEN
         BLOOM=(188.261-10.704*AT+0.1784*AT**2)*CALBRT(30)
      ELSE
         BLOOM=(252.591-15.321*AT+0.2531*AT**2)*CALBRT(30)
      ENDIF
      RETURN
      END


	subroutine First_Bloom_Today
! *****************************************************
! *** Set the date for first bloom                  ***
! *****************************************************
 
	use common_block

! *** Set the date counter for the calculation of ave. temp to FBL

	FBLOOM = KDAY
	iflday = daynum
 	ave_fbl_tavg = sum_fbl_tavg / float(kday)
	sum_fob_tavg = sum_fbl_tavg
	ifbldae = kday
	ifob = kday-1
	numfob = 1
      RETURN
      END


	subroutine Update_Boll_and_Square_Weight(k,l,m)
! **************************************************************
! *** Estimate boll weight and update square weight          ***
! **************************************************************
 
	use common_block

      BLUM(KDAY) = BLUM(KDAY) + FFRUT(K,L,M)
      BOLWGT(K,L,M) = .31 * SQRWT(K,L,M)
      PQFLR = PQFLR + .69 * SQRWT(K,L,M)
      GBOLWT = GBOLWT + BOLWGT(K,L,M)
      SQWT = SQWT - SQRWT(K,L,M)
      SQRWT(K,L,M) = 0.0
      BOLTMP(K,L,M) = TAVG
      AGEBOL(K,L,M) = 1.0
      FCODE(K,L,M) = 7
      RETURN
      END


	subroutine Loss_due_Pollination(k,l,m)
! **************************************************************
! *** Estimate the number of blooms affected by pollination  ***
! **************************************************************
 
	use common_block

      IF(POLYNA.EQ.0) THEN
         MCODE(K,L,M) = 8
         POLOSS = POLOSS + FFRUT(K,L,M)*.25
      ENDIF
      RETURN
      END


	subroutine Boll_Safe_Age(k,l,m)
! **************************************************************
! *** Calculate the boll safe time from abscission           ***
! **************************************************************
 
	use common_block

      FSTAVG(K,L,M) = (FSTAVG(K,L,M)*AGEBOL(K,L,M)+FSTRES)/		&
                                  (AGEBOL(K,L,M)+1)
      IF((PRPDAY.GT.0).AND.(DAYNUM.GT.PRPDAY)) THEN
          AGEABZ(K,L,M) = 12.82-3.229*PRPKGH-.4987*PRPKGH**2+	&
                              1.723*AVGTSP-.03181*AVGTSP**2+	&
                              0.09911*AVGTSP*PRPKGH
      ELSE
         AGEABZ(K,L,M) = CALBRT(11)+15.0-15.0*FSTAVG(K,L,M)
      ENDIF

      IF(AGEBOL(K,L,M).GE.AGEABZ(K,L,M)) FCODE(K,L,M) = 2
	if(iobday.eq.0.and.ifob.lt.kday) then
	   sum_fob_tavg = sum_fob_tavg + tavg
	   ifob = kday
	   numfob = numfob + 1
	endif

      RETURN
      END


	subroutine Time_to_Open_Bolln(k,l,m)
! **************************************************************
! *** Calculate the time from emergence to open boll, DEHISS ***
! **************************************************************
 
	use common_block

      atn = boltmp(k,l,m)
      if(atn.lt.17.0) atn = 17.0
      if(atn.gt.35.0) atn = 35.0
      if(prpkgh.gt.0.01) then
         dehiss(k,l,m) = (446.51 - 27.65*atn + 0.4807*atn**2-127.03*	&
                          prpkgh+10.28*atn*prpkgh-0.2105*prpkgh*atn**2)
      else
         dehiss(k,l,m) = (1.0/(-0.0058298+0.0009947*atn))*calbrt(31)
      endif
      if(dehiss(k,l,m).gt.70.) dehiss(k,l,m) = 70.

	if(iobday.eq.0.and.ifob.lt.kday) then
	   sum_fob_tavg = sum_fob_tavg + tavg
	   ifob = kday
	   numfob = numfob + 1
	endif

      return
      end


	subroutine Time_to_Open_Boll(k,l,m)
! **************************************************************
! *** Calculate the time from emergence to open boll, DEHISS ***
! **************************************************************
 
	use common_block

      ATN=BOLTMP(K,L,M)
      IF(ATN.GT.35.0) ATN = 35.0
      IF(ATN.LT.17.0) ATN = 17.0
      IF(PRPKGH.GT.0.01) THEN
         DEHISS(K,L,M)=(446.51-27.65*ATN+.4807*ATN**2-127.03*		  &
                       PRPKGH+10.28*ATN*PRPKGH-.2105*PRPKGH*ATN**2)
      ELSEIF((VARITY(IVARTY).EQ.'PIMA').OR.							  &
                      (VARITY(IVARTY).EQ.'pima')) THEN
         DEHISS(K,L,M)=(358.264-19.294*ATN+0.291*ATN**2)*CALBRT(31)
      ELSE
         DEHISS(K,L,M)=(327.396-17.251*ATN+0.255*ATN**2)*CALBRT(31)
      ENDIF
      IF(DEHISS(K,L,M).GT.70.) DEHISS(K,L,M)=70.

	if(iobday.eq.0.and.ifob.lt.kday) then
	   sum_fob_tavg = sum_fob_tavg + tavg
	   ifob = kday
	   numfob = numfob + 1
	endif

      RETURN
      END


	subroutine Seed_and_Burr_Ncontent(k,l,m)
! ****************************************************
! *** Calculate the N content of the burr and seed ***
! ****************************************************
 
	use common_block

      FCODE(K,L,M) = 3
      IF(iobday.eq.0)then
	   iobday = daynum
	   ifob = kday
	   ifobdae = kday
 	   ave_fob_tavg = sum_fob_tavg / float(kday)
	endif

      IF(GBOLWT.GT.0.0) THEN
         SEEDN = SEEDN - BOLWGT(K,L,M)*0.416*SEEDCN
         BURRN = BURRN - BOLWGT(K,L,M)*0.278*BURCN
         NLOSS = NLOSS + BOLWGT(K,L,M)*0.416*SEEDCN
         NLOSS = NLOSS + BOLWGT(K,L,M)*0.278*BURCN
      ENDIF
      RETURN
      END


	subroutine Fiber_Length_and_Strength(avtnod)
! **************************************************************
! *** Calculate the time from first square to bloom, BLOOM   ***
! **************************************************************
 
	use common_block

! *** FS=FIBER STRENGTH : (G / TEX * 1/8 INCH)
! *** FL IS FIBER LENGTH ( INCHES, 2.5% SPUN )
 
      FSX = 56.603 + AVTNOD *(-2.921+0.059*AVTNOD)
      IF(NOPEN.GT.0.01) FS=(FSX+FS)/NOPEN
      FLX=1.219-.0065*AVTNOD
      IF(NOPEN.GT.0.01) FL=(FL+FLX)/NOPEN
      RETURN
      END


      SUBROUTINE ADDPFNOD
! ***************************************************************
! *** Increment prefruiting node count and initiate leaf area ***
! ***************************************************************

	use common_block

      NUMPFN = NUMPFN + 1
      AVTPFN(NUMPFN) = TAVG
      AGEPFN(NUMPFN) = 1

! *** Calculate running average of leaf nitrogen 

	cntpflfcn(NUMPFN) = cntpflfcn(NUMPFN) + 1.0
	sumpflfcn(NUMPFN) = sumpflfcn(NUMPFN) + leafcn
!	avelfcn = sumpflfcn(NUMPFN-1) / cntpflfcn(NUMPFN-1)

! *** Calculate the age, weight and prefruiting leaf area at unfolding

      PfLfAge(numpfn) = PfLfAge(numpfn) 							 &
	                  + Duration_leaf_expansion(tavg)*calbrt(40) 
      if(PfLfAge(numpfn).lt.1.0)agepflf(numpfn) = agepflf(numpfn) + 1.0

! *** Assume that at node initiation, the leaf area is negligible

      PFWL(NUMPFN) = 0.001
      STEMWT = STEMWT - PFWL(NUMPFN)
      SLEAFN = SLEAFN + (PFWL(NUMPFN)* 0.035)
      STEMN  = STEMN  - (PFWL(NUMPFN)* 0.035)

      RETURN
      END


      SUBROUTINE ADDMSNOD(K)
! ************************************************************
! *** Increment mainstem node count and initiate leaf area ***
! *** NFBR(1) is set to 1 at first square                  ***
! ************************************************************

 	use common_block

      IF((K.GT.1).OR.(NFBR(K).GT.0)) THEN
        NFBR(K) = NFBR(K) + 1
	  nodeknt = numpfn + nfbr(k)
        NEWBR = NFBR(K)
        NNOD(K,NEWBR) = 1
        FFRUT(K,NEWBR,1) = 1.
        FCODE(K,NEWBR,1) = 1
        IF(PINHED.GE.1.) THEN
          PINHED = PINHED - 1.
          MCODE(K,NEWBR,1) = 6
        ENDIF

! *** Calculate running average of leaf nitrogen 

	  cntmslfcn(nodeknt) = cntmslfcn(nodeknt) + 1.0
	  summslfcn(nodeknt) = summslfcn(nodeknt) + leafcn
!	  avelfcn = summslfcn(NEWBR-1) / cntmslfcn(NEWBR-1)

! *** Calculate the age, weight and fruiting branch leaf area at unfolding 

        FRLfAge(k,NEWBR,1) = FRLfAge(k,NEWBR,1)						  &
			               + Duration_leaf_expansion(tavg)*calbrt(40)
        if(FRLfAge(k,NEWBR,1).lt.1.0) 								  &
	                    agefrlf(k,NEWBR,1) = agefrlf(k,NEWBR,1) + 1.0

! *** Assume that at node initiation, the leaf area is negligible

	  LEAFW(K,NEWBR,1) = 0.001
        STEMWT = STEMWT -  LEAFW(K,NEWBR,1)
        SLEAFN = SLEAFN + (LEAFW(K,NEWBR,1) * 0.035)
        STEMN  = STEMN  - (LEAFW(K,NEWBR,1) * 0.035)

! *** Calculate the age, weight and leaf area on main stem node at unfolding

        XMSLfAge(k,NEWBR) = XMSLfAge(k,NEWBR) 						  &
		              + Duration_leaf_expansion(tavg)*calbrt(40) 	  
        if(XMSLfAge(k,NEWBR).lt.1.0) 								  &
	                  agemslf(k,NEWBR) = agemslf(k,NEWBR) + 1.0

! *** Assume that at node initiation, the leaf area is negligible

	  MLEAFW(K,NEWBR) = 0.001
        STEMWT = STEMWT -  MLEAFW(K,NEWBR)
        SLEAFN = SLEAFN + (MLEAFW(K,NEWBR) * 0.035)
        STEMN  = STEMN  - (MLEAFW(K,NEWBR) * 0.035)
        NBRCH = NFBR(K)
        VDELAY(K) = 0.
        LEFCNT=LEFCNT+1
        LEFSRT(LEFCNT)=K*1000+NEWBR*10+1
      ELSE
        CALL FSTSQTODAY
      ENDIF

      RETURN
      END


	subroutine ADDFBNOD(k,l,nnid)
! *******************************************************************
! *** Increment fruiting branch node count and initiate leaf area ***
! *******************************************************************

 	use common_block

      NEWNOD = NNID + 1
      NNOD(K,L) = NNOD(K,L) + 1
      FFRUT(K,L,NEWNOD) = 1.
      FCODE(K,L,NEWNOD) = 1
	IF(PINHED.GE.1.) THEN
         MCODE(K,L,NEWNOD) = 6
         PINHED = PINHED - 1.
	endif
	nodeknt = numpfn + nfbr(k)

! *** Calculate running average of leaf nitrogen 

	  cntmslfcn(nodeknt) = cntmslfcn(nodeknt) + 1.0
	  summslfcn(nodeknt) = summslfcn(nodeknt) + leafcn
!	  avelfcn = summslfcn(nodeknt-1) / cntmslfcn(nodeknt-1)

! *** Calculate the age, weight and fruiting branch leaf area at unfolding 

      FRLfAge(k,L,NEWNOD) = FRLfAge(k,L,NEWNOD)							&
	    	               + Duration_leaf_expansion(tavg)*calbrt(40)
      if(FRLfAge(k,L,NEWNOD).lt.1.0) 									&
	                  agefrlf(k,L,NEWNOD) = agefrlf(k,L,NEWNOD) + 1.0

! *** Assume that at node initiation, the leaf area is negligible

      LEAFW(K,L,NEWNOD) = 0.001
      STEMWT = STEMWT -  LEAFW(K,L,NEWNOD)
      SLEAFN = SLEAFN + (LEAFW(K,L,NEWNOD) * 0.035)
      STEMN  = STEMN  - (LEAFW(K,L,NEWNOD) * 0.035)
      NNID = NNOD(K,L)
      DELAY(K,L) = 0.
      LEFCNT=LEFCNT+1
      LEFSRT(LEFCNT)=K*1000+L*10+NEWNOD
      RETURN
      END


      SUBROUTINE FSTSQTODAY

	use common_block

      ISQ=KDAY
      FCODE(1,1,1) = 1
      FFRUT(1,1,1) = 1.

! *** Calculate running average of leaf nitrogen 

	cntmslfcn(numpfn) = cntmslfcn(numpfn) + 1.0
	summslfcn(numpfn) = summslfcn(numpfn) + leafcn
!	avelfcn = summslfcn(numpfn-1) / cntmslfcn(numpfn-1)

! *** Calculate the age, weight and fruiting branch leaf area at unfolding 

      FRLfAge(1,1,1) = FRLfAge(1,1,1)									&
	    	               + Duration_leaf_expansion(tavg)*calbrt(40)
      if(FRLfAge(1,1,1).lt.1.0) 										&
	                  agefrlf(1,1,1) = agefrlf(1,1,1) + 1.0

! *** Assume that at node initiation, the leaf area is negligible

	LEAFW(1,1,1) = 0.001
      STEMWT = STEMWT -  LEAFW(1,1,1)
      SLEAFN = SLEAFN + (LEAFW(1,1,1) * 0.035)
      STEMN  = STEMN  - (LEAFW(1,1,1) * 0.035)
      NFBR(1) = 1
      NNOD(1,1) = 1
      FSTRES = 1.0
      LEFCNT = LEFCNT+1
      LEFSRT(LEFCNT)=1000+10+1
 
! *** DROP COTYLEDONS AT TIME OF FIRST SQUARE. 0.2 IS THE INITIAL LEAFW

 	daelfwt = 0.2
      LEFABS = LEFABS + daelfwt
	leafwt = leafwt - lefabs
      SLEAFN = SLEAFN - daelfwt*0.01
      NLOSS  = NLOSS + (daelfwt*0.01)
      PIXLOS = daelfwt * PIXCON
      PIXPLT = PIXPLT - PIXLOS
 
! *** Set the date counter for the calculation of ave. temp to FSQ

	isqday = daynum
	ifsqdae = kday
	ave_fsq_tavg = sum_fsq_tavg / float(isq)
	sum_fbl_tavg = sum_fsq_tavg
	ifbl = kday - 1

      RETURN
      END


!      SUBROUTINE QUALITY(K,L,M) 
!  *************************************************************  
!  *                                                           *  
!  *              FIBER QUALITY SUBROUTINE                     *  
!  *                                    from JLandivar's ICEMM  *  
!  *************************************************************  
!     
!      use common_block  

!      INCLUDE 'COTCOM.FOR'    
     
! FSX = FIBER STRENGTH : (G / TEX * 1/8 INCH)
! FLX = FIBER LENGTH ( INCHES, 2.5% SPUN )
!
!kit commented out 6/13/97

!      IF(IQUAL.EQ.0)THEN
!       WRITE(60,119)
!       WRITE(60,120)
! 119   FORMAT(/'              ----------- BOLL ---------  ----------'
!     .        ,' CROP -----------')     
! 120   FORMAT( 'DAE  K  L  M  LENGTH  MIKE  WEIGHT MATUR  LENGTH  '
!     .        ,'MIKE  WEIGHT  MATUR  YIELD',/)
!       IQUAL = 1
!      ENDIF

! ESTIMATE MICRONAIRE AS A FUNCTION OF POTENTIAL GROWTH POTENTIAL FOR
! FIBER AND FINAL FIBER WEIGHT PER BOLL

!      FMICRON(K,L,M) = 0.50 + .06*FMICRON(K,L,M)+0.80*FIBWGT(K,L,M)       

!      BOLLT = BOLTMP(K,L,M) + 9.0    
!      FS = 56.603 + BOLLT*(-2.921+0.059*BOLLT)    
!      FL = 1.219 -  0.0065*BOLTMP(K,L,M)  

!      GIN = GIN + (GINP * FFRUT(K,L,M))
!      FSX = FSX + FIBWGT(K,L,M)
!      FLX = FLX + (FLENGTH(K,L,M) * FFRUT(K,L,M))
!      FMK = FMK + (FMICRON(K,L,M) * FFRUT(K,L,M))
!      FMT = FMT + (FMATUR(K,L,M)  * FFRUT(K,L,M))

!      IF(NOPEN.GT.0.0) THEN 
!	   GIN     = GIN/NOPEN
!	   FWEIGHT = FSX/NOPEN  
!	   ALENGTH = FLX/NOPEN
!	   AMICRON = FMK/NOPEN
!	   AMATUR  = FMT/NOPEN
!      ENDIF
!     
!      IF(K.EQ.1.AND.M.EQ.1)THEN
!	 WRITE(60,100)DAE,K,L,M,FLENGTH(K,L,M)/25.4,FMICRON(K,L,M),
!     &             FIBWGT(K,L,M)/FFRUT(K,L,M),FMATUR(K,L,M),
!     &             ALENGTH/25.4,AMICRON,FWEIGHT,AMATUR,YIELD*500.0
!100      FORMAT(4I3,8F7.2,F8.2)
!      ENDIF
!
!      write(*,111)burwgt(k,l,m)/bolwgt(k,l,m), 
!     .            sedwgt(k,l,m)/bolwgt(k,l,m), 
!     .            fibwgt(k,l,m)/bolwgt(k,l,m),
!     .            bolwgt(k,l,m),k,l,m
!  111 format(4f6.3,3I3)
!
!      RETURN
!      END   


      SUBROUTINE ABCISE
! *************************************************************
! *                                                           *
! *                    ABSCISE SUBROUTINE                     *
! *                                                           *
! *************************************************************
 
	use common_block

      BOLLOS = 0.
      GBZ2   = 0.
      ABZ    = 0.
      SUMSQR = 0.
      SUMBOL = 0.
      SQRZ   = 0.
      SQRLOS = 0.
      SUMFRU = 0.
      FDROPS = 0.
      FDROPB = 0.
      SQABZ  = 0.
      BOLABZ = 0.
      DPSMX  = CALBRT(6)
      DPBMX  = CALBRT(7)
      DUM = LAI
      IF(DUM.GT.3.5) DUM = 3.5
      DROPLF = (90.0 - (5.0 * DUM)) * CALBRT(39)

! *** ABSCISE ALL LEAVES OLDER THAN DROPLF PHYSIOLOGICAL DAYS AND
! *** CALCULATE THE AMOUNT OF N REMOVED FROM THE PLANT SYSTEM BY THE
! *** ABSCISED LEAF. IT IS ASSUMED THAT DEAD LEAVES CONTAIN 1 % N.
! 
! *** ABSCISE PRE-FRUTING LEAVES
! 
      DO 4 J=1,NUMPFN
          IF((AGEPFN(J).GE.DROPLF).AND.(PFAL(J).GT.0.)			&
                                  .AND.(LAI.GT.0.3)) THEN
             AREA = AREA - PFAL(J)
             LAI = AREA/POPFAC
             PFAL(J) = 0.
             LEFABS = LEFABS + PFWL(J)
             LEAFWT = LEAFWT - PFWL(J)
             PIXLOS =  PFWL(J) * PIXCON
             PIXPLT = PIXPLT - PIXLOS
             SLEAFN = SLEAFN - PFWL(J) * 0.01
             NLOSS  = NLOSS + (PFWL(J) * 0.01)
             PFWL(J) = 0.
             IF((DEFBGN.GT.0).AND.(DAYNUM.GT.DEFBGN))			 &
                        LVSLOS = LVSLOS+1
          ENDIF
    4 CONTINUE

! *** ABSCISE MAIN STEM LEAVES
 
      DO 52 K=1,NVBRCH
           NBRCH = NFBR(K)
           DO 52 L=1,NBRCH
              NNID = NNOD(K,L)
              DO 52 M=1,NNID
                IF((M.EQ.1).AND.(LAGE(K,L,M).GE.DROPLF).AND.	 &
                       (MLAREA(K,L).GT.0.).AND.(LAI.GT.0.3)) THEN
                   LEFABS = LEFABS + MLEAFW(K,L)
                   LEAFWT = LEAFWT - MLEAFW(K,L)
                   PIXLOS = MLEAFW(K,L) * PIXCON
                   PIXPLT = PIXPLT - PIXLOS
                   SLEAFN = SLEAFN - MLEAFW(K,L) * 0.01
                   NLOSS  = NLOSS + (MLEAFW(K,L) * 0.01)
                   AREA = AREA - MLAREA(K,L)
                   LAI = AREA/POPFAC
                   MLAREA(K,L) = 0.
                   MLEAFW(K,L) = 0.
                   IF((DEFBGN.GT.0).AND.(DAYNUM.GT.DEFBGN))		 &
                       LVSLOS = LVSLOS+1
                ENDIF

 ! *** ABSCISE FRUITING LEAVES
 
                IF((LAGE(K,L,M).GE.DROPLF).AND.(LAREA(K,L,M).GT.0.)    &
                                          .AND.(LAI.GT.0.3)) THEN
                   LEFABS = LEFABS + LEAFW(K,L,M)
                   LEAFWT = LEAFWT - LEAFW(K,L,M)
                   PIXLOS = LEAFW(K,L,M) * PIXCON
                   PIXPLT = PIXPLT - PIXLOS
                   SLEAFN = SLEAFN - LEAFW(K,L,M) * 0.01
                   NLOSS  = NLOSS + (LEAFW(K,L,M) * 0.01)
                   AREA = AREA - LAREA(K,L,M)
                   LAI = AREA/POPFAC
                   LAREA(K,L,M) = 0.
                   LEAFW(K,L,M) = 0.								   
                   IF((DEFBGN.GT.0).AND.(DAYNUM.GT.DEFBGN))			   &
                       LVSLOS = LVSLOS+1
                ENDIF
   52 CONTINUE

      IF((DEFBGN.GT.0).AND.(DAYNUM.GT.DEFBGN)) THEN
        KNT = LFATDF*PERDEF/100.-LVSLOS
        I = 0
   60   CONTINUE
          I = I+1
          IF((KNT.GT.0).AND.(I.LE.LEFCNT)) THEN
            IDUM = LEFSRT(I)
            IF(IDUM.GT.0) THEN
              K = IDUM/1000
              L = IDUM/10-K*100
              M = IDUM-IDUM/10*10
			   IF(LEAFW(K,L,M).LE.0.0001) THEN
                    IDUM = 0
                    LEFSRT(I) = 0
                 ENDIF
            ENDIF
            IF(IDUM.GT.0) THEN
              LEFABS = LEFABS + LEAFW(K,L,M)
              LEAFWT = LEAFWT - LEAFW(K,L,M)
              PIXLOS = LEAFW(K,L,M) * PIXCON
              PIXPLT = PIXPLT - PIXLOS
              SLEAFN = SLEAFN - LEAFW(K,L,M) * 0.01
              NLOSS  = NLOSS + (LEAFW(K,L,M) * 0.01)
              AREA = AREA - LAREA(K,L,M)
              LAREA(K,L,M) = 0.
              LEAFW(K,L,M) = 0.
              KNT = KNT-1
              LVSLOS = LVSLOS+1
              IF((M.EQ.1).AND.(MLAREA(K,L).GT.0.)) THEN
                 LEFABS = LEFABS + MLEAFW(K,L)
                 LEAFWT = LEAFWT - MLEAFW(K,L)
                 PIXLOS = MLEAFW(K,L) * PIXCON
                 PIXPLT = PIXPLT - PIXLOS
                 SLEAFN = SLEAFN - MLEAFW(K,L) * 0.01
                 NLOSS  = NLOSS + (MLEAFW(K,L) * 0.01)
                 AREA = AREA - MLAREA(K,L)
                 MLAREA(K,L) = 0.
                 MLEAFW(K,L) = 0.
                 KNT=KNT-1
                 LVSLOS=LVSLOS+1
              ENDIF
              LEFSRT(I) = 0
            ENDIF
            GO TO 60
          ENDIF
        LAI = AREA/POPFAC
      ENDIF

! *** FRUIT ABSCISION

      IF(FCODE(1,1,1).gt.0) then 


        IF(FRATIO.LT.CALBRT(8)) THEN 
          FLOSS= (CALBRT(9)-3.607*FSTRES+1.605*FSTRES**2)*CALBRT(54)
        ELSE
          FLOSS=(CALBRT(10)-3.607*FSTRES+1.605*FSTRES**2)*CALBRT(55)
        ENDIF
        IF(FLOSS.LT.0.0) FLOSS = 0.0

! *** FCODES   1=SQ  2=GB  3=OB  4  5  6  or 7=YGB
! *** Count the total number of squares/bolls susceptible to abscission

        DO K=1,NVBRCH
           NBRCH = NFBR(K)
           DO L=1,NBRCH
              NNID = NNOD(K,L)
              DO M=1,NNID
                 IF(FCODE(K,L,M).EQ.1) THEN
                    IF(AGE(K,L,M).GT.2.0) THEN
                       SUMSQR=SUMSQR+FFRUT(K,L,M)
                    ENDIF
                 ELSEIF(FCODE(K,L,M).EQ.7) THEN
                    IF(AGEBOL(K,L,M).LT.AGEABZ(K,L,M)) THEN
                       SUMBOL = SUMBOL + FFRUT(K,L,M)
                    ENDIF
                 ENDIF
              ENDDO
           ENDDO
        ENDDO 
	susceptible_bolls = sumbol
        SUMFRU = SUMBOL + SUMSQR
 
! *** DO CALCULATION PERTAINING TO YOUNG GREEN BOLLS. REMOVE BOLLS WITH
! *** A SMALL FRACTION REMAINING, SET FCODE TO 4 AND REDUCE ITS WEIGHT AND
! *** NITROGEN CONTENT

        IF(SUMBOL.GT.0.001) THEN
	       call Running_Ave_Temp
	       call TempFact_Boll

! *** Increment the number of bolls marked for abcission due to	
! *** nutritional and heat stress         Kit's original code  04/27/00
! *** Count the number of bolls marked for abcission due to heat stress
 
	       HeatLoss = HeatIndex * (SUMBOL/SUMFRU)

! *** Count the number of bolls marked for abcission due to nutritional stress
 
           BOLLOS = FLOSS * (SUMBOL/SUMFRU)

           bollos = HeatLoss + Bollos 
	       if(bollos.gt.sumbol) bollos = sumbol

           BOLLOSMX = SUMBOL * DPBMX
           IF(POLOSS.GT.BOLLOSMX) THEN
              FDROPB = 0.0
              BOLLOS = POLOSS
           ELSEIF((BOLLOS + POLOSS).GT.BOLLOSMX) THEN
              FDROPB = (BOLLOSMX - POLOSS) / SUMBOL
              BOLLOS = DPBMX * SUMBOL
           ELSE
              FDROPB = (BOLLOS - POLOSS) / SUMBOL
              IF(FDROPB.LT.0.) FDROPB = 0.
              BOLLOS = FDROPB * SUMBOL + POLOSS
           ENDIF
           SUMBOL = SUMBOL - SUMBOL*FDROPB
        ENDIF

        IF(SUMSQR.GT.0.00001) THEN
           SQRLOS = FLOSS - BOLLOS
           IF(SQRLOS.LT.0.0) SQRLOS = 0.0
           FDROPS = SQRLOS / SUMSQR
           IF(FDROPS.GT.DPSMX) THEN
              FDROPS = DPSMX
              SQRLOS = SUMSQR * DPSMX
           ENDIF
           SUMSQR = SUMSQR - SQRLOS
        ENDIF

        DO K=1,NVBRCH
           NBRCH = NFBR(K)
           DO L=1,NBRCH
              NNID = NNOD(K,L)
              DO M=1,NNID
                 WTLOS = 0.0
                 IF(FCODE(K,L,M).EQ.1) THEN
                    IF(AGE(K,L,M).GT.2.0) THEN
                       WTLOS = SQRWT(K,L,M) * FDROPS
                       SQRWT(K,L,M) = SQRWT(K,L,M) - WTLOS
                       PQFLR = PQFLR + WTLOS
                       SQWT = SQWT - WTLOS
                       SQABZ = SQABZ + FFRUT(K,L,M)*FDROPS
                       FFRUT(K,L,M) = FFRUT(K,L,M)-FFRUT(K,L,M)*FDROPS
                       IF(FFRUT(K,L,M).LE.0.001) THEN
                          FFRUT(K,L,M) = 0.0
                          PQFLR = PQFLR + SQRWT(K,L,M)
                          SQWT = SQWT - SQRWT(K,L,M)
                          SQRWT(K,L,M) = 0.0
                       ENDIF
                       SQRZ = SQRZ + FFRUT(K,L,M)
                    ELSEIF(AGE(K,L,M).EQ.2) THEN
                       PINDROP = .045
                       WTLOS = SQRWT(K,L,M) * PINDROP
                       SQRWT(K,L,M) = SQRWT(K,L,M) - WTLOS
                       PQFLR = PQFLR + WTLOS
                       SQWT = SQWT - WTLOS
                       SQABZ = SQABZ + FFRUT(K,L,M) * PINDROP
                       FFRUT(K,L,M)=FFRUT(K,L,M)-FFRUT(K,L,M)*PINDROP
                       SQRZ = SQRZ + FFRUT(K,L,M)
                    ELSE
                       SQRZ = SQRZ + FFRUT(K,L,M)
                    ENDIF
                 ENDIF
                 IF(FCODE(K,L,M).EQ.2) THEN
                    GBZ2 = GBZ2 + FFRUT(K,L,M)
                 ENDIF
                 IF(FCODE(K,L,M).EQ.7) THEN
                    IF(MCODE(K,L,M).EQ.8) THEN
                       FFRUT(K,L,M) = FFRUT(K,L,M)*.75
                       MCODE(K,L,M) = 7
                    ELSE
                       CALL ABCISBOL(BOLABZ,K,L,M,FDROPB)
                    ENDIF
                    IF(FFRUT(K,L,M).LE.0.001) THEN
                       SEEDN = SEEDN - BOLWGT(K,L,M)*0.416*SEEDCN
                       BURRN = BURRN - BOLWGT(K,L,M)*0.278*BURCN
                       NLOSS = NLOSS + BOLWGT(K,L,M)*0.416*SEEDCN
                       NLOSS = NLOSS + BOLWGT(K,L,M)*0.278*BURCN
                       GBLOS = GBLOS + BOLWGT(K,L,M)
                       GBOLWT = GBOLWT - BOLWGT(K,L,M)
                       PIXLOS = BOLWGT(K,L,M) * PIXCON
                       PIXPLT = PIXPLT - PIXLOS
                       BOLWGT(K,L,M) = 0.0
                       BOLABZ = BOLABZ + FFRUT(K,L,M)
                       FFRUT(K,L,M) = 0.0
                    ENDIF
                    GBZ2 = GBZ2 + FFRUT(K,L,M)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        SQLOSS(DAYNUM) = SQABZ
        BOLOSS(DAYNUM) = BOLABZ
      ABZ = SQABZ + BOLABZ

      endif

      RETURN
      END


	Subroutine Running_Ave_Temp
! ***********************************************************
! *** CALCULATE RUNNING AVERAGE OF TAVG FOR THE LAST FIVE ***  
! *** DAYS. THIS IS DONE TO SMOOTH THE SHARP FLUCTUATIONS *** 
! *** IN DAILY TAVG.MOVE ALL VALUES DOWN ONE SLOT, PUT  	***
! *** NEW VALUE IN FIRST SLOT, AVERAGE 					***
! ***********************************************************

	use common_block

      do j = 2,7
         Bloom_tavg(j) = Bloom_tavg(daynum-j+1)
      enddo
      Bloom_tavg(1) = Bloom_tavg(daynum)

!  Used only when considering the current effect of canopy 
!  temperature differential                      06/04/03  Kit
!      Bloom_tavg(1) = tavg + Canopy_Temp_Differential(psild)

      Bloom_tavg(1) = tavg
      Boll_tavg = 0.
      do J =1,7
         Boll_tavg = Boll_tavg + Bloom_tavg(J)
      enddo
      Boll_tavg = Boll_tavg/7.0
	return
	end

	   
	Subroutine TempFact_Boll
! ********************************************************
! *** Calculates heat stress index for boll abscission ***
! *** based on percent fruit retention vs temperature. ***
! *** TempFact_Node is a multiplier to increase boll   *** 
! *** loss. Based temperature of 26.0 oC.  HeatIndex   ***
! *** is the number of bolls to be dropped today.      ***
! *** (Source:  K.R. Reddy)             01/18/2000	 ***
! ********************************************************


	use common_block

	stdtemp = 26.0
	templimit = 33.5
	if(Boll_tavg.lt.stdtemp) Boll_tavg = stdtemp
	if(Boll_tavg.gt.templimit) Boll_tavg = templimit
	dum1 = Boll_tavg/30.7978
	dum2 = stdtemp/30.7978
	tempnum = -7.2531 + (92.2099/(1.0 + dum1**28.5048))
	tempden = -7.2531 + (92.2099/(1.0 + dum2**28.5048))
	ratio = tempnum/tempden
	if(ratio.gt.1.0) ratio = 1.0
	if(ratio.lt.0.0) ratio = 0.0
	HeatIndex = (1.0 - ratio)

!  Number of bolls lost due to heat injury
	HeatIndex0 = (1.0 - ratio) * susceptible_bolls * calbrt(50)

	return
	end


      SUBROUTINE ABCISBOL(BOLABZ,K,L,M,FDROPB)

	use common_block

      FDUM = FDROPB
      IF(FFRUT(K,L,M)-FFRUT(K,L,M)*FDUM.GT.1.) THEN
        FDUM = (FFRUT(K,L,M)-1.0)/FFRUT(K,L,M)
      ENDIF
      WTLOS = BOLWGT(K,L,M) * FDUM
      SEEDN = SEEDN - WTLOS*0.416*SEEDCN
      BURRN = BURRN - WTLOS*0.278*BURCN
      NLOSS = NLOSS + WTLOS*0.416*SEEDCN
      NLOSS = NLOSS + WTLOS*0.278*BURCN
      GBLOS = GBLOS + WTLOS
      GBOLWT = GBOLWT - WTLOS
      PIXLOS = WTLOS * PIXCON
      PIXPLT = PIXPLT - PIXLOS
      BOLWGT(K,L,M) = BOLWGT(K,L,M) - WTLOS
      BOLABZ = BOLABZ + FFRUT(K,L,M)*FDUM
      FFRUT(K,L,M) = FFRUT(K,L,M)-FFRUT(K,L,M)*FDUM

      RETURN
      END


      SUBROUTINE PMAPS
! *************************************************************
! *                                                           *
! *                PLANT MAPS SUBROUTINE                      *
! *                                                           *
! *************************************************************
! * SETS MCODES TO DISTRIBUTE ABSCISSION OVER ENTIRE PLANT.   *
! * A FRUIT WILL BE INCLUDED IN THE AVERAGE PLANT MAP IF THE  *
! * FRACTION OF FRUIT REMAINING AT EACH SITE IS EQUAL OR      *
! * LARGER TO THE AVERAGE SIZE FOR EACH FRUIT CLASS.          *
! * THIS SUBROUTINE ALSO CALCULATES SUMMARIES OF FRUIT PRESENT*
! * AT EACH SITE AND CALLS COTPLT TO PRINT PLANT MAPS         *
! *************************************************************
! *                                                           *
! * FCODE                          MCODE                      *
! *                                                           *
! * 1 SQUARE                       1 SQUARE                   *
! * 2 GREEN BOLL                   2 GREEN BOLL               *
! * 3 MATURE BOLL                  3 MATURE BOLL              *
! * 4 FRUIT ABSCISED AS A SQUARE   4 FRUIT ABZ AS A SQUARE    *
! * 5 FRUIT ABSCISED AS A BOLL     5 FRUIT ABZ AS A BOLL      *
! * 6 NOT USED                     6 ABZ AS PINHEAD OR POLYNA *
! * 7 YOUNG GREEN BOLL             7 YOUNG GREEN BOLL         *
! *************************************************************
 
      use common_block
 
         SUMFRU1 = 0.0
         SUMFRU2 = 0.0
         SUMFRU3 = 0.0
         SUMFRU7 = 0.0
         FRU1NO  = 0.0
         FRU2NO  = 0.0
         FRU3NO  = 0.0
         FRU7NO  = 0.0

!         DO 23 K = 1,450

         DO 23 K = 1,mxfruts
         IF(FCODES(K).NE.0) THEN
           IF(FCODES(K).EQ.1) THEN
             SUMFRU1 = SUMFRU1 + FRUITS(K)
             FRU1NO  = FRU1NO  + 1.0
           ENDIF
           IF(FCODES(K).EQ.2) THEN
             SUMFRU2 = SUMFRU2 + FRUITS(K)
             FRU2NO  = FRU2NO  + 1.0
           ENDIF
           IF(FCODES(K).EQ.3) THEN
             SUMFRU3 = SUMFRU3 + FRUITS(K)
             FRU3NO  = FRU3NO  + 1.0
           ENDIF
           IF(FCODES(K).EQ.7) THEN
             SUMFRU7 = SUMFRU7 + FRUITS(K)
             FRU7NO  = FRU7NO  + 1.0
           ENDIF
        ENDIF
 23     CONTINUE
        IF(FRU1NO.GT.0.0) AVGFRU1 = SUMFRU1/FRU1NO
        IF(FRU2NO.GT.0.0) AVGFRU2 = SUMFRU2/FRU2NO
        IF(FRU3NO.GT.0.0) AVGFRU3 = SUMFRU3/FRU3NO
        IF(FRU7NO.GT.0.0) AVGFRU7 = SUMFRU7/FRU7NO

!         DO 24 K = 1,450

         DO 24 K = 1,mxfruts
         IF(FCODES(K).NE.0) THEN
           IF(FCODES(K).EQ.1) THEN
             IF(FRUITS(K).GE.AVGFRU1) THEN
               MCODES(K) = 1
             ELSE
               MCODES(K) = 4
             ENDIF
           ENDIF
           IF(FCODES(K).EQ.2) THEN
             IF(FRUITS(K).GE.AVGFRU2) THEN
               MCODES(K) = 2
             ELSE
               MCODES(K) = 5
             ENDIF
           ENDIF
           IF(FCODES(K).EQ.3) THEN
             IF(FRUITS(K).GE.AVGFRU3) THEN
               MCODES(K) = 3
             ELSE
               MCODES(K) = 5
             ENDIF
          ENDIF
           IF(FCODES(K).EQ.7) THEN
             IF(FRUITS(K).GE.AVGFRU7) THEN
               MCODES(K) = 7
             ELSE
               MCODES(K) = 5
             ENDIF
           ENDIF
        IF(FCODES(K).EQ.4) MCODES(K) = 4
        IF(FCODES(K).EQ.5) MCODES(K) = 5
        ENDIF
 24     CONTINUE

! *** DISPLAYS OUTPUT OF BOLL AGE AND WEIGHT. **NOTE** FRACTION 
! *** OF FRUTP AND BSIZE ARE TRUNKATED IN COTPLT

           DO 100 K=1,NVBRCH
              NBRCH = NFBR(K)
              DO 100 L=1,NBRCH
                 NNID = NNOD(K,L)
                 DO 100 M=1,NNID
                 IF(FCODE(K,L,M).EQ.2.OR.FCODE(K,L,M).EQ.3.OR.	  &
                 FCODE(K,L,M).EQ.7) THEN
                   IF(FFRUT(K,L,M).GT.0.001) THEN
                     BSIZE(K,L,M) = BOLWGT(K,L,M)/FFRUT(K,L,M)
                   ELSE
                     BSIZE(K,L,M) = 0.0
                   ENDIF
                 IF(BSIZE(K,L,M).GT.CALBRT(12))BSIZE(K,L,M)=CALBRT(12)  
                 ENDIF
                 IF(FCODE(K,L,M).EQ.4.OR.FCODE(K,L,M).EQ.5) THEN
                   FRUTP(K,L,M) = 11.
                   BSIZE(K,L,M) = 11.
                   MATURE(K,L,M)= 11
                 ELSE
                   FRUTP(K,L,M) = FFRUT(K,L,M)*10.0
                   MATURE(K,L,M)= (DEHISS(K,L,M)-AGEBOL(K,L,M))/10.
                   IF(FRUTP(K,L,M).LT.1.)FRUTP(K,L,M)  = 10.
                   IF(BSIZE(K,L,M).LT.1.)BSIZE(K,L,M)  = 10.
                   IF(MATURE(K,L,M).LT.1)MATURE(K,L,M) = 10
                   IF(FCODE(K,L,M).EQ.1) BSIZE(K,L,M)  = 13.
                   IF(FCODE(K,L,M).EQ.1) MATURE(K,L,M) = 13
                   IF(FCODE(K,L,M).EQ.3) MATURE(K,L,M) = 12
                 ENDIF
 100       CONTINUE
           
        RETURN
        END


      SUBROUTINE OUTPUT
!                          O U T P U T
!  ****************************************************************
!  *   COLLECTION OF WRITE STATEMENTS FOR OUTPUT OF MODEL RESULTS *
!  ****************************************************************

      use common_block

 1000 FORMAT(' ** FIRST SQUARE ON ',I2.2,'/',I2.2,' **')
 1020 FORMAT(11X,' YIELD BALES/ACRE   YIELD LBS/ACRE')
 1030 FORMAT(17X,F4.2,13X,F7.0)
 1040 FORMAT(' ** FIRST BLOOM ON ',I2.2,'/',I2.2,' **')
 1050 FORMAT(1X,I3,1X,I2,'/',I2.2,F7.1,F6.2,I7,2X,I8,1X,3I8,F7.2,I8)
 1060 FORMAT('1',19X,'GENERAL OUTPUT IN ENGLISH UNITS')
 1065 FORMAT(17X,'PER/PLANT',27X,'PER/ACRE')
 1070 FORMAT(12X,18('-'),5X,45('-'))
 1080 FORMAT(' DAE  DATE  HEIGHT  LAI  NODES     SITES',		 &
             '  SQUARES   GREEN   OPEN   YIELD   FRUIT')
 1090 FORMAT(13X,'(IN)',35X,'BOLLS   BOLLS',10X,' SHED')
 1100 FORMAT(' ')
 1110 FORMAT(' **',F6.1,' LBS OF NITROGEN APPLIED ON ',			 &
                                          I2.2,'/',I2.2,' **')	 
 1120 FORMAT(' **',F6.2,' INCHES OF IRRIGATION APPLIED ON ',	 &
                                          I2.2,'/',I2.2,' **')	 
 1140 FORMAT(' FIRST SQUARE  ',2X,I2.2,'/',I2.2,I5,F7.1,I7,		 &
             2F6.1,3I7,F6.2)
 1180 FORMAT(' ** LAST ACTUAL WEATHER ',I2.2,'/',I2.2,'**')
 1200 FORMAT(' LAST ACT WTHER',2X,I2.2,'/',I2.2,I5,F7.1,I7,		 &
             2F6.1,3I7,F6.2)
 1220 FORMAT(' FIRST BLOOM   ',2X,I2.2,'/',I2.2,I5,F7.1,I7,		 &
             2F6.1,3I7,F6.2)
 1300 FORMAT(' ** PREP APPLIED ON ',I2.2,'/',I2.2,'**')
 1320 FORMAT(' ** DEFOLIANT APPLIED ON ',I2.2,'/',I2.2,'**')
 1340 FORMAT(' ** ',F6.1,' OZS PIX APPLIED ON ',I2.2,'/',I2.2,' **')
 1360 FORMAT(' ** PLANT MAP ADJUSTMENT ',I2.2,'/',I2.2,			 &
             ' USED',I5,' PLANTS **')
 1380 FORMAT(' ** PLANT HT. ADJUSTMENT ',I2.2,'/',I2.2,			 &
             ' USED',I5,' PLANTS **')
 1400 FORMAT(' ** SIMULATED CROP MAY BE IMMATURE. ',			 &
             'YOU MAY WANT TO TRY A LATER STOP DATE. **')
 1420 FORMAT(' ** TEMPERATURE FELL BELOW 30 F.  IT IS OVER. **')
 1440 FORMAT(' ** TEMPERATURE FELL BELOW 32 F. **')
 1460 FORMAT(' ** Oh Boy - Leaf Weight Went To Zero. **')

      C00 = Z / 2.54
      C01 = INT * 100.
!      C02 = PIXCON * 1000000
      C02 = PIXCON * 1000
      C03 = PSILD / .101325
      N00 = NUMPFN + NFBR(1)
      I00=SQRZ*POPPLT/1000+.5
      I01=GBZ2*POPPLT/1000+.5
      I02=NOPEN*POPPLT/1000+.5
      AVAILN = AVAILN/.0448
      DAILYN=SUPNO3+SUPNH4+(LEFABS-PLEFABS)*.005
      PLEFABS=LEFABS
      FRTNUSE=.045*.416*SDWBOL+.006*.278*SDWBOL
      VEGNUSE=.009*SDWSTM+.015*SDWLEF+.019*RCH2O
      DYNEXS=(POPPLT/454)*(DAILYN-FRTNUSE-VEGNUSE)

      WRITE(21,ERR=300) IDAY,KDAY,YIELD,NF,WSTRSD,C00,RI,TMAX,TMIN,SQRZ	   &
            ,RAIN,FERN,GBZ2,C03,NV,LAI,SITEZ,NOPEN,FLOSS,ABZ			   &
            ,CSTRES,PSIAVG,TAVG,TDAY,TNYT,C01,WIND,ES,EP				   &
            ,RN,PN,CUMES,CUMEP,CUMSOK,CUMRAN,TH2O,H2OBAL,SPN,PLANTW	       &
            ,CHOBAL,LEAFWT,STEMWT,ROOTWT,SQWT,COTXX,GBOLWT,DEADWT		   &
            ,XTRAC,POPPLT,LEAFCN,ISQ,DAZE,MO,N00,C02,NUMPFN			       &
            ,SDWSTM,SDWLEF,RCH2O,SDWBOL,AVAILN,SUPNO3,SUPNH4			   &
            ,LEFABS,RUNOFF(DAYNUM),AMTIRR(DAYNUM)
  300 CONTINUE

      IF(DAYNUM.GE.EMERGE) THEN
        N00=NUMPFN+NFBR(1)
        N01=IFIX(SITEZ*POPPLT)
        N02=IFIX(SQRZ*POPPLT)
        N03=IFIX(GBZ2*POPPLT)
        N04=IFIX(NOPEN*POPPLT)
        N05=IFIX(ABZ*POPPLT)
        A01=Z/2.54

! *** data for calibration and validation

        m01=IFIX(SITEZ)
        m02=IFIX(SQRZ)
        m03=IFIX(GBZ2)
        m04=IFIX(NOPEN)
        m05=IFIX(ABZ)
        m06=IFIX(abzb*POPPLT)
        YLBSAC=YIELD*500.
        TAVG=TAVG/5.*9.+32.
        RAIN=RAIN/25.4

        WRITE(77,1925) DAYNUM,MO,DAZE,A01,LAI,N00,m01,m02,m03,m04,	 &
                       SQWT,gbolwt,cotxx,leafwt,stemwt,ylbsac,		 &
	                   tavg,rain,ri,wstrsd,supno3,supnh4
 1925   FORMAT(1X,I3,1X,I2,'/',I2.2,2F6.2,5I4,5F7.2,F9.2,3F8.2,3f7.3)
 
! TEST BEGIN
      IF((IDAY.EQ.1).OR.(IDAY/NFRQ*NFRQ.EQ.IDAY)) THEN
          IF(FCODE(1,1,1).NE.0.AND.NPC.GT.0) CALL COTPLT(1)
          IF(NPN.GT.0) CALL OUT(VNO3C,TNNO3,1)
          IF(NPN.GT.0) CALL OUT(VNH4C,TNNH4,1)
          IF(NPW.GT.0) CALL OUT(VH2OC,TH2O,2)
          IF(NPR.GT.0) CALL OUT(ROOTSV,ROOTWT,3)
          IF(NPP.GT.0) CALL OUT(PSIS,PSIAVG,4)
          IF(NPT.GT.0) CALL OUT(SOILT,TSOLAV(2),5)
      ENDIF
! TEST END
 

        LINE=LINE+1
      ENDIF
      
      !WATER OUTPUT -- BY GR August 12, 2008 BEGIN-
      do 813 L = 1, NL
        NINDEX = 1
        VH2OC_A(L) = 0.0
        do 814 N  = 1, NK
            if(isNaN(VH2OC(L,N)) .NE. .TRUE.) then
               VH2OC_A(L) = VH2OC_A(L)+VH2OC(L,N)
               NINDEX = NINDEX +1
            endif
814     continue            
        VH2OC_A(L) = VH2OC_A(L) / NINDEX
        WRITE(88,2008) DAYNUM,MO,DAZE,L, VH2OC_A(L)
813   continue    
       

2008  FORMAT(1X,I3,1X,I2,'/',I2.2, I4,F6.2)
  !WATER OUTPUT --  BY GR August 12, 2008  END
    
      IF(DAYNUM.GE.JDSTPS) SEND=.TRUE.
      IF(LEAFWT.LE.0.) SEND=.TRUE.
      IF((CLIMAT(DAYNUM,3).LT.30.).AND.(DAYNUM.GT.EMERGE)) SEND=.TRUE.

      RETURN
      END


      SUBROUTINE COTPLT(MTYPE)
!  ************************************************************
!  *                                                          *
!  *  THIS SUBROUTINE PLOTS THE COTTON PLANT WITH THE         *
!  *  CORRESPONDING SYMBOLS IN THE CORRECT LOCATIONS.         *
!  *  SQUARE-X             OLDER GREEN BOLL-*    OPEN BOLL-$  *
!  *  YOUNG GREEN BOLL-B   ABSCISED-O                         *
!  *                                       Original Code      *
!  ************************************************************
! 
      use common_block
! 
 95   FORMAT(12X,'AVERAGE PLANT FRUIT CODE ',I2.2,'/',I2.2				 &
          /,5X,'*  X  SQUARE          B  BLOOM OR YOUNG GREEN BOLL*'/	 &
            5X,'*  *  SET GREEN BOLL  $  MATURE BOLL              *'/	 &
            5X,'*  A  ABSCISED FRUIT  I  MONOPODIAL NODE          *'/	 &
            5X,'***************************************************'/)
 96   FORMAT(19X,'FRUIT CODE MAP ',I2.2,'/',I2.2						 &
          /,5X,'*  X  SQUARE          B  BLOOM OR YOUNG GREEN BOLL*'/	 &
            5X,'*  *  SET GREEN BOLL  $  MATURE BOLL              *'/	 &
            5X,'*  A  ABSCISED FRUIT  I  MONOPODIAL NODE          *'/	 &
            5X,'***************************************************'/)
 97   FORMAT(10X,'PERCENTAGE OF FRUIT BY POSITION ',I2.2,'/',I2.2		 &
          /,5X,'*      A -->  ABSCISED                            *'/	 &
            5X,'*      0 -->  0 -  9 %      1 --> 10 -  19 %      *'/	 &
            5X,'*      2 --> 20 - 29 %      3 --> 30 -  39 %      *'/	 &
            5X,'*      4 --> 40 - 49 %      5 --> 50 -  59 %      *'/	 &
            5X,'*      6 --> 60 - 69 %      7 --> 70 -  79 %      *'/	 &
            5X,'*      8 --> 80 - 89 %      9 --> 90 - 100 %      *'/	 &
            5X,'***************************************************'/	 &
            '  ' / '  ')
 98   FORMAT(14X,'BOLL WEIGHT BY POSITION ',I2.2,'/',I2.2				 &
          /,5X,'*      A -->  ABSCISED      X -->  SQUARE         *'/	 &
            5X,'*      0 -->  0.0 - 0.9 g   1 -->  1.0 - 1.9 g    *'/	 &
            5X,'*      2 -->  2.0 - 2.9 g   3 -->  3.0 - 3.9 g    *'/	 &
            5X,'*      4 -->  4.0 - 4.9 g   5 -->  5.0 - 5.9 g    *'/	 &
            5X,'*      6 -->  6.0 - 6.9 g   7 -->  7.0 - 7.9 g    *'/	 &
            5X,'***************************************************'/)
 99   FORMAT(11X,'DAYS TO MATURITY BY POSITION ',I2.2,'/',I2.2			 &
          /,5X,'*         A -->  ABSCISED    X --> SQUARE         *'/	 &
            5X,'*         $ -->  MATURE      0 -->   0 -  9       *'/	 &
            5X,'*         1 -->  10 - 19     2 -->  20 - 29       *'/	 &
            5X,'*         3 -->  30 - 39     4 -->  40 - 49       *'/	 &
            5X,'*         5 -->  50 - 59     6 -->  60 - 69       *'/	 &
            5X,'***************************************************'/)
! 
  100 FORMAT(14X,A1,6A2,2(13X,A1,6A2))
  101 FORMAT(2X,6A2,A1,2(13X,6A2,A1))									 
  102 FORMAT(14X,A1,25X,A1,25X,A1/14X,A1,25X,A1,25X,A1/					 &
             14X,A1,25X,A1,25X,A1)
  103 FORMAT(//)
! 
        WRITE(24,*)'FROM COTPLT SUB'
      WRITE(24,103)
      
      DO 1 K=1,NVBRCH
        NBRCH=NFBR(K)
        DO 2 L=1,NBRCH,2
          NNID  =NNOD(K,L)
          DO 3 M=1,NNID
            IF(FCODE(K,L,M).EQ.0)GO TO 2
            IF(MTYPE.EQ.1) THEN
              ICODE = MCODE(K,L,M)
              IF(ICODE.GT.0) PRT(K,L,M) = CHAR1(ICODE)
            ENDIF
            IF(MTYPE.EQ.2) THEN
              ICODE = FCODE(K,L,M)
              IF(ICODE.GT.0) PRT(K,L,M) = CHAR1(ICODE)
            ENDIF
            IF(MTYPE.EQ.3) THEN
              ICODE = FRUTP(K,L,M)
              IF((FRUTP(K,L,M).GE.10.).AND.(FRUTP(K,L,M).LT.10.01))		 &
                                                              ICODE=9
              IF(ICODE.GT.0) PRT(K,L,M) = CHAR3(ICODE)
            ENDIF
            IF(MTYPE.EQ.4) THEN
              ICODE = BSIZE(K,L,M)
              IF(ICODE.GT.0) PRT(K,L,M) = CHAR3(ICODE)
            ENDIF
            IF(MTYPE.EQ.5) THEN
              ICODE = MATURE(K,L,M)
              IF(ICODE.GT.0) PRT(K,L,M) = CHAR3(ICODE)
            ENDIF
            PRI(K,L) = CHARI
 3        CONTINUE
 2      CONTINUE
        DO 4 L=2,NBRCH,2
          NNID = NNOD(K,L)
          DO 5 M=1,NNID
            M2 = 7-M
!            M2 = 6-M
            IF(FCODE(K,L,M).EQ.0) GO TO 4
            IF(MTYPE.EQ.1) THEN
              ICODE = MCODE(K,L,M)
              IF(ICODE.GT.0) PRT(K,L,M2) = CHAR2(ICODE)
            ENDIF
            IF(MTYPE.EQ.2) THEN
              ICODE = FCODE(K,L,M)
              IF(ICODE.GT.0) PRT(K,L,M2) = CHAR2(ICODE)
            ENDIF
            IF(MTYPE.EQ.3) THEN
              ICODE = FRUTP(K,L,M)
              IF((FRUTP(K,L,M).GE.10.).AND.(FRUTP(K,L,M).LT.10.01))	   &
                                                              ICODE=9
              IF(ICODE.GT.0) PRT(K,L,M2) = CHAR4(ICODE)
            ENDIF
            IF(MTYPE.EQ.4) THEN
              ICODE = BSIZE(K,L,M)
              IF(ICODE.GT.0) PRT(K,L,M2) = CHAR4(ICODE)
            ENDIF
            IF(MTYPE.EQ.5) THEN
              ICODE = MATURE(K,L,M)
              IF(ICODE.GT.0) PRT(K,L,M2) = CHAR4(ICODE)
            ENDIF
            PRI(K,L) = CHARI
 5        CONTINUE
 4      CONTINUE
 1    CONTINUE
      IF(MTYPE.EQ.1) WRITE(24,95)MO,DAZE
      IF(MTYPE.EQ.2) WRITE(24,96)MO,DAZE
      IF(MTYPE.EQ.3) WRITE(24,97)MO,DAZE
      IF(MTYPE.EQ.4) WRITE(24,98)MO,DAZE
      IF(MTYPE.EQ.5) WRITE(24,99)MO,DAZE
      NBRCH = NFBR(1)
      IF(NBRCH.EQ.(NBRCH/2*2)) GO TO 20
          DO 10 L=1,NBRCH,2
              LX = NBRCH +1 - L
              LX1 = LX-1
              WRITE(24,100) PRI(1,LX),(PRT(1,LX,M),M=1,6),PRI(2,LX),	 &
              (PRT(2,LX,M),M=1,6),PRI(3,LX),(PRT(3,LX,M),M=1,6)
              IF(LX1.EQ.0) GO TO 6
                  WRITE(24,101) (PRT(1,LX1,M),M=1,6),PRI(1,LX1),		 &
                  (PRT(2,LX1,M),M=1,6),PRI(2,LX1),(PRT(3,LX1,M),M=1,6),  &
                   PRI(3,LX1)
 10       CONTINUE
      GO TO 6
 20       DO 21 L=1,NBRCH,2
              LX = NBRCH + 1 - L
              LX1 = LX - 1
              WRITE(24,101) (PRT(1,LX,M),M=1,6),PRI(1,LX),				 &
              (PRT(2,LX,M),M=1,6),PRI(2,LX),(PRT(3,LX,M),M=1,6),		 &
             PRI(3,LX)
              IF(LX1.EQ.0) GO TO 6
                  WRITE(24,100) PRI(1,LX1),(PRT(1,LX1,M),M=1,6),		 &
                  PRI(2,LX1),(PRT(2,LX1,M),M=1,6),PRI(3,LX1),			 &
                 (PRT(3,LX1,M),M=1,6)
 21       CONTINUE
 6    CONTINUE
      WRITE(24,102) PRI(1,1),PRI(2,1),PRI(3,1),PRI(1,1),PRI(2,1),		 &
                    PRI(3,1),PRI(1,1),PRI(2,1),PRI(3,1)
      WRITE(24,103)
      RETURN
      END


      SUBROUTINE OUT(ARRAY,TOTAL,IGO)
!  ************************************************************
!  *                                                          *
!  *  THIS SUBROUTINE PLOTS THE SOIL SLAB AND THE DENSITIES   *
!  *  OF THE ARRAY ELEMENTS IN EACH CELL.                     *
!  *                                                          *
!  ************************************************************

      use common_block

      DIMENSION CAPSCA(11), PSISCA(11), VNOSCA(11),						&
                ROOSCA(11), RANGE(11), ARRAY(40,20),TEMPSCA(11)

      CHARACTER TTL1*40, TTL3*40, TTL4*40, TTL5*40,TTL8*40,				&
                TTL7*40, UNITST*16, VNOUNI*24, VH2UNI*24, PSIUNI*24,	&
                NITUNT*16, TTL1R*40, TTL2R*40, UNITS*24, UNITSR*16,		&
                UNTS*24, UNTST*16, TL1*40, TL2*40, TEMPUNT*16,TTL9*40,	&
                TUNITS*24, UNITPS*16

      DATA ROOSCA/1.0E-25,.0001,.0005,.005,.01,.015,.02,.025,.03,		&
                 .035,.04/
      DATA PSISCA/-15.,-10.,-6.,-3.,-1.5,-1.,-.6,-.4,-.2,-.1,0./
      DATA VNOSCA/1.0E-25,.01,.02,.03,.04,.05,.06,.07,.08,.09,.1/
      DATA CAPSCA/0.0,.05,.1,.15,.2,.25,.3,.35,.4,.45,.5/
      DATA TEMPSCA/0.00,5.00,10.00,15.00,20.00,25.00,30.00,35.00,		&
                  40.00,45.00,50.00/ 
      DATA TTL1R  /'ROOTS IN EACH CELL, TOTAL               '/
      DATA TTL2R  /'                                        '/
      DATA TTL1   /'VOLUMETRIC WATER CONTENT OF SOIL        '/
      DATA TTL3   /'                                        '/
      DATA TTL4   /'SOIL WATER POTENTIAL FOR          '/
      DATA TTL5   /'VOLUMETRIC NITRATE CONTENT OF SOIL      '/
      DATA TTL7   /'VOLUMETRIC AMMONIA CONTENT OF SOIL      '/
      DATA TTL8   /'AVERAGE SOIL TEMPERATURE                '/
      DATA TTL9   /'AT THE END OF THE DAY                   '/
      DATA UNITS  /'G/CM**3 SOIL            '/
      DATA UNITPS /' BARS           '/
      DATA PSIUNI /' BARS IN ROOT ZONE      '/
      DATA VNOUNI /' MG/N PER CM**3         '/
      DATA VH2UNI /' CM**3/CM**3 SOIL       '/
      DATA TUNITS /' DEGREES CELSIUS        '/
      DATA UNITSR /' GM. DRY WEIGHT '/
      DATA UNITST /' MM WATER       '/
      DATA NITUNT /' LBS N PER ACRE '/
      DATA TEMPUNT/' DEGREES CELSIUS'/
 
 
 100  FORMAT(///6X,A40,15X,'DAY ',I2.2,'/',I2.2/6X,A40/					 &
        6X,'UNITS - ',A24,23X,'LEGEND'//24X,'1 1 1 1 1 1 1 1 1 1 2'/	 &
        6X,'1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0',18X,A1,			 &
           '  <= ',F8.4//53X,F8.4,' < ',A1,' <= ',F8.4)
 102  FORMAT(3(1X,I2,3X,20A2/),1X,I2,3X,20A2,7X,F8.4,' < ',				 &
             A1,' <= ',F8.4)
 104  FORMAT(3(1X,I2,3X,20A2/),1X,I2,3X,20A2,7X,F8.4,' < ',A1//			 &
             6X,'TOTAL = ',F11.2,1X,A16)
 200 FORMAT(///6X,A40,15X,'DAY ',I2.2,'/',I2.2/6X,A40/					 &
        6X,'UNITS - ',A24,23X,'LEGEND'//24X,'1 1 1 1 1 1 1 1 1 1 2'/	 &
        6X,'1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0',18X,A1,			 &
           '  < ',E8.2//53X,E8.2,' < ',A1,' < ',E8.2)
 202  FORMAT(3(1X,I2,3X,20A2/),1X,I2,3X,20A2,7X,E8.2,' < ',				 &
            A1,' < ',E8.2)
 204  FORMAT(3(1X,I2,3X,20A2/),1X,I2,3X,20A2,7X,E8.2,' < ',A1//			 &
               6X,'TOTAL = ',F11.2,1X,A16)
      IF(IGO.EQ.1) THEN
        IF(TOTAL.EQ.TNNH4) THEN
           TL1=TTL7
        ELSE
           TL1=TTL5
        ENDIF
        TL2=TTL3
        UNTS=VNOUNI
        UNTST=NITUNT
        DO 20 I=1,11
          RANGE(I)=VNOSCA(I)
   20   CONTINUE
      ENDIF
      IF(IGO.EQ.2) THEN
        TL1=TTL1
        TL2=TTL3
        UNTS=VH2UNI
        UNTST=UNITST
        DO 40 I=1,11
          RANGE(I)=CAPSCA(I)
   40   CONTINUE
      ENDIF
      IF(IGO.EQ.3) THEN
        TL1=TTL1R
        TL2=TTL2R
        UNTS=UNITS
        UNTST=UNITSR
        DO 60 I=1,11
          RANGE(I)=ROOSCA(I)
   60   CONTINUE
      ENDIF
      IF(IGO.EQ.4) THEN
        TL1=TTL4
        TL2=TTL3
        UNTS=PSIUNI
        UNTST=UNITPS
        DO 80 I=1,11
          RANGE(I)=PSISCA(I)
   80   CONTINUE
      ENDIF
      IF (IGO.EQ.5) THEN
         TL1 = TTL8
         TL2 = TTL9
         UNTS = TUNITS
         UNTST=TEMPUNT
         DO 70 I=1,11
            RANGE(I)=TEMPSCA(I)
   70    CONTINUE
      ENDIF     
      DO 1 K=1, 20
      DO 1 L=1, 40
      ARAYLK = ARRAY(L,K)
      DO 2 I=1, 11
      RANGE1 = RANGE(I)
      IF(ARAYLK.LE.RANGE1) GO TO 1
 2    CONTINUE
      I = 12
 1    KHAR(L,K) = KA(I)
      IF(IGO.EQ.2) GO TO 15
      RANGE1 = RANGE(1)
      WRITE(24,100) TL1,MO,DAZE,TL2,UNTS,KA(1),RANGE1,RANGE1,KA(2),		 &
                    RANGE(2)
      INDX=0
      DO 14 L=1, 33, 4
      INDX=INDX+1
      L1 = L+1
      L2=L+2
      L3=L+3
 14   WRITE(24,102)L,(KHAR(L,K),K=1,20),L1,(KHAR(L+1,K),K=1,20),		 &
                   L2,(KHAR(L+2,K),K=1,20),L3,(KHAR(L+3,K),K=1,20),		 &
                   RANGE(INDX+1),KA(INDX+2),RANGE(INDX+2)
      L37=37
      L38=38
      L39=39
      L40=40
      WRITE(24,104) L37,(KHAR(37,K),K=1,20),L38,(KHAR(38,K),K=1,20),	 &
                    L39,(KHAR(39,K),K=1,20),L40,(KHAR(40,K),K=1,20),	 &
                    RANGE(11),KA(12),TOTAL,UNTST
      RETURN
 15   RANGE1 = RANGE(1)
      WRITE(24,200) TL1,MO,DAZE,TL2,UNTS,KA(1),RANGE1,RANGE1,KA(2),		 &
                    RANGE(2)
      INDX=0
      DO 16 L=1, 33, 4
      L1 = L+1
      L2=L+2
      L3=L+3
      INDX=INDX+1
 16   WRITE(24,202)L,(KHAR(L,K),K=1,20),L1,(KHAR(L+1,K),K=1,20),		 &
                   L2,(KHAR(L+2,K),K=1,20),L3,(KHAR(L+3,K),K=1,20),		 &
                   RANGE(INDX+1),KA(INDX+2),RANGE(INDX+2)
      L37 = 37
      L38=38
      L39=39
      L40=40
      WRITE(24,204) L37,(KHAR(L37,K),K=1,20),L38,(KHAR(L38,K),K=1,20),	  &
                    L39,(KHAR(L39,K),K=1,20),L40,(KHAR(L40,K),K=1,20),	  &
                    RANGE(11),KA(12),TOTAL,UNTST
      RETURN
      END


      SUBROUTINE MDSEA(I)

      use common_block

      CHARACTER PMAPFLE*12,FILEPATH*80
      DIMENSION IDUMCDE(40,15)

 1000 FORMAT(A16,F6.1,2X,F4.2)
 1020 FORMAT(15I5)
 1040 FORMAT(I2,1X,I2,1X,I2)
 1060 FORMAT(2F10.1,3I10)
 1080 FORMAT(6I10)
 
      IF(KULKNT.EQ.1) THEN
        CALL NUMCHR(PRONAM,J)
        IF(I.GT.1) THEN
          J = J+2
          IF(J.GT.8) J = 8
          K = I
          L = I/10
          IF(I.GT.9) K = I - L*10
          PMAPFLE(1:1) = CHAR(L+48)
          PMAPFLE(2:2) = CHAR(K+48)
          PMAPFLE(3:J) = PRONAM(1:J-2)
        ELSE
          PMAPFLE(1:J) = PRONAM(1:J)
        ENDIF
        PMAPFLE(J+1:J+4) = '.MAP'
        FILEPATH(1:1) = CHAR(92)
        FILEPATH(2:5) = 'PMAP'
        FILEPATH(6:6) = CHAR(92)
        FILEPATH(7:J+10) = PMAPFLE(1:J+4)
        OPEN(11,FILE=FILEPATH,STATUS='UNKNOWN',ERR=120)
        DUM = 0.0
        DUM0 = 3.0
        WRITE(11,1000) PRONAM,DUM,DUM0
        WRITE(11,1040) MO,DAZE,IYEAR
        DO II=1,50
          DO KK=1,mxvbrch
            IF(FCODE(KK,1,1).GT.0) THEN
            DO L=1,mxfbrch
              DO M=1,mxfsite
                IDUMCDE(L,M) = FCODE(KK,L,M)
                IF(IDUMCDE(L,M).EQ.2) IDUMCDE(L,M) = 3
                IF(IDUMCDE(L,M).EQ.7) IDUMCDE(L,M) = 2
                IF((IDUMCDE(L,M).GT.0).AND.						 &
                   (FFRUT(KK,L,M)*100..LT.II*2-1)) 				 &
                                                IDUMCDE(L,M) = 4
              ENDDO
            ENDDO
            WRITE(11,1000) 'REP 1           '
            DUM0 = II+(KK-1)/10.
            IF(KK.GT.1) THEN
              IDUM1 = 0
            ELSE
              IDUM1 = NUMPFN
            ENDIF
            IF(NFBR(3).GT.0) THEN
              IDUM2 = 2
            ELSEIF(NFBR(2).GT.0) THEN
              IDUM2 = 1
            ELSE
              IDUM2 = 0
            ENDIF
            DUM1 = Z/2.54
            WRITE(11,1060) DUM0,DUM1,NFBR(KK),IDUM1,IDUM2
            IDUM1 = 0
            WRITE(11,1080) IDUM1,IDUM1,IDUM1,IDUM1,IDUM1
            WRITE(11,1080) IDUM1,IDUM1,IDUM1,IDUM1,IDUM1,IDUM1
            WRITE(11,1080) IDUM1,IDUM1,IDUM1,IDUM1,IDUM1
            WRITE(11,1020) (IDUMCDE(1,M),M=1,5),IDUM1,					&
               (IDUMCDE(2,M),M=1,5),IDUM1,(IDUMCDE(3,M),M=1,3)
            WRITE(11,1020) (IDUMCDE(3,M),M=4,5),IDUM1,					&
               (IDUMCDE(4,M),M=1,5),(IDUMCDE(5,M),M=1,5),				&
               (IDUMCDE(6,M),M=1,2)
            WRITE(11,1020) (IDUMCDE(6,M),M=3,5),(IDUMCDE(7,M),M=1,5),	&
               (IDUMCDE(8,M),M=1,5),(IDUMCDE(9,M),M=1,2)
            WRITE(11,1020) (IDUMCDE(9,M),M=3,4),(IDUMCDE(10,M),M=1,4),	&
               (IDUMCDE(11,M),M=1,4),(IDUMCDE(12,M),M=1,4),				&
                IDUMCDE(13,1)
            WRITE(11,1020) (IDUMCDE(13,M),M=2,4),(IDUMCDE(14,M),M=1,3),	&
               (IDUMCDE(15,M),M=1,3),(IDUMCDE(16,M),M=1,3),				&
               (IDUMCDE(17,M),M=1,3)
            WRITE(11,1020) (IDUMCDE(18,M),M=1,3),(IDUMCDE(19,M),M=1,2),	&
               (IDUMCDE(20,M),M=1,2),(IDUMCDE(21,M),M=1,2),				&
               (IDUMCDE(22,M),M=1,2),(IDUMCDE(23,M),M=1,2),				&
                IDUMCDE(24,1),IDUMCDE(25,1)
            ENDIF
          ENDDO
        ENDDO
        WRITE(11,1000) 'END             '
        CLOSE(11)
      ENDIF
  120 CONTINUE
      Z = PLTHT(I)*2.54
      IF(MSATYP(I).GT.0) THEN
        DO 300 K=1,mxvbrch
          DO 301 L=1,mxfbrch
            DO 302 M=1,mxfsite
              IF(FCODE(K,L,M).GT.0) THEN
                J=K*150-L*5+M
                FFRUT(K,L,M)=AVGPLT(I,J)
              ENDIF
  302		continue
  301	  continue
  300   CONTINUE
      ENDIF
      IF((RUNMODE.EQ.'pmap').OR.(RUNMODE.EQ.'PMAP')) THEN
        DO 320 K=1,mxvbrch
          DO 320 L=1,mxfbrch
            DO 320 M=1,mxfsite
              FRUTP(K,L,M)=FFRUT(K,L,M)*10.+.001
              IF(FRUTP(K,L,M).LT.1.) FRUTP(K,L,M)=10.1
  320   CONTINUE
        CALL COTPLT(3)
      ELSE
        CALL PMAPS
!        IF(IOUTFG(13).GT.0) THEN  ' BY GR
          CALL COTPLT(2)
          CALL COTPLT(3)
 !       ENDIF  ' BY GR
      ENDIF
      RETURN
      END


      SUBROUTINE PREPRINT

      use common_block

      INTEGER*2 IDUM

 1140 FORMAT(' FIRST SQUARE  ',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1200 FORMAT(' LAST ACT WTHER',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1220 FORMAT(' FIRST BLOOM   ',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1700 FORMAT(A80)
 1760 FORMAT(' 1ST OPN BOLL  ',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1780 FORMAT(' 60% OPN BOLL  ',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1800 FORMAT(' MAX YIELD-100#',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1820 FORMAT(' MAX YIELD- 50#',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1840 FORMAT(' MAX YIELD     ',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1860 FORMAT(' NITR. STRESS *',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1880 FORMAT(' WATER STRESS *',2X,I2.2,'/',I2.2,I4,F5.1,I6,		 &
             F5.1,F6.1,I6,2I7,F6.2)
 1900 FORMAT(/ '* First simulated water or nitrogen stress after ',	   &
      'last day of actual weather.' //)

      IF(ABEND) RETURN
      IFG1=0
      IFG2=0
      IFG3=0
      IFG4=0
      IFG5=0
      IFG6=0
      IFG7=0
      YMX100=YIELD*500.-100.
      YMX50=YIELD*500.-50.
      YLDMAX=YIELD
      REWIND(21)
  400 CONTINUE
      READ(21,END=420,ERR=500) IDAY,KDAY,YIELD,NF,WSTRSD,Z,RI,TMAX			&
            ,TMIN,SQRZ,RAIN,FERN,GBZ2,PSILD,NV,LAI,SITEZ,NOPEN,FLOSS		&
            ,ABZ,CSTRES,PSIAVG,TAVG,TDAY,TNYT,C01,WIND,ES,EP				&
            ,RN,PN,CUMES,CUMEP,CUMSOK,CUMRAN,TH2O,H2OBAL,SPN,PLANTW		    &
            ,CHOBAL,LEAFWT,STEMWT,ROOTWT,SQWT,COTXX,GBOLWT,DEADWT			&
            ,XTRAC,POPPLT,LEAFCN,ISQ,DAZE,MO,N00,PTSRED,NUMPFN				&
            ,SDWSTM,SDWLEF,RCH2O,SDWBOL,AVAILN,SUPNO3,SUPNH4				&
            ,LEFABS,RUNOFF(DAYNUM),AMTIRR(DAYNUM)
      I00=SQRZ*POPPLT/1000+.5
      I01=GBZ2*POPPLT/1000+.5
      I02=NOPEN*POPPLT/1000+.5
      IDUM=KDAY+EMERGE-1

      IF((KDAY.GT.0).AND.(KDAY.EQ.ISQ)) THEN
	   WRITE(22,1140) MO,DAZE,KDAY,Z,N00,LAI,C01,I00,I01,I02,YIELD
	   WRITE(LSTNG,1140) MO,DAZE,KDAY,Z,N00,LAI,C01,I00,I01,I02,YIELD
	endif
      IF((KDAY.GT.0).AND.(KDAY.EQ.FBLOOM)) THEN 
	   WRITE(22,1220) MO,DAZE,KDAY,Z,N00,LAI,C01,I00,I01,I02,YIELD
	   WRITE(LSTNG,1220) MO,DAZE,KDAY,Z,N00,LAI,C01,I00,I01,I02,YIELD
	endif
      IF(DAYNUM.EQ.LDAYAW) THEN
	   WRITE(22,1200) MO,DAZE,KDAY,Z,N00,LAI,C01,I00,I01,I02,YIELD
	   WRITE(LSTNG,1200) MO,DAZE,KDAY,Z,N00,LAI,C01,I00,I01,I02,YIELD
	endif

      IF(IDUM.GT.LDAYAW) THEN
         IF((IFG1.EQ.0).AND.(PSIAVG.LE.PSICMX)) THEN
            WRITE(22,1880,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                          I00,I01,I02,YIELD
            IFG1=1
            WRITE(lstng,1880,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,		   &
                          I00,I01,I02,YIELD
         ENDIF


         IF((GBZ2+NOPEN).GT.0.) THEN
            IF((IFG2.EQ.0).AND.(NF.LT.1.)) THEN
               WRITE(22,1860,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,		   &
                          I00,I01,I02,YIELD
               IFG2=1
               WRITE(lstng,1860,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,	   &
                          I00,I01,I02,YIELD
            ENDIF
         ELSE
            IF((IFG2.EQ.0).AND.(NV.LT.1.).AND.(KDAY.GT.0)) THEN			   
               WRITE(22,1860,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,		   &
                          I00,I01,I02,YIELD
               IFG2=1
               WRITE(lstng,1860,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,	   &
                          I00,I01,I02,YIELD
             ENDIF
         ENDIF
      ENDIF
      IF((IFG3.EQ.0).AND.(NOPEN.GT..1)) THEN
          WRITE(22,1760,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                           I00,I01,I02,YIELD
          IFG3=1
          WRITE(lstng,1760,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                           I00,I01,I02,YIELD
        ENDIF
        SIXPCT=(NOPEN+GBZ2)*.6
        IF((IFG4.EQ.0).AND.(NOPEN.GE.SIXPCT).AND.(NOPEN.GT.0.)) THEN
          WRITE(22,1780,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                          I00,I01,I02,YIELD
          IFG4=1
          WRITE(lstng,1780,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                          I00,I01,I02,YIELD
        ENDIF
        IF((IFG5.EQ.0).AND.(YMX100.LE.YIELD*500.).AND.					   &
                                                (YIELD.GT.0.)) THEN
          WRITE(22,1800,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                          I00,I01,I02,YIELD
          IFG5=1
          WRITE(lstng,1800,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                          I00,I01,I02,YIELD
        ENDIF
        IF((IFG6.EQ.0).AND.(YMX50.LE.YIELD*500).AND.					   &
                                                (YIELD.GT.0.)) THEN
          WRITE(22,1820,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                          I00,I01,I02,YIELD
          IFG6=1
          WRITE(lstng,1820,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                          I00,I01,I02,YIELD
        ENDIF
        IF((IFG7.EQ.0).AND.(YIELD.GE.YLDMAX).AND.						   &
                                                (YIELD.GT.0.)) THEN
          WRITE(22,1840,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                          I00,I01,I02,YIELD
          IFG7=1
          WRITE(lstng,1840,ERR=520) MO,DAZE,KDAY,Z,N00,LAI,C01,			   &
                          I00,I01,I02,YIELD
        ENDIF
      GO TO 400
  420 CONTINUE
      WRITE(22,1900,ERR=520)
      WRITE(LSTNG,1900)
      RETURN
  500 CONTINUE
        PRINTBUF=' Error reading binary output file in preprint.'
        CALL WRITERR
      RETURN
  520 CONTINUE
        PRINTBUF=' Error writing summary file.'
        CALL WRITERR
      RETURN
      END


      SUBROUTINE PRINTOUT
!  ****************************************************************
!  *   COLLECTION OF WRITE STATEMENTS FOR OUTPUT OF MODEL RESULTS *
!  ****************************************************************

      use common_block

 1000 FORMAT(' ** FIRST SQUARE ON ',I2.2,'/',I2.2,' **')
 1020 FORMAT(11X,' YIELD BALES/ACRE   YIELD LBS/ACRE')
 1030 FORMAT(17X,F4.2,13X,F7.0)
 1040 FORMAT(' ************ FIRST BLOOM ON DAY ',I3,' *************')
 1200 FORMAT('TABLE 1:',17X,'GENERAL OUTPUT IN ENGLISH UNITS')
 1220 FORMAT(17X,'PER/PLANT',27X,'PER/ACRE')
 1240 FORMAT(12X,18('-'),4X,43('-'))								 
 1260 FORMAT(' DAE  DATE  HEIGHT  LAI  NODES   SITES',			 &
        '  SQUARES   GREEN    OPEN  YIELD  FRUIT')
 1280 FORMAT(13X,'(IN)',33X,'BOLLS   BOLLS',9X,' SHED')
 1290 FORMAT(1X,I3,1X,I2,'/',I2.2,F7.1,F6.2,I5,2X,I8,1X,3I8,F6.2,I8)
 1300 FORMAT('TABLE 2:',22X,'STRESS FACTORS OUTPUT')
 1320 FORMAT('          %                         LEAF         ',	 &
       '                LEAF    SOIL')
 1340 FORMAT('  DATE LIGHT  CSTRES  NITR-STRESS    N    WATER ',	 &
       '   ET  PIXCON   WATER   WATER')
 1360 FORMAT('        INT           FRUIT   VEG  CONC  STRESS ',	 &
       '  (in)  (ppm)    POT     POT')
 1380 FORMAT(1X,I2.2,'/',I2.2,F6.1,F7.3,2F7.2,F6.3,F8.3,F6.3,F8.3,	 &
            2F8.3)
 1400 FORMAT('TABLE 5:',15X,'WEATHER VARIABLES IN METRIC UNITS')
 1420 FORMAT(12X,'TEMPERATURE DEGREES (C)')
 1440 FORMAT('                    DAILY  DAY  NIGHT (LNLYS)   (MM)',  &
       '   (MM)     (MM)')
 1460 FORMAT(12X,'VARIETY PARAMETERS FOR ',A9)
 1480 FORMAT(5X,5E12.4)
 1500 FORMAT('TABLE 5:',15X,'WEATHER VARIABLES IN ENGLISH UNITS')
 1520 FORMAT(12X,'TEMPERATURE DEGREES (F)')
 1540 FORMAT(9X,28('-'))
 1560 FORMAT('  DATE   MAX   MIN   AVG   AVG   AVG   SOLAR    RAIN',  &
       '    IRR    RUNOFF  WIND') 
 1580 FORMAT('                    DAILY  DAY  NIGHT (LNLYS)   (IN)',  &
       '   (IN)     (IN)')
 1590 FORMAT(1X,I2.2,'/',I2.2,5F6.1,5F8.2)
 1600 FORMAT('TABLE 4:',21X,'WATER AND ET VARIABLES')
 1620 FORMAT('  DATE    ES    EP     RN     PN  CUMES  CUMEP',	  &
       ' CUMSOK  CUMRAN    TH2O  H2OBAL')
 1640 FORMAT(1X,I2.2,'/',I2.2,2F6.2,F7.1,F7.3,3F7.1,3F8.1)
 1700 FORMAT('TABLE 3:',21X,'WEIGHTS IN METRIC UNITS')
 1720 FORMAT('  DATE    SPN PLANTW CO2BL  LFWT STMWT RUTWT',		  &
       ' SQRWT  COTXX GBOLWT DEADWT XTRAC')
 1740 FORMAT(1X,I2.2,'/',I2.2,2F7.2,F6.2,3F6.2,F6.2,3F7.2,F6.1)
 1800 FORMAT(' ')
 1820 FORMAT('1')
 1840 FORMAT(12X,'SOIL PARAMETERS FOR ',A13)
 1860 FORMAT(1X,I5,5E12.3)
 1880 FORMAT(6X,3E12.3,2I12)

      IF(IOUTFG(1).GT.0) THEN
         WRITE(lstng,1800)
         WRITE(lstng,1800)
         WRITE(lstng,1200)
         WRITE(lstng,1220)
         WRITE(lstng,1240)
         WRITE(lstng,1260)
         WRITE(lstng,1280)
         WRITE(lstng,1800)
         REWIND(21)
  180    CONTINUE
         READ(21,END=200,ERR=200) IDAY,KDAY,YIELD,NF,WSTRSD,Z,RI,TMAX	 &
            ,TMIN,SQRZ,RAIN,FERN,GBZ2,PSILD,NV,LAI,SITEZ,NOPEN,FLOSS	 &
            ,ABZ,CSTRES,PSIAVG,TAVG,TDAY,TNYT,C01,WIND,ES,EP			 &
            ,RN,PN,CUMES,CUMEP,CUMSOK,CUMRAN,TH2O,H2OBAL,SPN,PLANTW	     &
            ,CHOBAL,LEAFWT,STEMWT,ROOTWT,SQWT,COTXX,GBOLWT,DEADWT		 &
            ,XTRAC,POPPLT,LEAFCN,ISQ,DAZE,MO,N00,PIXCON,NUMPFN			 &
            ,SDWSTM,SDWLEF,RCH2O,SDWBOL,AVAILN,SUPNO3,SUPNH4			 &
            ,LEFABS,RUNOFF(DAYNUM),AMTIRR(DAYNUM)
         IF(KDAY.GT.0) THEN
            IF(KDAY.EQ.ISQ) WRITE(lstng,1000) MO,DAZE
            N01=IFIX(SITEZ*POPPLT)
            N02=IFIX(SQRZ*POPPLT)
            N03=IFIX(GBZ2*POPPLT)
            N04=IFIX(NOPEN*POPPLT)
            N05=IFIX(ABZ*POPPLT)
	        n06=IFIX(ABZ0*POPPLT)
            WRITE(lstng,1290) KDAY,MO,DAZE,Z,LAI,N00,					 &
                            N01,N02,N03,N04,YIELD,N05
         ENDIF
         GO TO 180
  200    CONTINUE
         YLBSAC=YIELD*500.
         WRITE(lstng,1800)
         WRITE(lstng,1020)
         WRITE(lstng,1030) YIELD,YLBSAC
         WRITE(lstng,1800)
         WRITE(lstng,1800)
	endif

      IF(IOUTFG(2).GT.0) THEN
         WRITE(lstng,1300)
         WRITE(lstng,1320)
         WRITE(lstng,1340)
         WRITE(lstng,1360)
         WRITE(lstng,1800)
         REWIND(21)
  220    CONTINUE
         READ(21,END=240,ERR=240) IDAY,KDAY,YIELD,NF,WSTRSD,Z,RI,TMAX	  &
            ,TMIN,SQRZ,RAIN,FERN,GBZ2,PSILD,NV,LAI,SITEZ,NOPEN,FLOSS	  &
            ,ABZ,CSTRES,PSIAVG,TAVG,TDAY,TNYT,C01,WIND,ES,EP			  &
            ,RN,PN,CUMES,CUMEP,CUMSOK,CUMRAN,TH2O,H2OBAL,SPN,PLANTW	      &
            ,CHOBAL,LEAFWT,STEMWT,ROOTWT,SQWT,COTXX,GBOLWT,DEADWT		  &
            ,XTRAC,POPPLT,LEAFCN,ISQ,DAZE,MO,N00,PIXCON,NUMPFN			  &
            ,SDWSTM,SDWLEF,RCH2O,SDWBOL,AVAILN,SUPNO3,SUPNH4			  &
            ,LEFABS,RUNOFF(DAYNUM),AMTIRR(DAYNUM)
         IF(KDAY.GT.0) THEN
            XDUM=(ES+EP)/25.4
            IF(KDAY.EQ.ISQ) WRITE(lstng,1000) MO,DAZE
            WRITE(lstng,1380) MO,DAZE,C01,CSTRES,NF,NV,LEAFCN,			  &
                           WSTRSD,XDUM,PIXCON,PSILD,PSIAVG
         ENDIF
         GO TO 220
  240    CONTINUE
         YLBSAC=YIELD*500.
         WRITE(lstng,1800)
         WRITE(lstng,1020)
         WRITE(lstng,1030) YIELD,YLBSAC
         WRITE(lstng,1800)
         WRITE(lstng,1800)
	endif

      IF(IOUTFG(3).GT.0) THEN
         WRITE(lstng,1700)
         WRITE(lstng,1720)
         WRITE(lstng,1800)
         REWIND(21)
  260    CONTINUE
         READ(21,END=280,ERR=280) IDAY,KDAY,YIELD,NF,WSTRSD,Z,RI,TMAX	  &
            ,TMIN,SQRZ,RAIN,FERN,GBZ2,PSILD,NV,LAI,SITEZ,NOPEN,FLOSS	  &
            ,ABZ,CSTRES,PSIAVG,TAVG,TDAY,TNYT,C01,WIND,ES,EP			  &
            ,RN,PN,CUMES,CUMEP,CUMSOK,CUMRAN,TH2O,H2OBAL,SPN,PLANTW	      &
            ,CHOBAL,LEAFWT,STEMWT,ROOTWT,SQWT,COTXX,GBOLWT,DEADWT		  &
            ,XTRAC,POPPLT,LEAFCN,ISQ,DAZE,MO,N00,PIXCON,NUMPFN			  &
            ,SDWSTM,SDWLEF,RCH2O,SDWBOL,AVAILN,SUPNO3,SUPNH4			  &
            ,LEFABS,RUNOFF(DAYNUM),AMTIRR(DAYNUM)
         IF(KDAY.GT.0) THEN
            IF(KDAY.EQ.ISQ) WRITE(lstng,1000) MO,DAZE
            WRITE(lstng,1740) MO,DAZE,SPN,PLANTW,CHOBAL,LEAFWT,			  &
                       STEMWT,ROOTWT,SQWT,COTXX,GBOLWT,DEADWT,XTRAC
         ENDIF
         GO TO 260
  280    CONTINUE
         YLBSAC=YIELD*500.
         WRITE(lstng,1800)
         WRITE(lstng,1020)
         WRITE(lstng,1030) YIELD,YLBSAC
         WRITE(lstng,1800)
         WRITE(lstng,1800)
	endif

      IF(IOUTFG(4).GT.0) THEN
         WRITE(lstng,1600)
         WRITE(lstng,1620)
         WRITE(lstng,1800)
         REWIND(21)
  300    CONTINUE
         READ(21,END=320,ERR=320) IDAY,KDAY,YIELD,NF,WSTRSD,Z,RI,TMAX	  &
            ,TMIN,SQRZ,RAIN,FERN,GBZ2,PSILD,NV,LAI,SITEZ,NOPEN,FLOSS	  &
            ,ABZ,CSTRES,PSIAVG,TAVG,TDAY,TNYT,C01,WIND,ES,EP			  &
            ,RN,PN,CUMES,CUMEP,CUMSOK,CUMRAN,TH2O,H2OBAL,SPN,PLANTW	      &
            ,CHOBAL,LEAFWT,STEMWT,ROOTWT,SQWT,COTXX,GBOLWT,DEADWT		  &
            ,XTRAC,POPPLT,LEAFCN,ISQ,DAZE,MO,N00,PIXCON,NUMPFN			  &
            ,SDWSTM,SDWLEF,RCH2O,SDWBOL,AVAILN,SUPNO3,SUPNH4			  &
            ,LEFABS,RUNOFF(DAYNUM),AMTIRR(DAYNUM)
         IF(KDAY.GT.0) THEN
            IF(KDAY.EQ.ISQ) WRITE(lstng,1000) MO,DAZE
            WRITE(lstng,1640) MO,DAZE,ES,EP,RN,PN,CUMES,CUMEP,			  &
                       CUMSOK,CUMRAN,TH2O,H2OBAL
         ENDIF
         GO TO 300
  320    CONTINUE
         YLBSAC=YIELD*500.
         WRITE(lstng,1800)
         WRITE(lstng,1020)
         WRITE(lstng,1030) YIELD,YLBSAC
         WRITE(lstng,1800)
         WRITE(lstng,1800)
	endif

      IF(IOUTFG(5).GT.0) THEN
         WRITE(lstng,1400)
         WRITE(lstng,1420)
         WRITE(lstng,1540)
         WRITE(lstng,1560)
         WRITE(lstng,1440)
         WRITE(lstng,1800)
         REWIND(21)
  340    CONTINUE
         READ(21,END=360,ERR=360) IDAY,KDAY,YIELD,NF,WSTRSD,Z,RI,TMAX	 &
            ,TMIN,SQRZ,RAIN,FERN,GBZ2,PSILD,NV,LAI,SITEZ,NOPEN,FLOSS	 &
            ,ABZ,CSTRES,PSIAVG,TAVG,TDAY,TNYT,C01,WIND,ES,EP			 &
            ,RN,PN,CUMES,CUMEP,CUMSOK,CUMRAN,TH2O,H2OBAL,SPN,PLANTW	     &
            ,CHOBAL,LEAFWT,STEMWT,ROOTWT,SQWT,COTXX,GBOLWT,DEADWT		 &
            ,XTRAC,POPPLT,LEAFCN,ISQ,DAZE,MO,N00,PIXCON,NUMPFN			 &
            ,SDWSTM,SDWLEF,RCH2O,SDWBOL,AVAILN,SUPNO3,SUPNH4			 &
            ,LEFABS,RUNOFF(DAYNUM),AMTIRR(DAYNUM)
         IF(KDAY.GT.0) THEN
            IF(KDAY.EQ.ISQ) WRITE(lstng,1000) MO,DAZE
            WIND=WIND*1.609344
            C00=RAIN+RUNOFF(DAYNUM)-AMTIRR(DAYNUM)*25.4
            C01=AMTIRR(DAYNUM)*25.4
            WRITE(lstng,1590) MO,DAZE,TMAX,TMIN,TAVG,TDAY,TNYT,			 &
                  RI,C00,C01,RUNOFF(DAYNUM),WIND
         ENDIF
         GO TO 340
  360    CONTINUE
         YLBSAC=YIELD*500.
         WRITE(lstng,1800)
         WRITE(lstng,1020)
         WRITE(lstng,1030) YIELD,YLBSAC
         WRITE(lstng,1800)
	endif

      IF(IOUTFG(6).GT.0) THEN
         WRITE(lstng,1500)
         WRITE(lstng,1520)
         WRITE(lstng,1540)
         WRITE(lstng,1560)
         WRITE(lstng,1580)
         WRITE(lstng,1800)
         REWIND(21)
  380    CONTINUE
         READ(21,END=400,ERR=400) IDAY,KDAY,YIELD,NF,WSTRSD,Z,RI,TMAX	 &
            ,TMIN,SQRZ,RAIN,FERN,GBZ2,PSILD,NV,LAI,SITEZ,NOPEN,FLOSS	 &
            ,ABZ,CSTRES,PSIAVG,TAVG,TDAY,TNYT,C01,WIND,ES,EP			 &
            ,RN,PN,CUMES,CUMEP,CUMSOK,CUMRAN,TH2O,H2OBAL,SPN,PLANTW		 &
            ,CHOBAL,LEAFWT,STEMWT,ROOTWT,SQWT,COTXX,GBOLWT,DEADWT		 &
            ,XTRAC,POPPLT,LEAFCN,ISQ,DAZE,MO,N00,PIXCON,NUMPFN			 &
            ,SDWSTM,SDWLEF,RCH2O,SDWBOL,AVAILN,SUPNO3,SUPNH4			 &
            ,LEFABS,RUNOFF(DAYNUM),AMTIRR(DAYNUM)
         IF(KDAY.GT.0) THEN
            IF(KDAY.EQ.ISQ) WRITE(lstng,1000) MO,DAZE
            TMAX=TMAX/5.*9.+32.
            TMIN=TMIN/5.*9.+32.
            TAVG=TAVG/5.*9.+32.
            TDAY=TDAY/5.*9.+32.
            TNYT=TNYT/5.*9.+32.
            RAIN=RAIN/25.4
            RUNOFF(DAYNUM)=RUNOFF(DAYNUM)/25.4
            C00=RAIN+RUNOFF(DAYNUM)-AMTIRR(DAYNUM)
            WRITE(lstng,1590) MO,DAZE,TMAX,TMIN,TAVG,TDAY,TNYT,			 &
                  RI,C00,AMTIRR(DAYNUM),RUNOFF(DAYNUM),WIND
         ENDIF
         GO TO 380
  400    CONTINUE
         YLBSAC=YIELD*500.
         WRITE(lstng,1800)
         WRITE(lstng,1020)
         WRITE(lstng,1030) YIELD,YLBSAC
         WRITE(lstng,1800)
      ENDIF
      write(lstng,*) "SPN  PN   ES"
      write(lstng,*) SPN, PN, ES
      write(lstng,*)
      WRITE(lstng,1840) hydfle
      WRITE(lstng,1800)
      DO 375 I=1,LYRSOL
         WRITE(lstng,1860) LDEPTH(I),DIFF0(I),THETA0(I),BETA(I),		 &
               THETAS(I),FCININ(I)
         WRITE(lstng,1880) THETAR(I),AIRDR(I),BD(I),IPSAND(I),			 &
               IPCLAY(I)
  375 CONTINUE
      WRITE(lstng,1800)
      WRITE(lstng,1800)

      call PrintSummary
 	  call preprint

  500 CONTINUE
	  return
	  end
 