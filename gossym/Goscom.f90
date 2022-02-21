! *** COMMON BLOCK FOR GOSSYM *** DATE MARCH 30,1984.
! *
! *** COLLECTION OF ALL FORTRAN TYPE DECLARATIONS FOR GOSSYM
! *
      module common_block

      DOUBLE PRECISION C1,DEC,PHI,DEGRAD,XLAT        
      LOGICAL ABEND,FULPRO,RTEXNT(40),SEND,SKPFLG,TUPF(41,21),TTUPF(41,21) 
      CHARACTER*1 PLTMAP(6000)
      CHARACTER*2 KA(12),KHAR(40,20),PRI(10,40),PRT(10,40,15)
      CHARACTER*4 CHAR1(7),CHAR2(7),CHAR3(13),CHAR4(13),CHARI
      CHARACTER*5 OPSYS
      CHARACTER*6 FMTHOD(7),HMTHOD(7),IMTHOD(7),PGRUNT(7),RUNMODE
      CHARACTER*7 VARNAME
 	  CHARACTER*8	VTYNAM
      CHARACTER*9 VERSION,VARITY(50)
      CHARACTER*10 MSDATE
      CHARACTER*11 FILFRM
      CHARACTER*12 SNAME
      CHARACTER*18 PROFLE,ACTwea,FURwea,IRRfle,HYDfle,INTfle,PGRfle, &
	               INSfle,FNGfle,PMAfle,gcmfle
      CHARACTER*15 RUNDATE
	  Character*25 yldfle,wtsfle
!!      CHARACTER*40 PROFILE,PRONAM,ACTWTH,FURWTH,PRDWTH,vldfle, &
!!                   IRRFRT,SOLHYD,INTSOL,PGRHRB,INSCDE,FNGCDE,PMAFIL, &
!!	               plotfle,listfle,sumryfle,tablefle,binfle
      CHARACTER*45 ROOTDIR
      CHARACTER*51 PDESCP
      CHARACTER*80 ERRFLE
      CHARACTER*80 PRINTBUF
      CHARACTER*120 PROFILE,PRONAM,ACTWTH,FURWTH,PRDWTH,vldfle,gcmfil, &
                   IRRFRT,SOLHYD,INTSOL,PGRHRB,INSCDE,FNGCDE,PMAFIL, &
	               plotfle,listfle,sumryfle,tablefle,binfle,solimpd, &
	               variety, plotfle2
      INTEGER*2 IOUTFG,IPCLAY,IPSAND,LATUDE,LDEPTH,LYRSOL,MFREQ, &
                MNDAYB,MTHIRR,NAPS,NFRQ,NUMRUN,SEASON,IFGIRR,IFGRAIN, &
                IVARTY,WATTBL,DEFMTH,PRPMTH,PIXMTH,NODPMAP, &
                NDMSMS,NDMSPF,NDV1MS,NDV1PF,NDV2MS,NDV2PF
      INTEGER CO2,DAYNUM,DAZE,DEFBGN,DEFDATE,DEFDAY, &
              EMERGE,FCODE,PIXDAY,PRPDATE,PRPDAY,POLYNA, dae
      REAL INT, KSAT, KSATC, KSATW, KWIDTH, LAGE, LAI, LAMDA, &
           LAMDAC, LAMDAS, LAREA, LEAFCN,LEAFRS, LEAFR1, LEAFW, &
		   LEAFWT, LEFABS, LMAX, LYTRES, &
           MH2O,MMUPN1,MMUPN2,MMUPN3,NDLAY,NEWEP,NEWES,MLAREA, &
           MLEAFW,NOPEN,NPOOL,NR,NV,NF,NYTTYM,NYTWTF,NLOSS,	ntop
 
      COMMON/ ARCOM /ABZ,ACELLDW,ACTIRRG,ACTRAIN,ADDEDN,ADPG(20), &
                     AGE(10,40,15),AGEABZ(10,40,15),AGEBOL(10,40,15), &
                     AGEPFN(10),AGETOP,AIRDR(9),AIRDRC,	&
                     AIRDRI,AIRDRW,ALPHA,AMTIRR(365),angboll(30),APRES, &
                     ARDRCN(40),AREA,AT,AVAILN,AVGPLT(30,6000),	&
                     AVGT(10,40,15),AVGTSD,AVGTSP,AVTEMP,AVTPFN(10), &
	                 add60,amicron,abzb,abz0
      COMMON/ BRCOM /BD(9),BDC,BDL(40),BDI,BDRATO,BDSLOP,BDW, &          
	                 BDELAY(10,40),BETAK(20), &
	                 BETA(9),BETAC,BETAI,BETAW,BOLL1,BOLOSS(366), &
	                 BOLTMP(10,40,15),BOLWGT(10,40,15),BURCN,BURMIN, &
	                 BURR1,BURRN,BSIZE(10,40,15),BLUM(366), bloom
      COMMON/ CICOM /CO2 
      COMMON/ CRCOM /CALBRT(60),CD,CDBOLL,CDLEAF,CDROOT, &
	                 CDSQAR,CDSTEM,CLIMAT(366,10), CHOBAL, &           
					 CMXIRR,CONAMM,CONNIT,CONSK(20),CONURA,COTXX, CPOOL, &           
					 CSTORE,CSTRES,CUMEP,CUMES,CUMRAN,CUMNSOK,CUMSOK, &           
					 cdlayf,cdlayv
      COMMON/ DICOM /DAYNUM,DAZE,DEFBGN,DEFDATE(5),DEFDAY
      COMMON/ DRCOM /DAY1PN,DAY1SN,DAYLNG,DAYTYM,DAYWTF, DCELL, &           
	                 DEAD2DAY,DEADWT,DEFKGH,DEFPPA(5),DEHISS(10,40,15), &           
					 DELAY(10,40),DELT,DIFF(40,20), DIFF0(9), &           
					 DIFF0C,DIFF0I,DIFF0W,DTAVG(7),dwrt(40,20),DZ, &
					 dae,dd60,d60avg(7),dayfb,daymt,daysnf,daysnv,daysq
      COMMON/ EICOM /EMERGE
      COMMON/ ERCOM /EP,EPAVG(5),ES,ETA(9),ETAC,ETAW
      COMMON/ FICOM /FCODE(10,40,15)
      COMMON/ FRCOM /F2,FBLOOM,FC(40),FCININ(9),FCINIC,FCFCTI, &               
	                 FCINIW,FERN,FFRUT(10,40,15),FL,FLOSS, &           
	                 FLNMIN,FLXMAX(9),FLXMIN(9),FNH4,FNL(40,21),FNO3, &           
					 FNU(41,20),FOLIARN,FRATIO,FRUTP(10,40,15),FS,FSQ, &           
					 FSTAVG(10,40,15),FSTRES,FWL(40,21),FWU(41,20),	&
					 fibwgt(10,40,15),flength(10,40,15),fmatur(10,40,15), &           
					 frstsq,frstbl,frstob, fsqfra,fblfra,fobfra
      COMMON/ GRCOM /GEOTR,GAMMA,GBLOS,GBOLWT,GBZ2,GH2OC(9),GIN,GINP, GSUBR
      COMMON/ HRCOM /H2OBAL,H2OINT(14)
      COMMON/ IICOM /IADJDY,IADJMO,IDAY,IMAP,INRIM,IPIX,IPLT,ISCRN,ISQ,ITER,IYEAR
      COMMON/ IRCOM /INT
      COMMON/ JICOM /JDAY,JDAYLW,JDSTRM,JDSTPM,JDSTRS,JDSTPS
      COMMON/ KICOM /KDAY,KRAIN,KLL(40),KRL(40),KULCLF,KULCRT, &           
	                 KULDAY(5),KULKNT,KUPPER,KUPPR2,KUPPR3,KUPPR4
      COMMON/ KRCOM /KSAT(9),KSATC,KSATW,KWIDTH
      COMMON/ LICOM /LAYVAL,LDAYAW,LDAYIR,LDAYFW,LDAYPW, &
	                 LEFCNT,LEFSRT(6000),LFATDF,LINE, LPLOW,LR,	&           
					 LSTNG,LTYPE,LVSLOS,LVS2LOS,LYRDPH(40)
      COMMON/ LRCOM /LAGE(10,40,15),LAI,LAMDA,LAMDAC,LAMDAS, &
                     LAREA(10,40,15),LEAFCN,LEAFR1,LEAFRS,LEAFW(10,40,15), &
                     LEAFWT,LEFABS,LMAX,LYTRES
      COMMON/ MICOM /MATURE(10,40,15),MCODE(10,40,15),MLDAY,MMAP,MO, &           
	                 MSADTE(30),MSANODE(30),MSATYP(30),MSDAY
      COMMON/ MRCOM /MH2O,MLAREA(10,40),MLEAFW(10,40), &           
	                 MMUPN1,MMUPN2,MMUPN3
      COMMON/ NICOM /NCURVE,NDAYS,NFBR(10),NFERT(365,7),NK,NL,NNOD(10,40), &
                     NOITR,NPLANT(30),NUMPFN,NVBRCH,n06
      COMMON/ NRCOM /NDLAY,NEWEP,NEWES,NF,NLOSS,NOPEN,NPOOL, &
                     NR,NV,NYTTYM,NYTWTF
      COMMON/ ORCOM /OMA(14),ORGN
      COMMON/ PICOM /PIXDAY(10),POLYNA,PRPDATE(5),PRPDAY
      COMMON/ PRCOM /PDADAY(10,40,15),PDANYT(10,40,15),PDAMLD(10,40), &          
	                 PDAMLN(10,40),PDBOLL,PDLEAF,pdwrt(40,20),PDROOT, &           
					 PDSQ,PDSTEM,PDWBOD(10,40,15),PDWBON(10,40,15),	&
					 PDWFLD(10,40,15),PDWFLN(10,40,15),PDWMLD(10,40), &           
					 PDWMLN(10,40),PDWSQ(10,40,15),PERDEF,PFAL(10),PFAREA, &
                     PFDAL,PFDALD(10),PFDALN(10),PFDWL,PFDWLD(10), &             
	                 PFDWLN(10),PFNODAGE(20),PFNODLTH(20), &          
					 PFWL(10),PI,PIN,PIXCON,PIXLOS,	& 
					 PIXPLT,PIXPPA(10),PIXDA,PIXDN,PIXDZ,PIXDPN,&           
					 PLANTN,PLANTW,PLEFABS,PLTHT(30),PLTN, &           
					 PLTPFT,PN,PNETCOR,PNBAL,POPFAC,POLOSS, &           
					 POPPLT,PQFLR,PRPKGH,PRPPPA(5),PSIAVG, &          
					 PSICMX,PSILD,PSILN,PSIMAX,PSINUM,PSIS(40,20),PSISAT, &           
					 PSISFC,PTSRED,PUPF(41,21),punits
      COMMON/ RRCOM /RAIN,RCH2O,RECDAT(24),REQ1,reqv,RESC,RESN, &          
	                 RI,RN,RNNH4(14),RNNO3(14),ROOTCN,ROOTN,ROOTR1,ROOTRS, &           
					 ROOTS,ROOTSV(40,20),ROOTWT,ROWSAVE,ROWSP,RSUBO, &           
					 RTIMPD(40,20),RTP1,RTP2,RTWT(40,20,3),RTWTCU(40,20), &           
					 RUNOFF(366),RUTOFF,rtwtcg(40,20),rnfactor
      COMMON/ SRCOM /SBOLL,SDWBOL,SDWLEF,SDWSQR, &           
	                 SDWSTM,SEEDCN,SEEDN,SEEDR1, &           
					 SESI,SESII,SITES,SITEZ,SKIPWD,	&           
					 SLEAF,SLEAFN,SLF,SOAKN(20), SOILT(40,20), &           
					 SPDWBO,SPDWLD,SPDWLN,SPDWRT,SPDWSQ,SPN,SOILN,SNBAL, &           
					 SQLOSS(366),SQUAR,SQRWT(10,40,15),SQRZ,SQWT,SROOT,	&           
					 SSTEM,STEMCN,STEMN,STEMRS,stemr1,STEMWT, &           
	                 STMWT(366),SUMES,SUMEP,SUMSTRS,SUPNO3,SUPNH4, &           
					 SUMSUB,SUBIRR,SUPF,str01,str02,str03,str04,str05,str06
      COMMON/ TRCOM /T,TAIR(15),TAVG,TCELL,TD,TDAY,TDFKGH,TEMP1C,TEMP1G, &
                     TEMP1R,TH2O,THETA0(9),THETAI,THETAR(9),THETAS(9), &
                     THTS(40),THTR(40),THAD(40),THRLN,THTA0C,THTA0I, &
                     THTA0W,THTARC,THTARI,THTARW,THTASC,THTASI,THTASW, &
                     TIH2OC,TMAX,TMIN,TNNH4,TNNO3,TNYT,TSMN(40),TSMX(40), &
                     TSOILD(40),TSOILN(40),TSOLAV(20),TSQ,TSTBD(9,40), &
                     TSTIMP(9,40),TNO3UP,TNH4UP
      COMMON/ URCOM /UPNO3,UPNH4,UPTAKEN,UPF(40,20)
      COMMON/ VRCOM /VCELL,VDELAY(10),VEGWT,VH2OC(40,20), &           
	                 VNC(40,20),VNH4C(40,20),VNO3C(40,20),VSTRES
      COMMON/ WRCOM /WATTSM,WCELL,WIND,WSTRLF,WSTRST,WSTRS,	&           
	                 WSTRSD,WSTRSN,WTDAY1,WTSLFD
      COMMON/ XRCOM /XMNODAGE(40),XMNODLTH(40),XMXIRR,XTRAC,XTRAN
      COMMON/ YRCOM /YIELD
      COMMON/ ZRCOM /Z,ZPIXD,ZUPT(40,20)
      COMMON/ CHCOM /ACTWTH,CHAR1,CHAR2,CHAR3,CHAR4,CHARI,ERRFLE,FILFRM, &
	                 FMTHOD,FNGCDE,FURWTH,profile,HMTHOD,IMTHOD,IRRFRT, &
					 INTSOL,INSCDE,KA,KHAR,MSDATE,OPSYS,PDESCP,PGRHRB, &
					 PGRUNT,PLTMAP,PMAFIL,PRDWTH,PRI,PRINTBUF,PROFLE, &
					 PRONAM,PRT,ROOTDIR,RUNDATE,VTYNAM,	&
                     RUNMODE,SNAME,SOLHYD,VARNAME,VARITY,VERSION
      COMMON/ LGCOM /ABEND,FULPRO,RTEXNT,SEND,SKPFLG,TUPF,TTUPF 
      COMMON/ DPCOM /C1(9),DEC,DEGRAD,PHI,XLAT        
      COMMON/ I2COM /DEFMTH(5),IFGIRR,IFGRAIN,IOUTFG(23),IPCLAY(9), &
	                 IPSAND(9),IVARTY,LATUDE,LDEPTH(9),LYRSOL, &
                     MFREQ,MNDAYB,MTHIRR(365),NAPS,NDMSMS,NDMSPF, & 
	                 NDV1MS,NDV1PF,NDV2MS,NDV2PF,NFRQ,NODPMAP(30,6), &
					 NUMRUN,PIXMTH(10),PRPMTH(5),SEASON,WATTBL

	common/ sumcom/RUNNO,SUMYLD,SUMISQ,SUMFBL,SUMOBL,SUMMAT,SUMMIK, &
	               SUMHT,SUMLAI,ISUMNOD,SUMLFWT,SUMSTMWT,SUMSQRZ, &			
				   SUMGBOL,isqday,iflday,iobday,imtday
    common/ avglyr/ VH2OC_A(40)

      INTEGER FCODES(6000)
      DIMENSION BOLAGE(6000), FRUITS(6000), MCODES(6000)
      DIMENSION WTSQ(6000), WTBO(6000), SQAGE(6000)
 
      EQUIVALENCE (FCODES(1),FCODE(1,1,1)), (MCODES(1),MCODE(1,1,1))
      EQUIVALENCE (BOLAGE(1),AGEBOL(1,1,1)),(FRUITS(1),FFRUT(1,1,1))
      EQUIVALENCE (WTSQ(1),SQRWT(1,1,1)),   (WTBO(1),BOLWGT(1,1,1))
      EQUIVALENCE (SQAGE(1),AGE(1,1,1))

! *** Common block for rutgro   kit 1999

      integer*2 ilc,irc,kmod,kl1,km1,kp1,kr1,ld1,ldc,lp1

	common/ integ1/ilc,irc,kmod,km1,kp1,ld1,ldc,lp1, &             
	               mxfbrch,mxvbrch,mxfsite,mxfruts

! *** Common block for plantmaps   kit 02/01/2000

	common/  ticom/	pfti,vbti,xmsti,fbnti,flx,fsx

! *** Common block for plantheight   kit 10/4/1999

	common/ nitcom/	cntlfcn(40),sumlfcn(40), xnodage(40), &                
	                cntpflfcn(10), sumpflfcn(10), &               
					cntmslfcn(40), summslfcn(40), idanoflg

! *** Common block for leaf area development   kit 10/4/1999

	common/  lfcom/ pflfage(10),pflfarea(10),agepflf(10), &                
	                xmslfage(10,40),xmslfarea(10,40),agemslf(10,40), &                
					frlfage(10,40,15),frlfarea(10,40,15), &				
					agefrlf(10,40,15),areapflf,areamslf,areafblf, &                
					day_lfstress,eve_lfstress

! *** Common block for fruit development   kit 12/14/1999

	common/ frtcom/ day_expfac,eve_expfac,day_lotemp,eve_lotemp, &
	                day_hitemp,eve_hitemp,eve_water_index, &                
					day_water_index,Bloom_tavg(365),Boll_tavg, &				
					HeatIndex,susceptible_bolls
 
! *** Common block for stem development   kit 12/14/1999

	common/ stmcom/ PDSTMD,PDSTMN

! *** Common block for root development   kit 12/14/1999

	common/ rutcom/ SPDWRD,SPDWRN,dtop,ntop
   
	common/ datcom/ sum_fsq_tavg,sum_fbl_tavg,sum_fob_tavg,	&
	                ave_fsq_tavg,ave_fbl_tavg,ave_fob_tavg,	&                
					ifsqdae,ifbldae,ifobdae,ifbl,ifob,matday, &                
					sumxirr,sumirr,ipdays,irrflag,	&                
	                iniFertDate,iFDayAW

	common/ gcmcom/ igcmonth,radfactor(12),tmaxfactor(12),tminfactor(12),	&
					rainfactor(12),relhumfactor(12),windfactor(12),			&
					gcm(6),igcmflg
	
	end module common_block
