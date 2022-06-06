from enum import IntEnum, IntFlag

class SPCFileType(IntFlag):
    """
    Describes the various file format flags.
    See the old specification for greater detail on what each setting changes.
    """
    DEFAULT = 0 # Will be 32 bit single file, one trace, with even x values
    SIXTEENPREC = 1	# 32 bit float or 16 bit spc format, currently only support 32 
    TCGRAM = 2	# Enables fexper in older software (not used accroding to old format doc)
    TMULTI = 4	# Multiple traces format (set if more than one subfile), Y vs YYYY, XY vs XYYY/XYXY
    TRANDM = 8	# If TMULTI and TRANDM=1 then arbitrary time (Z) values, must pass full z array
    TORDRD = 16	# If TMULTI and TORDRD=1 then ordered but uneven z, must pass full z array if true, else only z_inc and single z value in z array 
    TALABS = 32	# If true, use the specified custom axes
    TXYXYS = 64	# If true the file is XYXYXY, flase is considered XYYYYYY 
    TXVALS = 128	# Non even X values, must pass full x array 

# From GRAMSDDE.h
# Currently only support PPNONE
class SPCProcessCode(IntFlag):
    PPNONE	= 0    # No post processing 
    PPCOMP	= 1    # Compute (run PPCOMP?.ABP) 
    PPDLLC	= 2    # Compute with DLL (run PPCOMP?.DLL) 
    PPTRANS = 4    # Transmission (run PPTRANS?.ABP) 
    PPABSRB = 8    # Absorbance (run PPABSRB?.ABP) 
    PPKMUNK = 12  # Kuebelka-Munk (run PPKMUNK?.ABP) 
    PPPEAK	= 32   # GRAMS built-in peak picking and reporting 
    PPSRCH	= 64   # Library Search associated w/experiment's LIB driver 
    PPUSER	= 128  # General user-written post program (run PPUSER?.ABP) 

class SPCModFlags(IntFlag):
    UNMOD = 0 # unmodified
    A = 2**1 # Averaging (from multiple source traces)
    B = 2**2 # Baseline correction or offset functions
    C = 2**3 # Interferogram to spectrum Computation
    D = 2**4 # Derivative (or integrate) functions
    E = 2**6 # Resolution Enhancement functions (such as deconvolution)
    I = 2**9 # Interpolation functions
    N = 2**14 # Noise reduction smoothing
    O = 2**15 # Other functions (add, subtract, noise, etc.)
    S = 2**19 # Spectral Subtraction
    T = 2**20 # Truncation (only a portion of original X axis remains)
    W = 2**23 # When collected (date and time information) has been modified
    X = 2**24 # X units conversions or X shifting
    Y = 2**25 # Y units conversions (transmission->absorbance, etc.)
    Z = 2**26 # Zap functions (features removed or modified)

# looking at example parsers,
# axes z and w use the same units
class SPCXType(IntEnum):
    SPCXArb	        = 0
    SPCXWaven	    = 1
    SPCXUMetr	    = 2
    SPCXNMetr	    = 3
    SPCXSecs	    = 4
    SPCXMinuts	    = 5
    SPCXHertz	    = 6
    SPCXKHertz	    = 7
    SPCXMHertz	    = 8
    SPCXMUnits	    = 9
    SPCXPPM	        = 10
    SPCXDays	    = 11
    SPCXYears	    = 12
    SPCXRamans	    = 13
    SPCXeV	        = 14
    SPCZTextL	    = 15
    SPCXDiode	    = 16
    SPCXChanl	    = 17
    SPCXDegrs	    = 18
    SPCXDegrF	    = 19
    SPCXDegrC	    = 20
    SPCXDegrK	    = 21
    SPCXPoint	    = 22
    SPCXMSec	    = 23
    SPCXUSec	    = 24
    SPCXNSec	    = 25
    SPCXGHertz	    = 26
    SPCXCM	        = 27
    SPCXMeters	    = 28
    SPCXMMetr	    = 29
    SPCXHours	    = 30
    SPCXAngst	    = 31
    SPCXDblIgm	    = 255

class SPCTechType(IntEnum):
    SPCTechGen    = 0
    SPCTechGC     = 1
    SPCTechCgm    = 2
    SPCTechHPLC   = 3
    SPCTechFTIR   = 4
    SPCTechNIR    = 5
    SPCTechUV     = 7
    SPCTechXry    = 8
    SPCTechMS     = 9
    SPCTechNMR    = 10
    SPCTechRmn    = 11
    SPCTechFlr    = 12
    SPCTechAtm    = 13
    SPCTechDAD    = 14
    SPCTechThrm   = 15
    SPCTechCD     = 16
    SPCTechCNMR   = 20
    SPCTechHNMR   = 21
    SPCTechDNMR   = 22
    SPCTechANMR   = 23

class SPCYType(IntEnum):
    SPCYArb	        = 0
    SPCYIgram	    = 1
    SPCYAbsrb	    = 2
    SPCYKMonk	    = 3
    SPCYCount	    = 4
    SPCYVolts	    = 5
    SPCYDegrs	    = 6
    SPCYAmps	    = 7
    SPCYMeters      = 8
    SPCYMVolts      = 9
    SPCYLogdr	    = 10
    SPCYPercnt      = 11
    SPCYIntens      = 12
    SPCYRelInt      = 13
    SPCYEnergy      = 14
    SPCYDecbl	    = 16
    SPCYAbund	    = 17
    SPCYRelAbn      = 18
    SPCYDegrF	    = 19
    SPCYDegrC	    = 20
    SPCYDegrK	    = 21
    SPCYIdxRf	    = 22
    SPCYExtCf	    = 23
    SPCYReal	    = 24
    SPCYImag	    = 25
    SPCYCmplx	    = 26
    SPCYMgram	    = 27
    SPCYGram	    = 28
    SPCYKGram	    = 29
    SPCYSRot	    = 30
    SPCYTrans	    = 128
    SPCYReflec      = 129
    SPCYValley      = 130
    SPCYEmisn	    = 131

# From GRAMSDDE.h
# currently only supporting default PSTDEFT
class SPCPostDisposition(IntEnum):
    PSTDEFT = 0    # Use default setting (xfrpost is "unspecified") 
    PSTSAVE = 1    # Save file to disk (but remove from memory) 
    PSTAPPD = 2    # Append to end of #R database or xfrmtd .GDB database 
    PSTMERG = 3    # Merge into current #R database row 
    PSTBACK = 4    # Save as new background for experiment 
    PSTNONE = 5    # Do not save after processing (if any) 
    PSTKEEP = 6    # Do not save and keep in memory as #S 
    PSTBOTH = 7    # Both disk save & keep in #S memory (ABC Driver Only!) 
    PSTASK	= 8    # Query user: save, keep, or both (ABC Driver Only!) 

class SPCSubfileFlags(IntEnum):
    SUBNONE = 0 # No changes flag, this will likely be the only one used. The others seem to be for software that regularly change .spc files
    SUBCHGD = 1	# Subflags bit if subfile changed 
    SUBNOPT = 8	# Subflags bit if peak table file should not be used 
    SUBMODF = 128	# Subflags bit if subfile modified by arithmetic 