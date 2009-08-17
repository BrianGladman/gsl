;; Produces output for gsl_const header files using GNU Calc.
;;
;; Generate output with
;;
;;   emacs -batch -l const.el -f run 
;;

(setq calc-display-working-message t) ;; display short working messages
(setq calc-float-format '(sci 20))

(calc-eval "")
(load-library "calc/calc-units.el")
;;(calc-extensions)

(add-to-list 'load-path (expand-file-name "."))
(require 'calc-units-update)


(setq  gsl-dimensionless-constants
       '(("fsc"           "FINE_STRUCTURE")
         ("Nav"           "AVOGADRO")

         ("1e24"          "YOTTA")
         ("1e21"          "ZETTA")
         ("1e18"          "EXA")
         ("1e15"          "PETA")
         ("1e12"          "TERA")
         ("1e9"           "GIGA")
         ("1e6"           "MEGA")
         ("1e3"           "KILO")
         ("1e-3"          "MILLI")
         ("1e-6"          "MICRO")
         ("1e-9"          "NANO")
         ("1e-12"         "PICO")
         ("1e-15"         "FEMTO")
         ("1e-18"         "ATTO")
         ("1e-21"         "ZEPTO")
         ("1e-24"         "YOCTO")
         )
       )

(setq  gsl-constants
       '(("c"             "SPEED_OF_LIGHT")
         ("G"             "GRAVITATIONAL_CONSTANT")
         ("h"             "PLANCKS_CONSTANT_H")
         ("hbar"          "PLANCKS_CONSTANT_HBAR")

         ("au"            "ASTRONOMICAL_UNIT")
         ("float(lyr)"    "LIGHT_YEAR")
         ("pc"            "PARSEC")

         ("ga"            "GRAV_ACCEL")

         ("ev"            "ELECTRON_VOLT")
         ("me"            "MASS_ELECTRON")
         ("mmu"           "MASS_MUON")
         ("mp"            "MASS_PROTON")
         ("mn"            "MASS_NEUTRON")

         ("Ryd"           "RYDBERG")
         ("k"             "BOLTZMANN")
         ("R0"            "MOLAR_GAS")
         ("V0"            "STANDARD_GAS_VOLUME")

         ("min"           "MINUTE")
         ("hr"            "HOUR")
         ("day"           "DAY")
         ("wk"            "WEEK")

         ("in"            "INCH")
         ("ft"            "FOOT")
         ("yd"            "YARD")
         ("mi"            "MILE")
         ("nmi"           "NAUTICAL_MILE")
         ("fath"          "FATHOM")

         ("mil"           "MIL")
         ("point"         "POINT")
         ("texpt"         "TEXPOINT")

         ("mu"            "MICRON")
         ("Ang"           "ANGSTROM")

         ("hect"          "HECTARE")
         ("acre"          "ACRE")
         ("b"             "BARN")

         ("l"             "LITER")
         ("gal"           "US_GALLON")
         ("qt"            "QUART")
         ("pt"            "PINT")
         ("cup"           "CUP")
         ("ozfl"          "FLUID_OUNCE")
         ("tbsp"          "TABLESPOON")
         ("tsp"           "TEASPOON")
         ("galC"          "CANADIAN_GALLON")
         ("galUK"         "UK_GALLON")
         
         ("mph"           "MILES_PER_HOUR")
         ("kph"           "KILOMETERS_PER_HOUR")
         ("knot"          "KNOT")

         ("lb"            "POUND_MASS")
         ("oz"            "OUNCE_MASS")
         ("ton"           "TON")
         ("t"             "METRIC_TON")
         ("tonUK"         "UK_TON")
         ("ozt"           "TROY_OUNCE")
         ("ct"            "CARAT")
         ("amu"           "UNIFIED_ATOMIC_MASS")

         ("gf"            "GRAM_FORCE")
         ("lbf"           "POUND_FORCE")
         ("kip"           "KILOPOUND_FORCE")
         ("pdl"           "POUNDAL")

         ("cal"           "CALORIE")
         ("Btu"           "BTU")
         ("therm"         "THERM")

         ("hp"            "HORSEPOWER")
         
         ("bar"           "BAR")
         ("atm"           "STD_ATMOSPHERE")
         ("Torr"          "TORR")
         ("mHg"           "METER_OF_MERCURY")
         ("inHg"          "INCH_OF_MERCURY")
         ("inH2O"         "INCH_OF_WATER")
         ("psi"           "PSI")

         ("P"             "POISE")
         ("St"            "STOKES")
         
         ("sb"            "STILB")
         ("lm"            "LUMEN")
         ("lx"            "LUX")
         ("ph"            "PHOT")
         ("fc"            "FOOTCANDLE")
         ("lam"           "LAMBERT")
         ("flam"          "FOOTLAMBERT")
         
         ("Ci"            "CURIE")
         ("R"             "ROENTGEN")
         ("rd"            "RAD")

         ("1.98892e30 kg"       "SOLAR_MASS")
         ("0.5291772083e-10 m"  "BOHR_RADIUS")

         ("N"                     "NEWTON")
         ("1e-5 N"                "DYNE")
         ("J"                     "JOULE")
         ("1e-7 J"                "ERG")

         ("pi^2 k^4 / (60 hbar^3 c^2)"       "STEFAN_BOLTZMANN_CONSTANT")
         ("8 pi fsc^2 hbar^2/(3*c^2*me^2)" "THOMSON_CROSS_SECTION")

         )
       )

(setq gsl-em-constants 
      '(("muB"           "BOHR_MAGNETON")
        ("muN"           "NUCLEAR_MAGNETON")
        ("mue"           "ELECTRON_MAGNETIC_MOMENT")
        ("mup"           "PROTON_MAGNETIC_MOMENT")
        ("Fdy"           "FARADAY")
        ("e"             "ELECTRON_CHARGE")))

(setq gsl-special-em-constants 
      '(("8.854187817e-12 F/m" "VACUUM_PERMITTIVITY")
        ("mu0"           "VACUUM_PERMEABILITY")
        ("(1e-21/c) C/m" "DEBYE")
        ("Gs"            "GAUSS")))

;;; work around bug in calc 2.02f
(defun math-extract-units (expr)
  (if (memq (car-safe expr) '(* /))
      (cons (car expr)
	    (mapcar 'math-extract-units (cdr expr)))
    (if (math-units-in-expr-p expr nil) expr 1))
)

(defun fn (prefix system expr name)
  (let* ((x (calc-eval expr 'raw))
         (y (math-to-standard-units x system))
         (z (math-simplify-units y))
         (q (calc-eval (math-remove-units z)))
         (qq (format "evalv(%s + 0.0)" q))
         (quantity (calc-eval qq))
         (units (calc-eval (math-extract-units z)))
         )
    ;;(print x)
    ;;(print y)
    ;;(print z)
    ;;(print (math-extract-units z))
    ;;(print quantity)
    ;;(print units)
    (princ (format "#define %s_%s (%s) /* %s */\n" prefix name quantity units))
    )
  )

(setq cgsm (nth 1 (assq 'cgsm math-standard-units-systems)))
(setq mksa (nth 1 (assq 'mksa math-standard-units-systems)))

(setq cgs (nth 1 (assq 'cgs math-standard-units-systems)))
(setq mks (nth 1 (assq 'mks math-standard-units-systems)))

(defun display (prefix system constants)
  (princ (format "#ifndef __%s__\n" prefix))
  (princ (format "#define __%s__\n\n" prefix))
  (mapcar (lambda (x) (apply 'fn prefix system x)) constants)
  (princ (format "\n#endif /* __%s__ */\n" prefix))
)

(defun run-cgs ()
  (display "GSL_CONST_CGS" cgs gsl-constants)
)

(defun run-cgsm ()
  (display "GSL_CONST_CGSM" cgsm (append gsl-constants gsl-em-constants))
)

(defun run-mks ()
  (display "GSL_CONST_MKS" mks (append gsl-constants gsl-em-constants gsl-special-em-constants))
)

(defun run-mksa ()
  (display "GSL_CONST_MKSA" mksa (append gsl-constants gsl-em-constants gsl-special-em-constants))
)


(defun run-num ()
  (display "GSL_CONST_NUM" mks gsl-dimensionless-constants)
)
